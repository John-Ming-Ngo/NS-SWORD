"""Generate the complete raw data set for the parameter grid in Table III."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "examples" / "paper_grid.json"
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

MODELS = (
    ("RNS", "Tabulated", ""),
    ("MLCB", "Morsink", "oblate"),
    ("BBPO", "BaubockEtAlBL", "Whatever"),
    ("AM", "AlGendy", "Original"),
    ("SPYY_s", "SilvaEtAl", "slow-eliptical-fit"),
    ("SPYY_f", "SilvaEtAl", "fast-eliptical-fit"),
    ("PVLS", "Papigkiotis", "Whatever"),
)

PARAM_HEADER = ["-Smn", "-Sm", "-Sf", "-Sr", "-Ss", "-Sq", "-Oi", "-Gm", "-Od", "-Gn", "-St"]


@dataclass(frozen=True)
class Case:
    eos: str
    eos_file: Path
    density_axis: str
    density_cgs: str
    input_spin: str
    name: str


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--output", default="paper_data_261x2081")
    parser.add_argument("--rns", help="Path to the RNS executable")
    parser.add_argument("--atmosphere", help="Path to nsx_H_v200804.out")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--shapes-only", action="store_true")
    parser.add_argument("--plan", action="store_true")
    return parser.parse_args()


def absolute(path_text: str, base: Path = ROOT) -> Path:
    path = Path(path_text)
    return (path if path.is_absolute() else base / path).resolve()


def find_executable(configured: str | None, names: tuple[str, ...], description: str) -> Path:
    candidates = (absolute(configured),) if configured else tuple(ROOT / name for name in names)
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(f"{description} was not found. See rns/README.md.")


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def axis_text(value: object) -> str:
    return format(float(value), ".17g")


def slug(value: str) -> str:
    return value.replace("+", "").replace("-", "m").replace(".", "p")


def cases(config: dict, rns_root: Path) -> list[Case]:
    result = []
    for (eos, eos_path), density, spin in itertools.product(
        config["eos"].items(), config["central_density_1e14_g_cm3"], config["spin_hz"]
    ):
        density_axis = axis_text(density)
        density_cgs = format(float(density) * 1.0e14, ".17g")
        spin_text = axis_text(spin)
        name = f"rho{slug(density_axis)}e14_spin{slug(spin_text)}"
        result.append(Case(eos, (rns_root / eos_path).resolve(), density_axis, density_cgs, spin_text, name))
    return result


def parse_rns_output(text: str) -> tuple[str, str, str, list[tuple[str, str]]]:
    exact = re.search(r"^RNS_EXACT,(.+)$", text, re.MULTILINE)
    if exact:
        fields = exact.group(1).split(",")
        if any("=" not in field for field in fields):
            raise ValueError("RNS_EXACT line is malformed")
        values = dict(field.split("=", 1) for field in fields)
        missing = [name for name in ("mass_msun", "radius_km", "spin_hz") if name not in values]
        if missing:
            raise ValueError(f"RNS_EXACT line is missing: {', '.join(missing)}")
        mass, radius, spin = values["mass_msun"], values["radius_km"], values["spin_hz"]
    else:
        summary = re.search(r"e15\s+Msun\s+km\s+Hz\s*\n\s*\S+\s+(\S+)\s+(\S+)\s+(\S+)", text)
        if not summary:
            raise ValueError("mass, radius, and spin were not found")
        mass, radius, spin = summary.groups()
    header = re.search(r"^\s*Mu\s*,\s*(?:Radius|Radii).*?$", text, re.MULTILINE | re.IGNORECASE)
    if not header:
        raise ValueError("surface table was not found")
    rows = []
    row_pattern = re.compile(r"^\s*(" + NUMBER + r")\s*,\s*(" + NUMBER + r")\s*$")
    for line in text[header.end():].splitlines():
        if not line.strip():
            continue
        match = row_pattern.match(line)
        if not match:
            break
        rows.append(match.groups())
    if not rows:
        raise ValueError("surface table contained no rows")
    return mass, radius, spin, rows


def write_csv(path: Path, header: list[str] | tuple[str, ...], rows: list[list[str]] | list[tuple[str, ...]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    temporary.replace(path)


def save_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "eos", "central_density_g_cm3", "input_spin_hz", "status", "mass_msun",
        "equatorial_radius_km", "output_spin_hz", "surface_file", "error",
    ]
    write_csv(path, fields, [[row.get(field, "") for field in fields] for row in rows])


def run_rns(case: Case, executable: Path, output: Path, timeout: int) -> dict[str, str]:
    row = {
        "eos": case.eos, "central_density_g_cm3": case.density_cgs, "input_spin_hz": case.input_spin,
        "status": "failed", "mass_msun": "", "equatorial_radius_km": "", "output_spin_hz": "",
        "surface_file": "", "error": "",
    }
    try:
        if not case.eos_file.is_file():
            raise FileNotFoundError(f"EOS file not found: {case.eos_file}")
        eos_argument = os.path.relpath(case.eos_file, executable.parent)
        process = subprocess.run(
            [str(executable), "-f", eos_argument, "-e", case.density_cgs, "-s", case.input_spin],
            cwd=executable.parent, capture_output=True, text=True, timeout=timeout,
        )
        if process.returncode != 0:
            message = process.stderr.strip() or process.stdout.strip()
            if message:
                raise RuntimeError(f"RNS failed: {message.splitlines()[-1]}")
            raise RuntimeError(f"RNS returned exit code {process.returncode}")
        mass, radius, spin, surface = parse_rns_output(process.stdout)
        surface_path = output / "rns_surfaces" / case.eos / f"{case.name}.csv"
        write_csv(surface_path, ("Mu", "Radius (KM)"), surface)
        row.update({
            "status": "valid", "mass_msun": mass, "equatorial_radius_km": radius,
            "output_spin_hz": spin, "surface_file": surface_path.relative_to(output).as_posix(),
        })
    except subprocess.TimeoutExpired:
        row["error"] = f"RNS timed out after {timeout} seconds"
    except Exception as exc:
        row["error"] = str(exc).replace("\r", " ").replace("\n", " ")[:1000]
    return row


def common_values(config: dict, model: dict[str, str], library: str, model_name: str, inclination: object) -> list[str]:
    return [
        model_name, model["mass_msun"], model["output_spin_hz"], model["equatorial_radius_km"],
        library, "0", axis_text(inclination), str(config["latitude_bins"]),
        axis_text(config["observer_distance_m"]), str(config["longitude_bins"]), axis_text(config["temperature_kev"]),
    ]


def program_tasks(config: dict, output: Path, valid: list[dict[str, str]], shapes_only: bool) -> tuple[list[tuple[list[str], tuple[Path, ...]]], list[list[str]], list[list[str]]]:
    tasks: list[tuple[list[str], tuple[Path, ...]]] = []
    shape_params: list[list[str]] = []
    spot_params: list[list[str]] = []
    for star in valid:
        case_name = Path(star["surface_file"]).stem
        tabulated_name = f"{star['eos']}/{case_name}"
        for label, library, fixed_name in MODELS:
            model_name = tabulated_name if library == "Tabulated" else fixed_name
            shape_path = output / "shape_outputs" / star["eos"] / case_name / f"{label}.csv"
            values = common_values(config, star, library, model_name, 0)
            shape_params.append(values + [shape_path.relative_to(output).as_posix()])
            tasks.append((values + ["-Fsh", str(shape_path)], (shape_path,)))
            if shapes_only:
                continue
            for inclination in config["inclination_deg"]:
                spot_path = output / "spot_grids" / star["eos"] / case_name / label / f"i{axis_text(inclination)}.csv"
                model_path = output / "model_outputs" / star["eos"] / case_name / label / f"i{axis_text(inclination)}.csv"
                spot_values = common_values(config, star, library, model_name, inclination)
                spot_params.append(spot_values + [model_path.relative_to(output).as_posix(), spot_path.relative_to(output).as_posix()])
                tasks.append((spot_values + ["-Fr", str(model_path), "-Fg", str(spot_path)], (model_path, spot_path)))
    return tasks, shape_params, spot_params


def run_program(command_values: list[str], expected: tuple[Path, ...], executable: Path, environment: dict[str, str]) -> tuple[str, str]:
    for path in expected:
        path.parent.mkdir(parents=True, exist_ok=True)
    command = [str(executable)]
    for flag, value in zip(PARAM_HEADER, command_values[:len(PARAM_HEADER)]):
        command.extend((flag, value))
    command.extend(command_values[len(PARAM_HEADER):])
    process = subprocess.run(command, cwd=ROOT, env=environment, capture_output=True, text=True)
    if process.returncode != 0 or not all(path.is_file() and path.stat().st_size > 0 for path in expected):
        message = process.stderr.strip() or process.stdout.strip() or f"NS-SWORD returned {process.returncode}"
        return "failed", message.replace("\r", " ").replace("\n", " ")[:1000]
    return "ok", ""


def main() -> int:
    args = arguments()
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    config_path = absolute(args.config)
    config = load_config(config_path)
    nominal = len(config["eos"]) * len(config["central_density_1e14_g_cm3"]) * len(config["spin_hz"])
    print(f"RNS models: {nominal}")
    print(f"shape files, at most: {nominal * len(MODELS)}")
    print(f"spot-grid files, at most: {0 if args.shapes_only else nominal * len(MODELS) * len(config['inclination_deg'])}")
    if args.plan:
        return 0

    output = absolute(args.output)
    if output.exists():
        raise FileExistsError(f"output folder already exists: {output}")
    output.mkdir(parents=True, exist_ok=True)
    rns_exe = find_executable(args.rns, ("rns/rns.exe", "rns/rns"), "RNS")
    rns_root = rns_exe.parent
    all_cases = cases(config, rns_root)
    manifest_path = output / "model_manifest.csv"
    manifest: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(run_rns, case, rns_exe, output, int(config["rns_timeout_seconds"])): case
            for case in all_cases
        }
        for index, future in enumerate(as_completed(futures), 1):
            manifest.append(future.result())
            print(f"RNS {index}/{len(all_cases)}", flush=True)
    manifest.sort(key=lambda row: (list(config["eos"]).index(row["eos"]), float(row["central_density_g_cm3"]), float(row["input_spin_hz"])))
    save_manifest(manifest_path, manifest)
    valid = [row for row in manifest if row["status"] == "valid"]
    print(f"valid RNS models: {len(valid)}")
    if not valid:
        raise RuntimeError("RNS produced no valid models; inspect model_manifest.csv")

    ns_sword = find_executable(None, ("NS-SWORD.exe", "NS-SWORD"), "NS-SWORD")
    atmosphere = absolute(args.atmosphere) if args.atmosphere else ROOT / "inputs/atmosphere/nsx_H_v200804.out"
    if not args.shapes_only and not atmosphere.is_file():
        raise FileNotFoundError("the hydrogen atmosphere table was not found; see inputs/atmosphere/README.md")
    environment = os.environ.copy()
    environment["NS_SWORD_RNS_DIR"] = str((output / "rns_surfaces").resolve())
    if atmosphere.is_file():
        environment["NS_SWORD_ATMOSPHERE"] = str(atmosphere.resolve())

    tasks, shape_params, spot_params = program_tasks(config, output, valid, args.shapes_only)
    write_csv(output / "RunParams_shapes.csv", PARAM_HEADER + ["-Fsh"], shape_params)
    if not args.shapes_only:
        write_csv(output / "RunParams_spots.csv", PARAM_HEADER + ["-Fr", "-Fg"], spot_params)
    errors: list[list[str]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_program, values, paths, ns_sword, environment): paths for values, paths in tasks}
        for index, future in enumerate(as_completed(futures), 1):
            status, error = future.result()
            if status == "failed":
                errors.append([";".join(path.relative_to(output).as_posix() for path in futures[future]), error])
            if index % 100 == 0 or index == len(tasks):
                print(f"NS-SWORD {index}/{len(tasks)}", flush=True)
    write_csv(output / "run_errors.csv", ("output_file", "error"), errors)
    print(f"errors: {len(errors)}")
    print(output)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
