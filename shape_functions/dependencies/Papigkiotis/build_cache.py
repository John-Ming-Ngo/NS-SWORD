"""Build strict, file-backed PVLS inference tables.

The C++ model invokes the Python entry points with default iostream precision
(six significant digits).  This builder deliberately performs inference from
those rendered values so a cache hit is numerically identical to the existing
Python fallback.  Out-of-domain cases are never written.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
CANONICAL_DIR = Path(__file__).resolve().parent
LEGACY_DIR = CANONICAL_DIR.parent / "Papigkiotis_old"
DEFAULT_CACHE_DIR = REPO_ROOT / ".pvls_migration" / "cache" / "PVLS"
MAGIC = "PVLS_INFERENCE_CACHE_V1"

C_MIN, C_MAX = 0.0876346858172578, 0.3094541325480277
SIGMA_MIN, SIGMA_MAX = 0.0, 0.9612274013913829
ECCENTRICITY_MIN, ECCENTRICITY_MAX = 0.0, 0.7797886226038347
RP_REQ_MIN, RP_REQ_MAX = 0.6260428931452016, 1.0

C_LIGHT = 2.99792458e10
G_GRAV = 6.6742e-8
M_SUN = 1.9891e33
PI = 3.14159265358979323846
KAPPA = 1.0e-15 * C_LIGHT * C_LIGHT / G_GRAV

DEPENDENCIES = (
    "shape_functions/dependencies/Papigkiotis/ns_radius.py",
    "shape_functions/dependencies/Papigkiotis/ns_log_derivative.py",
    "shape_functions/dependencies/Papigkiotis_old/ns_radius.py",
    "shape_functions/dependencies/Papigkiotis_old/ns_log_derivative.py",
    "shape_functions/dependencies/Papigkiotis/DNN.py",
    "shape_functions/dependencies/Papigkiotis/Model/Surface/Surface-model.pth",
    "shape_functions/dependencies/Papigkiotis/Model/Derivative/Derivative-model.pth",
)


@dataclass(frozen=True, order=True)
class Case:
    req: str
    compactness: str
    sigma: str

    @property
    def key(self) -> str:
        return f"{self.req}\t{self.compactness}\t{self.sigma}"


def cxx_default(value: float) -> str:
    """Match defaultfloat with the C++ default precision of six."""
    return format(value, ".6g")


def polar_ratio(compactness: float, sigma: float) -> float:
    coefficients = {
        (0, 0): 0.942328,
        (0, 1): -0.617711,
        (0, 2): 0.544639,
        (0, 3): -0.440968,
        (0, 4): 0.196118,
        (1, 0): 1.296632,
        (1, 1): -1.458921,
        (1, 2): -0.226904,
        (1, 3): 0.527775,
        (2, 0): -10.45611,
        (2, 1): 8.668382,
        (2, 2): -2.506686,
        (3, 0): 36.131881,
        (3, 1): -7.524662,
        (4, 0): -45.301523,
    }
    return sum(
        coefficients.get((n, m), 0.0) * compactness**n * sigma**m
        for n in range(5)
        for m in range(5 - n)
    )


def domain_values(compactness: float, sigma: float) -> tuple[bool, float, float]:
    ratio = polar_ratio(compactness, sigma)
    eccentricity_squared = 1.0 - ratio * ratio
    eccentricity = math.sqrt(eccentricity_squared) if eccentricity_squared >= 0.0 else math.nan
    valid = (
        math.isfinite(compactness)
        and math.isfinite(sigma)
        and math.isfinite(ratio)
        and math.isfinite(eccentricity)
        and C_MIN <= compactness <= C_MAX
        and SIGMA_MIN <= sigma <= SIGMA_MAX
        and ECCENTRICITY_MIN <= eccentricity <= ECCENTRICITY_MAX
        and RP_REQ_MIN <= ratio <= RP_REQ_MAX
    )
    return valid, ratio, eccentricity


def case_from_values(req: float, compactness: float, sigma: float) -> Case:
    if not domain_values(compactness, sigma)[0]:
        raise ValueError(f"out of trained domain: C={compactness}, sigma={sigma}")
    case = Case(cxx_default(req), cxx_default(compactness), cxx_default(sigma))
    rendered_valid, ratio, eccentricity = domain_values(
        float(case.compactness), float(case.sigma)
    )
    if not rendered_valid:
        raise ValueError(
            "six-significant-digit inference key is out of trained domain: "
            f"{case.key}; Rp/Req={ratio}; eccentricity={eccentricity}"
        )
    return case


def physical_case(mass_msun: float, frequency_hz: float, radius_km: float) -> Case:
    length_factor = math.sqrt(KAPPA)
    mass_factor = math.sqrt(KAPPA) * (G_GRAV**-1.0) * C_LIGHT**2.0
    inverse_time_factor = KAPPA**-0.5 * C_LIGHT
    req = radius_km * 1.0e5 / length_factor
    mass = mass_msun * M_SUN / mass_factor
    omega = 2.0 * PI * frequency_hz / inverse_time_factor
    compactness = mass / req
    sigma = omega**2.0 * req**3.0 / mass
    return case_from_values(req, compactness, sigma)


def fnv1a(data: bytes, initial: int = 14695981039346656037) -> int:
    value = initial
    for byte in data:
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value


def key_hash(case: Case) -> str:
    return f"{fnv1a(case.key.encode('utf-8')):016x}"


def dependency_fingerprint() -> str:
    value = 14695981039346656037
    for relative in DEPENDENCIES:
        value = fnv1a(relative.encode("utf-8"), value)
        value = fnv1a(b"\0", value)
        with (REPO_ROOT / relative).open("rb") as source:
            while chunk := source.read(8192):
                value = fnv1a(chunk, value)
    return f"{value:016x}"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def iter_run_params(paths: Iterable[Path]) -> tuple[set[Case], int]:
    cases: set[Case] = set()
    excluded = 0
    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                if row.get("-Ss") != "Papigkiotis":
                    continue
                try:
                    cases.add(
                        physical_case(
                            float(row["-Sm"]), float(row["-Sf"]), float(row["-Sr"])
                        )
                    )
                except ValueError:
                    excluded += 1
    return cases, excluded


def iter_key_file(path: Path) -> set[Case]:
    cases: set[Case] = set()
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            req = row.get("req") or row.get("req_nounits")
            compactness = row.get("compactness") or row.get("C")
            sigma = row.get("sigma")
            if req is None or compactness is None or sigma is None:
                raise ValueError(f"{path} must contain req, compactness, and sigma columns")
            cases.add(case_from_values(float(req), float(compactness), float(sigma)))
    return cases


def cache_text(case: Case, radius_module, derivative_module, fingerprint: str) -> str:
    req = float(case.req)
    compactness = float(case.compactness)
    sigma = float(case.sigma)

    radius_pole = radius_module.legacy.R_pole(req, compactness, sigma)
    surface_mu, surface_values = radius_module.r_mu(
        compactness, sigma, radius_pole, req
    )
    derivative_pole = derivative_module.legacy.R_pole(req, compactness, sigma)
    derivative_mu, derivative_values = derivative_module.dlogR_dmu(
        compactness, sigma, derivative_pole, req
    )
    if len(surface_mu) != 20 or len(derivative_mu) != 500:
        raise RuntimeError("PVLS inference entry points changed table dimensions")

    lines = [MAGIC, case.key, "trained-domain", f"fingerprint\t{fingerprint}", "20 500"]
    lines.extend(f"S {mu} {value}" for mu, value in zip(surface_mu, surface_values))
    lines.extend(f"D {mu} {value}" for mu, value in zip(derivative_mu, derivative_values))
    return "\n".join(lines) + "\n"


def write_atomic(path: Path, content: str) -> bool:
    encoded = content.encode("ascii")
    if path.exists() and path.read_bytes() == encoded:
        return False
    temporary = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    temporary.write_bytes(encoded)
    os.replace(temporary, path)
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-params", type=Path, nargs="*", default=[])
    parser.add_argument("--keys-file", type=Path, action="append", default=[])
    parser.add_argument(
        "--key",
        nargs=3,
        type=float,
        action="append",
        metavar=("REQ", "COMPACTNESS", "SIGMA"),
        default=[],
    )
    parser.add_argument(
        "--physical-key",
        nargs=3,
        type=float,
        action="append",
        metavar=("MASS_MSUN", "FREQUENCY_HZ", "RADIUS_KM"),
        default=[],
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="classify and deduplicate keys without loading PyTorch or writing cache files",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cases, excluded = iter_run_params(args.run_params)
    for path in args.keys_file:
        cases.update(iter_key_file(path))
    for req, compactness, sigma in args.key:
        cases.add(case_from_values(req, compactness, sigma))
    for mass_msun, frequency_hz, radius_km in args.physical_key:
        cases.add(physical_case(mass_msun, frequency_hz, radius_km))
    if not cases:
        raise SystemExit("no in-domain PVLS inference keys were supplied")
    if args.dry_run:
        print(f"PVLS cache dry run: keys={len(cases)}, excluded_rows={excluded}")
        return 0

    os.chdir(REPO_ROOT)
    sys.path.insert(0, str(LEGACY_DIR))
    radius_module = load_module("pvls_cache_radius", CANONICAL_DIR / "ns_radius.py")
    derivative_module = load_module(
        "pvls_cache_derivative", CANONICAL_DIR / "ns_log_derivative.py"
    )
    fingerprint = dependency_fingerprint()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    written = 0
    for index, case in enumerate(sorted(cases), start=1):
        destination = args.output_dir / f"{key_hash(case)}.pvls"
        written += int(
            write_atomic(
                destination,
                cache_text(case, radius_module, derivative_module, fingerprint),
            )
        )
        if index % 50 == 0 or index == len(cases):
            print(f"PVLS cache: {index}/{len(cases)} keys complete", flush=True)
    print(
        f"PVLS cache ready: keys={len(cases)}, written={written}, "
        f"unchanged={len(cases) - written}, excluded_rows={excluded}, "
        f"fingerprint={fingerprint}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
