# NS-SWORD

NS-SWORD calculates the surface geometry and radiation of a rotating neutron star in the oblate Schwarzschild approximation. A surface may be supplied by RNS or by one of the included shape functions.

## Requirements

- GNU Make;
- a C++17 compiler;
- Python 3;
- the hydrogen atmosphere table described in `inputs/atmosphere/README.md`.

NumPy and PyTorch are only required for the PVLS shape function. They may be installed with:

```text
python -m pip install -r requirements.txt
```

## Build

Run:

```text
make
```

This produces `NS-SWORD` on Linux and macOS, or `NS-SWORD.exe` on Windows.

## Running NS-SWORD

For example, this runs a 1.4-solar-mass, 12-km star at 600 Hz using the Morsink `oblate` shape function:

```text
.\NS-SWORD.exe -Smn oblate -Sm 1.4 -Sf 600 -Sr 12 -Ss Morsink -Sq 0 -Oi 45 -Gm 128 -Od 6.1713552e18 -Gn 128 -St 0.0861733326 -Fr outputs/morsink_model.csv -Fs outputs/morsink_spectrum_ -Fg outputs/morsink_spots.csv -Fsh outputs/morsink_shape.csv
```

On Linux or macOS, use `./NS-SWORD` instead. `-Sm`, `-Sr`, `-Sf`, and `-St` set the mass in solar masses, equatorial radius in km, spin in Hz, and temperature in keV. `-Oi` and `-Od` set the observer inclination in degrees and distance in metres. `-Gm` and `-Gn` set the number of latitude and longitude bins.

The output flags are:

- `-Fr`: integrated model results;
- `-Fs`: spectra results (in beta), producing `morsink_spectrum_BB.csv` and `morsink_spectrum_H.csv`;
- `-Fg`: the surface spot grid;
- `-Fsh`: the sampled surface shape.

Without `-Fr`, the integrated model results are printed to the terminal. The hydrogen atmosphere table must be in `inputs/atmosphere/nsx_H_v200804.out`.

## Reproducing the paper data

The inputs from table 3 in the paper are in `examples/paper_grid.json`. This gives 882 RNS models. Only valid stellar models are passed to NS-SWORD.

RNS is not included. Prepare it as described in `rns/README.md`. For the surface, derivative, and differential-solid-angle data, build RNS with `MDIV=261` and `SDIV=2081`, then run:

```text
python scripts/run_paper_grid.py --rns PATH_TO_RNS --atmosphere PATH_TO_nsx_H_v200804.out --output paper_data_261x2081
```

/run_paper_grid.py runs rns and NS-SWORD with all of those inputs. For example, when running from the repository root with `rns.exe` in `rns/` and the atmosphere table in its standard input directory, run:

```text
python scripts/run_paper_grid.py --rns rns/rns.exe --atmosphere inputs/atmosphere/nsx_H_v200804.out --output paper_data_261x2081
```

On Linux or macOS, use `--rns rns/rns` instead.

The command produces:

- `model_manifest.csv`, detailing each RNS model and its status;
- `rns_surfaces/`, containing each outputted tabulated surface;
- `shape_outputs/`, containing the RNS, MLCB, BBPO, AM, SPYY(s), SPYY(f), and PVLS surfaces;
- `spot_grids/`, containing the same models at inclinations of 0, 45, and 90 degrees;
- `model_outputs/`, containing the integrated quantities for each spot grid;
- `RunParams_shapes.csv` and `RunParams_spots.csv`, containing every NS-SWORD input;
- `run_errors.csv`, containing any failed NS-SWORD calculation.

For the lower-resolution surface comparison, rebuild RNS with `MDIV=261` and `SDIV=521`, then run:

```text
python scripts/run_paper_grid.py --rns PATH_TO_RNS --shapes-only --output paper_data_261x521
```

Use `--plan` to check the number of calculations without running them. Plotting code is not included.

## Repository contents

- `src/` and `include/`: the main calculation;
- `shape_functions/`: the included surface models;
- `inputs/`: coefficient files and locations for external data;
- `examples/paper_grid.json`: the full input grid;
- `scripts/run_paper_grid.py`: a complete RNS-to-NS-SWORD calculation;
- `rns/`: the location for a separate RNS copy.

Read `NOTICE.md` before distributing this repository.
