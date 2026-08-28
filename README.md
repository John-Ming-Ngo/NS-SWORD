# NS-SWORD

NS-SWORD (Neutron Star - Schwarzschild With Oblate Rotational Deformations) calculates the surface geometry and radiation of a rotating neutron star in the oblate Schwarzschild approximation. A surface may be supplied by RNS or by one of the included shape functions.

This code was used in "The Impact of Errors in the Shape Function of Rotating Neutron Stars in the Oblate Schwarzschild Approximation," by John Ming Ngo, Charlee Amason, and Sharon M. Morsink. 

The solid angle calculations are the rigorously tested parts of NS-SWORD. Flux-related calculations, such as spectra, are in beta and have not been rigorously tested.

## Requirements

- GNU Make;
- a C++17 compiler;
- Python 3;

NumPy and PyTorch are only required for the PVLS shape function. They may be installed with:

```text
python -m pip install -r requirements.txt
```

## Optional

- the hydrogen atmosphere table described in `inputs/atmosphere/README.md`. If present, it enables the hydrogen spectrum (spectra generation is in beta).

## Build

Run:

```text
make
```

This produces `NS-SWORD` on Linux and macOS, or `NS-SWORD.exe` on Windows.

## Running NS-SWORD

For example, this runs a 1.4-solar-mass, 12-km star at 600 Hz using the Morsink `oblate` shape function:

```text
.\NS-SWORD.exe -Smn oblate -Sm 1.4 -Sf 600 -Sr 12 -Ss Morsink -Oi 45 -Gm 128 -Od 6.1713552e18 -Gn 128 -St 0.0861733326 -Fr outputs/morsink_model.csv -Fg outputs/morsink_spots.csv -Fsh outputs/morsink_shape.csv
```

On Linux or macOS, use `./NS-SWORD` instead. All the following flags are for inputs and are mandatory. `-Sm`, `-Sr`, `-Sf`, and `-St` set the mass in solar masses, equatorial radius in km, spin in Hz, and temperature in keV. `-Ss` and `-Smn` set the shape library and model name. `-Oi` and `-Od` set the observer inclination in degrees and distance in metres. `-Gm` and `-Gn` set the number of latitude and longitude bins.

The optional output flags are:

- `-Fr` (optional): Produce the model summary statistics;
- `-Fg` (optional): Produce the data for each point on the surface of the star, 2d data;
- `-Fsh` (optional): Produce the radial shape of the model, 1d data;
- `-Fs` (optional): Produce the flux spectra results (in beta), producing `morsink_spectrum_BB.csv` and `morsink_spectrum_H.csv`. If the hydrogen atmosphere table is not available, the latter file is not produced.

Without `-Fr`, the integrated model results are printed to the terminal. 

## Reproducing the paper data

The inputs from Table III in the paper are in `examples/paper_grid.json`. This gives 882 RNS models. Only valid stellar models are passed to NS-SWORD.

RNS is not included. Prepare it as described in `rns/README.md`. For the surface, derivative, and differential-solid-angle data, build RNS with `MDIV=261` and `SDIV=2081`, then run:

```text
python scripts/run_paper_grid.py --rns PATH_TO_RNS --output paper_data_261x2081
```

`scripts/run_paper_grid.py` runs RNS and NS-SWORD with all of those inputs. For example, when running from the repository root with `rns.exe` in `rns/`, run:

```text
python scripts/run_paper_grid.py --rns rns/rns.exe --output paper_data_261x2081
```

On Linux or macOS, use `--rns rns/rns` instead.

The command produces:

- `model_manifest.csv`, detailing each RNS model and its status;
- `rns_surfaces/`, containing each tabulated surface produced by RNS;
- `shape_outputs/`, containing the surfaces calculated by RNS, MLCB, BBPO, AM, SPYY(s), SPYY(f), and PVLS;
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

- `src/` and `include/`: the main code;
- `shape_functions/`: the implemented shape functions;
- `inputs/`: coefficient files and locations for external data;
- `examples/paper_grid.json`: the full input grid;
- `scripts/run_paper_grid.py`: a complete RNS-to-NS-SWORD calculation;
- `rns/`: the location for a separate RNS copy.

Read `NOTICE.md` before distributing this repository.
