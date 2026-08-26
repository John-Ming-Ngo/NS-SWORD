# RNS

RNS is not included. Obtain it separately from:

https://github.com/rns-alberta/rns

It may be placed in this directory or kept elsewhere. /run_paper_grid.py uses it through `--rns`.

The EOS files must have the names and locations given in `examples/paper_grid.json`, relative to the RNS directory. The RNS output must contain the mass, equatorial radius, spin frequency, and a table headed by `Mu, Radius (KM)` or `Mu, Radii (KM)`. The parser also accepts the full-precision `RNS_EXACT` line used in our calculations.

Two RNS builds were used in the paper:

- `MDIV=261`, `SDIV=2081` for the radius, derivative, and differential-solid-angle data;
- `MDIV=261`, `SDIV=521` for the lower-resolution surface comparison.

Run the first build with:

```text
python scripts/run_paper_grid.py --rns PATH_TO_RNS --atmosphere PATH_TO_nsx_H_v200804.out --output paper_data_261x2081
```

For example, from the NS-SWORD repository root on Windows, with the executable and atmosphere table in their standard directories, run:

```text
python scripts/run_paper_grid.py --rns rns/rns.exe --atmosphere inputs/atmosphere/nsx_H_v200804.out --output paper_data_261x2081
```

On Linux or macOS, use `--rns rns/rns` instead.

Run the second build with:

```text
python scripts/run_paper_grid.py --rns PATH_TO_RNS --shapes-only --output paper_data_261x521
```

`model_manifest.csv` details each RNS model and its status. Valid surfaces and their exact stellar parameters are in the output directory.
