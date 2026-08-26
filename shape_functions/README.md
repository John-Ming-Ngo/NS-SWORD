# Shape functions

This directory contains the shape functions as implemented for use by NS-SWORD. Running `make` builds each model as a shared library in `shape_functions/lib/`.

The included models are Morsink et al., AlGendy and Morsink, Baubock et al., Silva et al., PVLS, and a manually tabulated surface. PVLS requires NumPy, PyTorch, and the included trained weights.

The origin of the PVLS Python files is recorded in `NOTICE.md` and `dependencies/Papigkiotis/README.md`.
