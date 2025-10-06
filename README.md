# DM recast limits

This repository contains the Jupyter notebooks and input data used to reproduce part of the analysis presented in:

> *G. D’Amico, M. Doro, and M. De Caria, “Recasting and Forecasting Dark Matter Limits Without Raw Data: A Generalized Algorithm for Gamma-Ray Telescopes” (submitted to Physics of the Dark Universe, 2025).*

[**Browse the API Docs**](https://giacomodamico24.github.io/DM_recast_limits/html/)

---

## 📘 Description

The notebooks in this repository demonstrate the use of the **recasting framework** introduced in the paper “Recasting and Forecasting Dark Matter Limits Without Raw Data: A Generalized Algorithm for Gamma-Ray Telescopes”.
They allow users to reproduce **Figures 5 and 6** from the manuscript, corresponding to the validation of the recasting method using published MAGIC, Fermi-LAT, LHAASO upper limits o DM and the projected CTAO upper limits from the Galactic Center.

Each subfolder corresponds to a specific gamma-ray instrument:
- **`RECAST_MAGIC`** – Recasting MAGIC Collaboration limits for the Coma Berenices dSph.  
- **`RECAST_CTAO`** – Recasting CTAO projected limits for Galactic Center.  
- **`RECAST_FERMI`** – Example of recasting based on Fermi-LAT published upper limits.  
- **`RECAST_LHAASO`** – Example adaptation for LHAASO data.

---


## 🛠️ Dependencies

All notebooks require the following Python packages:

```bash
numpy (preferably version 1.26)
scipy (preferably version 5.11)
matplotlib (preferably version 3.8)
astropy (preferably version 5.1)
gammapy version 1.2


