# Thymic mesenchymal niche cells drive T cell immune regeneration

Code and notebooks for *Gustafsson et al.* scRNA-Seq analysis of thymic mesenchymal cells.

## Content
- `ipynb_notebooks` contains .ipynb Jupyter notebooks:
  - `00_Data_preparation.ipynb`: QC and creation of AnnData objects for the further analysis,
  - `01_Murine_thymus.ipynb`: analysis of murine thymus composition,
  - `02_Human_thymus.ipynb`: analysis of human thymus composition,
  - `03_Cross_species_alignment.ipynb`: comparative analysis of murine and human thymic mesenchymal cells and validation of MC cell types on public data,
  - `04_Velocity_analysis.ipynb`: RNA velocity analysis of Postn+ and Penk+ MCs from murine thymus,
  - `05_Bornstein_alignment.ipynb`: validation of MC cell types on public data,
  - `06_Beta_regression.ipynb`: some statistical computations,
  - `07_Figures.ipynb`: figures preparation;
- `html_notebooks` contains .html Jupyter notebooks:
  - *same as in `ipynb_notebooks`*;
- `data` contains some data that is needed to reproduce an analysis;
- `tools` contains additional functions that are used in Jupyter notebooks.
