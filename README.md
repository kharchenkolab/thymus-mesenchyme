# Code for the paper "Thymic mesenchymal niche cells drive T cell immune regeneration"

Code for *Gustafsson et al.* scRNA-Seq analysis of thymic mesenchymal cells. AnnData and conos objects are available from Zenodo repository (https://doi.org/10.5281/zenodo.6656930).

## Content
- `ipynb_notebooks` contains .ipynb Jupyter notebooks:
  - `00_Data_preparation.ipynb`: QC and creation of AnnData objects for the further analysis ([view](https://nbviewer.org/github/kharchenkolab/thymus-mesenchyme/blob/main/ipynb_notebooks/00_Data_preparation.ipynb)),
  - `01_Murine_thymus.ipynb`: analysis of murine thymus composition ([view](https://nbviewer.org/github/kharchenkolab/thymus-mesenchyme/blob/main/ipynb_notebooks/01_Murine_thymus.ipynb)),
  - `02_Human_thymus.ipynb`: analysis of human thymus composition ([view](https://nbviewer.org/github/kharchenkolab/thymus-mesenchyme/blob/main/ipynb_notebooks/02_Human_thymus.ipynb)),
  - `03_Cross_species_alignment.ipynb`: comparative analysis of murine and human thymic mesenchymal cells and validation of MC cell types on public data ([view](https://nbviewer.org/github/kharchenkolab/thymus-mesenchyme/blob/main/ipynb_notebooks/03_Cross_species_alignment.ipynb)),
  - `04_Velocity_analysis.ipynb`: RNA velocity analysis of Postn+ and Penk+ MCs from murine thymus ([view](https://nbviewer.org/github/kharchenkolab/thymus-mesenchyme/blob/main/ipynb_notebooks/04_Velocity_analysis.ipynb)),
  - `05_Bornstein_alignment.ipynb`: validation of MC cell types on public data ([view](https://nbviewer.org/github/kharchenkolab/thymus-mesenchyme/blob/main/ipynb_notebooks/05_Bornstein_alignment.ipynb)),
  - `06_Beta_regression.ipynb`: some statistical computations ([view](https://nbviewer.org/github/kharchenkolab/thymus-mesenchyme/blob/main/ipynb_notebooks/06_Beta_regression.ipynb)),
  - `07_Figures.ipynb`: figures preparation ([view](https://nbviewer.org/github/kharchenkolab/thymus-mesenchyme/blob/main/ipynb_notebooks/07_Figures.ipynb));
- `tools` contains additional functions that are used in Jupyter notebooks.
