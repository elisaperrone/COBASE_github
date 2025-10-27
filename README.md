# COBASE example — A copula-based shuffling for ensemble forecast postprocessing

This repository accompanies the paper:

> **COBASE: A new copula-based shuffling method for ensemble weather forecast postprocessing**  
> *Maurits Flos, Bastien François, Irene Schicker, Kirien Whan, and Elisa Perrone*  
> (Eindhoven University of Technology, KNMI, and GeoSphere Austria)

---

## About this project

Weather forecasts are often generated as **ensembles**, i.e., multiple runs of a numerical weather prediction model, to represent uncertainty.  
However, these ensembles usually need **statistical postprocessing** to correct biases and restore realistic dependence structures across variables and locations.

The **COBASE** framework introduced in the paper is a **copula-based postprocessing approach** that combines the strengths of parametric and nonparametric methods through a rank-shuffling mechanism. This design allows COBASE to preserve calibrated univariate forecasts while reconstructing realistic multivariate dependence.

This example uses **mock ensemble-style data** to demonstrate the workflow.

---

## Project Structure

```
COBASE_github/
├── COBASE_github.Rproj           # R Project file
├── main.r                        # Main script controlling the workflow
├── config.yml                    # Configuration file for directories and model settings
│
├── Creating Plots/               # Scripts for generating figures and tables
│   ├── create_plots_paper_selection.R
│   └── Utilities/
│       ├── dm_plots.R
│       ├── skill_scores.R
│       ├── univariate_pit.R
│       ├── score_boxplots.R
│       └── create_score_tables.R
│
├── Postprocessing/               # Postprocessing routines
│   ├── univariate.R
│   ├── multivariate.R
│   └── Utilities/
│       ├── mvpp_methods.R
│       ├── mvpp.R
│       ├── sim_util.R
│       ├── DM_util.R
│       ├── ensfc_util.R
│       ├── scores_util.R
│       ├── uvpp_util.R
│       └── load_data_util.R
│
├── Data/                         # Input and mock datasets
│   ├── Datasets/
│   │   └── Mock_data.csv         # Mock dataset used for demonstration
│   └── ENS/
│       └── transformed_data_Mock_data.Rdata
│
├── Results/                      # Generated results and figures
│   ├── Figures/
│   │   ├── flos_dmcrps_t2m.pdf
│   │   ├── flos_dmesvs_simssh.pdf
│   │   ├── flos_dmesvs_mvpp.pdf
│   │   ├── flos_dmesvs_parcop.pdf
│   │   ├── crps_table.tex
│   │   └── esvs_table.tex
│   ├── Scores/
│   │   └── score_env_Mock_data_mout_10.RData
│   ├── TestStatistic/
│   │   └── dm_statistics_Mock_data_mout_10.Rdata
│   ├── SimilarityMatrix/
│   │   └── simMatrix_Mock_data.Rdata
│   └── UVPP/
│       └── uvpp_Mock_data.Rdata
│
└── README.md
```

---

## ⚙️ Requirements

Tested with **R ≥ 4.2**.

Required R packages:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "xtable", "yaml"))
```

---

## ▶️ Running the example

1. **Open the R Project**
   - Open `COBASE_github.Rproj` in RStudio to automatically set the working directory.

2. **Check the configuration**
   - Review `config.yml` to ensure paths and directories match your local setup.

3. **Run the workflow**
   - Execute `main.r`, which:
     - Loads the mock dataset.
     - Performs univariate and multivariate postprocessing.
     - Generates the scores, diagnostics, and plots under `Results/`.

4. **Inspect results**
   - Figures and tables used in the paper are reproduced under:
     ```
     Results/Figures/
     ```

---

## Notes

- The included data are **mock examples** for demonstration only.  
- In the original study, COBASE was applied to real **ALADIN-LAEF** and **ECMWF** ensemble data.  
- The scripts and structure provided here are designed for clarity and reproducibility.

---

### Note on Code Origin

Portions of this repository are adapted from the [*multiv_pp* package by Lerch et al. (2020)](https://github.com/slerch/multiv_pp). These components have been modified and extended to support the specific case studies and methodological developments presented in the accompanying paper.

**Reference:**  
Lerch, S., Baran, S., Möller, A., Groß, J., Schefzik, R., Hemri, S., & Gräter, M. (2020). *Simulation-based comparison of multivariate ensemble post-processing methods.* *Nonlinear Processes in Geophysics Discussions.* [https://doi.org/10.5194/npg-2019-62](https://doi.org/10.5194/npg-2019-62)

---

## License

This repository is released under the [MIT License](LICENSE).

---

## Contact

For questions or collaboration, please contact:  
**Elisa Perrone** — [e.perrone@tue.nl](mailto:e.perrone@tue.nl)
