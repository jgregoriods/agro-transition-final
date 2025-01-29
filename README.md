# Modeling the Global Transition to Food-Producing Economies

## Overview
This repository contains the code and data for the paper *Modeling the Global Transition to Food-Producing Economies* by Jonas Gregorio de Souza, Javier Ruiz-Pérez, Abel Ruiz-Giralt, Carla Lancelotti, and Marco Madella. The study investigates the shift from foraging to plant cultivation during the Holocene, analyzing spatial and environmental factors that influenced the timing of agricultural adoption worldwide. A machine-learning approach using Random Survival Forests was employed to analyze radiocarbon dates and environmental predictors.

## Files in This Repository

### Data
- **Radiocarbon Data**: Radiocarbon dates linked to the first adoption of agriculture.
- **Environmental Predictors**: Data on bioclimatic, edaphic, terrain, and vegetation variables derived from paleoclimate simulations and other sources.

### Code
- **`main.R`**: Contains the full implementation of the Random Forest model used for survival analysis. The script includes:
  - Data preprocessing.
  - Model fitting and tuning.
  - Evaluation using Harrell's concordance index.
  - Variable importance analysis with SHAP values.

## How to Use
1. Clone or download this repository.
2. Download the QuickMEM source code from <a href="https://github.com/ajsmit/Quantitative_Ecology/blob/main/Num_Ecol_R_book_ed1/quickMEM.R" target="_blank">here</a> and place it in the same folder as the `main.R` file.
3. Run the `main.R` script to reproduce the results. The script will:
   - Install the required packages.
   - Load the data.
   - Preprocess the data.
   - Fit the Random Survival Forest model.
   - Evaluate the model.
   - Generate SHAP values for variable importance analysis.

## Citation
If you use this code or data, please cite:

> Jonas Gregorio de Souza, Javier Ruiz-Pérez, Abel Ruiz-Giralt, Carla Lancelotti, Marco Madella. *Modeling the global transition to food-producing economies*. 2024.

## Attribution

This repository is primarily licensed under the **MIT License** (see [`LICENSE`](LICENSE)).

However, some third-party content included in this repository is licensed under [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](https://creativecommons.org/licenses/by-nc-sa/4.0/):

- **Soil Rasters** in `rasters/edaphic` downloaded from the [Harmonized World Soil Database](https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v20/en/) by the Food and Agriculture Organization (FAO) and IIASA.
- **Paleoclimate Data** in `rasters/bioclim` was averaged over the Holocene time slices from [Beyer, R.M., Krapp, M. & Manica, A. "High-resolution terrestrial climate, bioclimate, and vegetation for the last 120,000 years", Scientific Data 7, 236 (2020)](https://www.nature.com/articles/s41597-020-0552-1). Accessed via the [pastclim R package](https://github.com/EvolEcolGroup/pastclim).

Users must comply with the relevant license terms when using these datasets.