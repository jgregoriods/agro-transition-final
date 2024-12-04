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
2. Install dependencies (see `main.R` for required packages).
3. Download the QuickMEM source code from <a href="https://github.com/ajsmit/Quantitative_Ecology/blob/main/Num_Ecol_R_book_ed1/quickMEM.R" target="_blank">here</a>.
4. Run the `main.R` script to reproduce the results.

## Citation
If you use this code or data, please cite:

> Jonas Gregorio de Souza, Javier Ruiz-Pérez, Abel Ruiz-Giralt, Carla Lancelotti, Marco Madella. *Modeling the global transition to food-producing economies*. 2024.
