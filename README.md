# Bridge-phyloregions

### Overview

This repository contains scripts and data for the paper "Latitudinal phylogenetic and diversity gradients are explained by a tropical-temperate transitional region." The scripts implement analyses exploring phylogenetic and diversity latitudinal gradients, as well as simulations, focusing on the role of bridge phyloregions in generating latitudinal gradients. 

The primary objective is to reproduce the results presented in the study, including:

1. Global phyloregions of tree species based on distribution and phylogenetic datasets.
2. Phylogenetic alpha and diversity analyses.
3. Environmental relationships with phyloregions using climate datasets.
4. Run and store simulations

## Repository Structure

Scripts
**load_functions.R**
Contains helper functions required for the phyloregions_1.R script.

**phyloregions_1.R**
Performs the empirical analyses, including:

Generating global phyloregions of tree species using:
Species distribution dataset (Lerner et al., 2023)
Phylogenetic tree dataset (Sanchez-Martinez et al., 2023)

Downstream analyses:
Phylogenetic diversity using Standardized Effect Size of Mean Phylogenetic Distance (SES-MPD).
Phylogenetic beta diversity using Rao's Quadratic Entropy (RaoD).

Environmental relationships with phyloregions via random forest models, incorporating:
WorldClim Global Climate data (Hijmans et al., 2005)
OlsenP dataset (McDowell et al., 2023)

**phyloregions_script.R**
Defines functions for running simulations with and without bridge phyloregions.

**phyloregions_script_load.R**
Loads, organizes, and processes simulation outputs for analysis.

