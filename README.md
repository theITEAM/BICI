# Bayesian Individual-based Compartmental Inference (BICI)

## Introduction

Compartmental models have long been used as a means of understanding the collective dynamics of interacting agents, with notable applications in epidemiology, chemistry and ecology. BICI allows for arbitrary compartmental model specification and performs simulation, inference and posterior simulation.
Models can contain multiple interacting species, which can each be treated at a population or individual level. BICI can be run entirely using a point-and-click interface. Alternatively, a scripted language, referred to as “BICI-script”, can be used to construct, store, and export complex models and allow for them to be run on HPC (high performance computing).
For inference BICI takes in a variety of individual and/or population-level data, and priors can be specified from a large range of possibilities. Posterior parameter outputs include trace plots, distributions, correlations, and summary statistics (means and 95% credible intervals) as well as diagnostics. State outputs include various visualisations for populations, transitions and individuals.

## Download

The following files can be downloaded:

* **Windows:** [BICI_v0.79_windows.zip](https://github.com/theITEAM/BICI/releases/download/v0.79/BICI_v0.79_windows.zip). Once unzipped BICI is run by clicking on the “BICI.exe” icon (if the error message “Windows protected your PC” appears, click on “More info” and “Run”).

* **Mac:** 
[BICI_v0.61_mac.zip](https://github.com/theITEAM/BICI/releases/download/v0.61/BICI_v0.61_mac.zip).
 Once unzipped BICI is run by clicking on the “BICI.app” icon (if the error message “BICI can’t be opened because it is from an unidentified developer…” appears, right click on “BICI.app” and select “Open” to give the option to run).

## Features

This multipurpose software allows for sophisticated Bayesian analysis of data using an easy point and click interface.
The current version contains the following features:

* **Arbitrary compartmental models.** The interface allows for easy model specification and can accommodate multiple classifications (e.g. disease status as well as location and sex of individuals). Furthermore multiple species (e.g. predetor-prey models, or environmental accumulation of pathogen).

* **Markovian and Non-Markovian transitions.** This allows for more realistic models (e.g. disease recovery can be modeled using a more adaptable gamma distribution instead of assuming an exponential distribution).

* **Multiple species** This allows for multiple interacting species, e.g. predetor-prey models or pathogen accumulation model.

* **Indivudual and population-based** Allows for mixing population and individual-based species together within the same model.

* **Individual variation** Allows for individual variation (on top of compartmental stratification) by incorporating individual covariates, or individual effects (e.g. for quantitative genetics models).

* **Parallel** An optional parallel implementation can be use to tackle big problems (e.g. on a Linux cluster).

* **Simulation.** The initial conditions can be set and the dynamic variation in the model can be graphically represented in a variety of ways.

* **Inference.** BICI can take a variety of different data types (e.g. state data, population estimates, event times or even uncertain data such as disease diagnostic test results) and infer model parameters as well as underlying model dynamics.

* **Posterior Simulation.** Future prediction, scenario analysis, counterfactual analysis and posterior predictive check. .

* **Variation in time.** The software allows for the the possibility of model parameters to varying in time according to linear spline functions. 

## Documentation

The following manual can be downloaded [here](BICI_Manual_v0.79.pdf)

## Build

To edit and rebuild this software the following instructions must be followed:

* The files from this repository are first downloaded onto your own computer.

* The C++ code in the "src" directory must be compiled on your platform of choice (Windows / Linux / Mac). 

* This software relies of NW.js to run the graphical user interface. This can be downloaded [here](https://github.com/nwjs/nw.js) for your platform of choice.  



