# Bayesian Individual-based Compartmental Inference (BICI)

## Introduction

Compartmental models have long been used as a means of understanding the collective dynamics of interacting agents, with notable applications in epidemiology, chemistry and ecology. BICI allows for arbitrary compartmental model specification and performs simulation, inference and posterior simulation.
BICI can be run entirely using a point-and-click interface. Alternatively, a scripted language, referred to as “BICI-script”, can be used to construct, store, and export complex models and allow for them to be run on HPC (high performance computing).

For inference BICI takes in a variety of individual and/or population-level data, and priors can be specified from a large range of possibilities. Posterior parameter outputs include trace plots, distributions, correlations, and summary statistics (means and 95% credible intervals) as well as diagnostics. State outputs include various visualisations for populations, transitions and individuals.

## Download

BICI can be working for you in just a few clicks! The following pre-compiled versions can be downloaded:

* **Windows:** [BICI_v0.91_windows.zip](https://github.com/theITEAM/BICI/releases/download/v0.91/BICI_v0.91_windows.zip). Once unzipped BICI is run by clicking on the “BICI.exe” icon. If the error message “Windows protected your PC” appears, click on “More info” and “Run” (if this still doesn’t work then open a command prompt within the BICI folder, e.g. by typing “cmd” into the Windows File Explorer address bar, and then execute the command ‘bici’).

* **Mac:** 
[BICI_v0.91_mac.zip](https://github.com/theITEAM/BICI/releases/download/v0.91/BICI_v0.91_mac.zip). Once unzipped, BICI is run by clicking on the “BICI.app” icon (if the error message “BICI can’t be opened because it is from an unidentified developer…” appears, right click on “BICI.app” and select “Open” to give the option to run).

* **Linux:** 
[BICI_v0.92_linux.tar.gz](https://github.com/theITEAM/BICI/releases/download/v0.92/BICI_v0.92_linux.tar.gz). Once downloaded, this file can be extracted by using the terminal command “tar -zxvf BICI_v0.92_linux.tar.gz”. The code is executed using “./BICI”.

## Documentation

The manual can be downloaded [here](BICI_Manual_v0.90.pdf). This explains all the features within the software and well as comprehensive documentation for the BICI-script language.  

## Features

* **Arbitrary compartmental models** — The interface allows for easy model specification and can accommodate an arbitrary number of compartments in multiple classifications (e.g. disease status, location and/or demographic grouping).

* **Markovian and Non-Markovian transitions** — Allows for more realistic models (e.g. disease recovery can be modelled using a more adaptable gamma distribution instead of assuming an exponential distribution).

* **Multiple species** — Allows for multiple interacting species (e.g. predator-prey models, or pathogen accumulation models, modelling of disease vectors).

* **Individual and population-based** — Allows for combining population and individual-based species together within the same model.

* **Individual variation** — Allows for individual variation (on top of compartmental stratification) by incorporating individual-based covariates, or individual effects (e.g. for quantitative genetics models).

* **Parallel** — A parallel implementation can be used to tackle big problems (e.g. on a Linux cluster).

* **Variation in time.** —Allows for the possibility of model parameters to vary in time according to spline functions. 

* **Scripted language** — Referred to as “BICI-script”, this can be used to construct, store, and export complex models and allow for them to be run on HPC.

## Modes of operation
* **Simulation** — The initial conditions can be set and the dynamic variation in the model can be graphically represented in a variety of ways.

* **Inference** — BICI can take a variety of different data types (e.g. state data, population estimates, event times, disease diagnostic test results, individual covariates, genomic relationship matrix data, pathogen genetic data) and infer model parameters as well as underlying model dynamics.

* **Posterior Simulation** — Future prediction, scenario analysis, counterfactual analysis and posterior predictive checking. 

## Build

Note, this is only recommended for developers of the code, and users are strong encourage to download the pre-compiled versions above.

To edit and rebuild this software the following instructions must be followed:

* The files from this repository are first downloaded onto your own computer.

* The C++ code in the "src" directory must be compiled on your platform of choice (Windows / Linux / Mac). 

* This software relies of NW.js to run the graphical user interface. This can be downloaded [here](https://github.com/nwjs/nw.js) for your platform of choice.  

## Contact

If you have any problems running this software (it is currently under development), please contact Chris Pooley at [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk).


