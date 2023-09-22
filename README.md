# Bayesian Individual-based Compartmental Inference (BICI)

## Introduction

Compartmental models have long been used as a means of understanding the collective dynamics of interacting agents, with notable applications in epidemiology and ecology. BICI allows for arbitrary compartmental model specification and performs simulation and/or inference on that model.
For inference BICI takes in a variety of individual-based or population-based data, priors can be specified from a large range of possibilities, and the outputs consist of posterior trace plots for model parameters, distributions, correlations, visualisations of transitions, dynamic population estimates, and summary statistics (means and 95% credible intervals) as well as diagnostics.

## Features

This multipurpose software allows for sophisticated Bayesian analysis of data using an easy point and click interface.
The current version contains the following features:

* **Arbitrary compartmental models.** The interface allows for easy model specification and can accommodate multiple classifications (e.g. disease status as well as location and sex of individuals). Furthermore multiple species (e.g. predetor-prey models, or environmental accumulation of pathogen).

* **Markovian and Non-Markovian transitions.** This allows for more realistic models (e.g. disease recovery can be modeled using a more adaptable gamma distribution instead of assuming an exponential distribution).

* **Simulation.** The initial conditions can be set and the dynamic variation in the model can be graphically represented in a variety of ways.

* **Inference.** BICI can take a variety of different data types (e.g. state data, population estimates, event times or even uncertain data such as disease diagnostic test results) and infer model parameters as well as underlying model dynamics.

* **Variation in time.** The software allows for the the possibility of model parameters to varying in time according to linear spline functions. 

## Download

The following files can be downloaded:

* **Windows:** [BICI_windows.zip](https://github.com/theITEAM/BICI/releases/download/1.2/BICI_windows.zip). Once unzipped SIRE is run by clicking on the “BICI.exe” icon (if you get the message “Windows protected your PC” you must click “More info” and “Run anyway”).

* **Mac:** [BICI_mac.zip](https://github.com/theITEAM/BICI/releases/download/v1.3/BICI_mac.zip) BICI is run by clicking on the “BICI” icon (if you get the message “Can’t be opened because it’s from an unidentified developer” you must press the control button and click the icon and then select “Open”).


## Documentation

A manual is currently be written.

## Build

To edit and rebuild this software the following instructions must be followed:

* The files from this repository are first downloaded onto your own computer.

* The C++ code in the "src" directory must be compiled on your platform of choice (Windows / Linux / Mac). 

* This software relies of NW.js to run the graphical user interface. This can be downloaded [here](https://github.com/nwjs/nw.js) for your platform of choice.  



