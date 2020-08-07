# Bayesian Individual-based Compartmental Inference (BICI)

## Introduction

Compartmental models have long been used as a means of understanding the collective dynamics of interacting agents, with notable applications in epidemiology and ecology. BICI allows for arbitrary compartmental model specification and performs simulation and/or inference on that model.
For inference BICI takes in a variety of individual-based or population-based data, priors can be specified from a large range of possibilities, and the outputs consist of posterior trace plots for model parameters, distributions, correlations, visualisations of transitions, dynamic population estimates, and summary statistics (means and 95% credible intervals) as well as diagnostics.

## Features

This multipurpose software allows for sophisticated Bayesian analysis of data using an easy point and click interface.
The current version (v1.0) contains the following features:

* **Arbitrary compartmental models.** The interface allows for easy model specification and can accommodate multiple classifications (e.g. disease status as well as location and sex of individuals).

* **Markovian and Non-Markovian transitions.** This allows for more realistic models (e.g. disease recovery can be modeled using a more adaptable gamma distribution instead of assuming an exponential distribution).

* **Simulation.** The initial conditions can be set and the dynamic variation in the model can be graphically represented in a variety of ways.

* **Inference.** BICI can take a variety of different data types (e.g. state data, population estimates, event times or even uncertain data such as disease diagnostic test results) and infer model parameters as well as underlying model dynamics.

* **Variation in time.** The software allows for the the possibility of model parameters discretely varying in time. Also variation in parameters with the age of individuals can be incorporated.

## Download

The following files can be downloaded:

* **Windows:** [BICI_v1.0_windows.zip](https://github.com/BioSS-EAT/SIRE/releases/download/v1.0/BICI_v1.0_windows.zip). Once unzipped SIRE is run by clicking on the “BICI.exe” icon.

* **Linux:** BICI_v1.0_linux.tar.gz(Currently not availalbe). This can then be extracted by using the terminal command “tar -zxvf BICI_v1.0_linux.tar.gz”. The code is executed using “./BICI”.

* **Mac:** [BICI_v1.0_Mac](https://github.com/BioSS-EAT/BICI/releases/download/v1.0/BICI_v1.0_Mac) BICI is run by clicking on the “BICI.app” icon.

## Documentation

Information about the software can be obtained from the [BICI manual](https://github.com/BioSS-EAT/BICI/blob/master/BICI_Manual_v1.0.pdf). A [website](https://bioss-eat.github.io/BICI.html) provides screenshots of the software.

## Build

To edit and rebuild this software the following instructions must be followed:

* The files from this repository are first downloaded onto your own computer.

* The C++ code in the "Execute" directory must be compiled on your platform of choice (Windows / Linux / Mac). For example "g++ bici.cc tinyxml2.cc -o a.exe -O3" can be used. The resulting "a.exe" executable file is placed into the "Execute" directory.

* This software relies of NW.js to run the graphical user interface. This can be downloaded [here](https://github.com/nwjs/nw.js) for your platform of choice.  

* All the files and folders downloaded from this repository are copies directly into the NW folder. 

* BICI is run by clicking on "NW.exe".


