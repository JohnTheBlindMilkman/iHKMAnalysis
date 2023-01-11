# iHKM Macros - analysis code for STAR BES and more

This repository contains CERN ROOT macros I used during my PhD when studying p+p correlations using the integrated Hydro-Kinetic Model (iHKM).

## Table of Contents
* [General Info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)
* [Bugs](#bugs)
* [To-Do](#to-do)

## General Info

The [iHKM](#http://cloud-5.bitp.kiev.ua/?p=1247&lang=en) provides a complex treatment  of description of the evolution of the heavy-ion collision (How do I make this sentence less weird?). Together with the Kiev group we (my supervisor and I) test the sutability of the model at lower energies (few GeV). Our testing ground was STAR's BES II FXT programe.

The repository showcases my journey towards that goal, focusing on p-p correlation to complement my other study done for STAR.

## Technologies

* Used OS: Ubuntu 20.04 LTS
* C++ version: 9.4.0
* CERN ROOT version: 6.26/00
* Black Magic (most likely)

## Setup

Setup couldn't be more simple. After cloning the repository cd to the working direcory and use make to complie the program:

```
$ cd iHKMMacros/
$ make
```

## Bugs
* The error bars get scaled by the applied weight, which means that the same function will have different error bars for different interactions included. Which is wrong of course. 

## To-Do

### testiHKM.cc: 
    - [ ] Remove rapidity distributions from the study
### therm2_femto.cxx: 
    - [x] Find the cause of Segmentation Violation
    - [ ] Move methods to other classes
    - [ ] Use std::vector instead of C arrays
    - [ ] Allow to turn on and off all three interactions form the femto.ini file
    - [ ] Move as many parameters as possible (and reasonable) to femto.ini
    - [ ] Store the parameters inside output .root file in a TTree (1 entry and multiple branches)
### therm2_hbtfit.cxx:
    - [ ] Save CF and FIT of Rinv to output hbtfitpipi0a.root files
    - [ ] Make it distinguishible which interactions are included in the correlation function of the output files
### figureHBT.cc:
    - [x] Create a macro for drawing HBT Radii rependence
### figureCorr.cc:
    - [x] Create a macro fro drawing correlation functions
    - [ ] Add 1D option
    - [ ] Add option for different interactions