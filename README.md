# iHKM Macros - analysis code for STAR BES and more

This repository contains CERN ROOT macros I used during my PhD when studying p+p correlations using the integrated Hydro-Kinetic Model (iHKM).

## Table of Contents
* [General Info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)
* [Usage](#usage)
* [Bugs](#bugs)
* [To-Do](#to-do)

## General Info
---
The [iHKM](#http://cloud-5.bitp.kiev.ua/?p=1247&lang=en) provides a complex treatment  of description of the evolution of the heavy-ion collision (How do I make this sentence less weird?). Together with the Kiev group we (my supervisor and I) test the sutability of the model at lower energies (few GeV). Our testing ground was STAR's BES II FXT programe.

The repository showcases my journey towards that goal, focusing on p-p correlation to complement my other study done for STAR.

## Technologies
---
* Used OS: Ubuntu 20.04 LTS
* C++ version: 9.4.0
* CERN ROOT version: 6.26/00
* Black Magic (most likely)

## Setup
---
Setup couldn't be more simple. After cloning the repository cd to the working direcory and use make to complie the program:

```
$ cd iHKMMacros/
$ make
```

## Usage
---
### therm2_event:
Program for data simulation.
Same usage as in the original THERMINATOR 2.

### therm2_femto:
Program for particle mixing for femtoscopic correlation analysis.
Currently the program is execuded as follows:
```
./therm2_femto [OPTIONS]
./therm2_femto <INI FILE> <PARENT PID>
```
Possible options contain `-h` / `--help` and `-v` / `--version`. I hope they are self explanatory :)

The .ini file contains all parameters, paths and options to be set before the analysis to obtain desired results. When making your own please make sure to keep the same format, otherwise the reader will not be able retreve the parameters. Defalut .ini file is `./femto.ini`.

The output file path is the same as input. Make sure you have write permission in the input directory.

### therm2_hbtfit:
Program for fitting of the femtoscopic correlation function.
Currently the program is execuded as follows:
```
./therm2_hbtfit [OPTIONS]
./therm2_hbtfit <DIR PATH> <PAIR TYPE> <KT MIN> <KT MAX> [INI FILE]
```
Options and the .ini file idea are the same as in `therm2_femto`. `therm2_hbtfit` uses different .ini file with different options, but stricture is the same.

Currently implemented pair types are *pion-pion* and *pionM-pionM*. Fit only includes quantum statistics. I do not plan to extend this in the fututre.
Defalut .ini file is `./hbtfit.ini`.

### femto_mixer (EXPERIMENTAL):
Program for the same purpose as `therm2_femto`.
With the idea to have it structured as a proper C++ code from 2010s, maybe even 2020s.

Usage is subject to change.

### macros:
Currently the `macro/` folder is a mess and is very much subject to change, feel free to use it, but do not expect it to work flawlessly.

### scripts:
The scripts are:
- divideiHKM.sh - script to divide large data set into separate folders and number them accordingly. Used to submit multiple jobs to cluster.
- jobScriptFemto_SL.sh - job script for GSI Virgo computing cluster (uses Slurm)
- submitFemto.py - submit script for GSI Virgo computing cluster
- wrap.sh - script used to submit jobs to GSI Virgo computing cluster for an older system version

## Bugs
---
## To-Do
---
### therm2_femto.cxx: 
    - [x] Find the cause of Segmentation Violation
    - [x] Move methods to other classes
    - [x] Use std::vector instead of C arrays
    - [ ] Allow to turn on and off all three interactions form the femto.ini file
    - [ ] Move as many parameters as possible (and reasonable) to femto.ini
    - [ ] Store the parameters inside output .root file in a TTree (1 entry and multiple branches)
    - [ ] Make an option to make as little graphs as possible (it seems that the less the histograms the quicker the code)
### therm2_hbtfit.cxx:
    - [ ] Save CF and FIT of Rinv to output hbtfitpipi0a.root files
    - [ ] Make it distinguishible which interactions are included in the correlation function of the output files
    - [ ] Add proper error calculation
### femto_mixer.cxx:
    - [x] Test of it is working
    - [ ] Revert to the old FemtoInteraction methods from therm2_femto
    - [ ] Optimise FemtoInteraction again
    - [ ] Implement all remainign classes
    - [ ] Rename it to therm2_femto.cxx
    - [ ] Clean & debug the code
    - [ ] Think about implementing RDataFrame
