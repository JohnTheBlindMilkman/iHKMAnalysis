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

1. The program breaks after saving the file. Reason is yet unknown.

## To-Do

1. testiHKM.cc: 
    - [ ] Remove rapidity distributions from the study
2. therm2_femto.cxx: 
    - [ ] Find the cause of Segmentation Violation
    - [ ] Move methods to other classes
    - [ ] Use std::vector instead of C arrays
    - [ ] Allow to turn on and off all three interactions form the femto.ini file