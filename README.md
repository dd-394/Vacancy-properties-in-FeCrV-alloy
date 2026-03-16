# Vacancy-properties-in-FeCrV-alloy

VASP Data Analysis Tools
This repository contains Python scripts for post-processing VASP (Vienna Ab initio Simulation Package) output files, specifically designed for analyzing vacancy formation energy, atomic displacements, charge density redistribution, and local atomic environments in high-entropy alloys.

Overview
The tools are organized into three modules:

atoms.py – Core classes and functions for handling atomic structures, VASP input/output (POSCAR, CONTCAR, OUTCAR), and basic calculations (distances, energies, forces).
charge.py – Functions to read and normalize CHGCAR charge density files, and compute charge redistribution between two systems.
envirCalc.py – Neighbor search and classification (first to fourth nearest neighbors) with export to Excel.

These scripts were used in the manuscript to process the raw data from 128 vacancy configurations.
