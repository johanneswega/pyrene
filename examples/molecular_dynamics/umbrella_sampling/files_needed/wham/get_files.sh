#!/bin/bash

find "/Users/wega/Library/CloudStorage/OneDrive-unige.ch/PhD/Perylene Bimolecular/Franck Project/MD Simulations/MD/Umbrella Sampling/BPe_DPA_IP_in_DCM/umbrella" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name "umbrella.tpr" \
    | sort -V > tpr-files.dat

find "/Users/wega/Library/CloudStorage/OneDrive-unige.ch/PhD/Perylene Bimolecular/Franck Project/MD Simulations/MD/Umbrella Sampling/BPe_DPA_IP_in_DCM/umbrella" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name "umbrella_pullx.xvg" \
    | sort -V > pullx-files.dat