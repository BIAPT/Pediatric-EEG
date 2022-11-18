# Pediatric-EEG
This repository is used for preprocessing and analyzing 26-channel EEG data 
- Data provided by Dr. Kevin Jones

Requires following add-ons: 
- Symbolic Math Toolbox
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox
- EEGapp

## Main pipeline scripts
- Spectrogram_check
- PD_EEG_pipeline_baseline: smallest baseline threshold acquired
- PD_EEG_pipeline_rest

## Preprocessing steps in MATLAB
- step-by-step instruction uploaded as the pdf file

## Analysis Pipeline in MATLAB
- derived from the NET_ICU Pipeline by Charlotte and Miriam
- Utils include csv files of left, right, whole brain electrode locations in 'BESA'
- Topographic Map measuring from 1-15Hz (delta, theta, alpha) frequency bands
