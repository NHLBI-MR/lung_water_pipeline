# lung_water_pipeline

This code is for the manuscript (in sub): Seemann, F., et al. Imaging gravity-induced lung water redistribution with automated inline processing at 0.55T MRI

## This repository contains 
- Gadgetron and Matlab code for inline reconstruction and analysis
- Python processing for lung segmentation for analysis, including a pkl file from human and pig data. 
- Segment plugin for offline analysis

# NHLBI TOOLBOX GADGETRON Installation
Repository with gadgetron image reconstruction code: https://github.com/NHLBI-MR/selfgated_noncartesian_reconstruction.git

% GADGETRON MATLAB
- Install application through matlab
- https://www.mathworks.com/matlabcentral/fileexchange/72715-gadgetron
- Make directory "+nhlbi" inside the +gadgetron class
- Copy the lung_water_pipeline.m in to the +nhlbi folder.

# Reconstruction & Inline Lung Water Analysis
- gadgetron configuration [file](https://github.com/NHLBI-MR/lung_water_pipeline/blob/main/gadgetron%20xml%20and%20matlab/spiral_inline_binning_ncc_lungwater.xml)

# Offline analysis using Segment 
Installation
- https://github.com/Cardiac-MR-Group-Lund/segment-open
- [Medviso](https://medviso.com/segment/)


