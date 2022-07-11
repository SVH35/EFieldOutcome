# EFO: EFieldOutcome

![EFO Frontpage](/GitHubFigure.jpg)

## Introduction
The main goal of EFO is to provide a tool for researchers to provide outcome measure homogeneity across studies. Therefore, the script provides both the E-field in a grey matter region of interest (ROI), as well as the 99th percentile (P99) E-field value across the whole grey matter. Thus, it supplies the two most-used outcome measures in the field. 
Moreovero, EFO also provides the volumetric properties (mm3) of the ROI and the P99 E-fields, as well as the overlap between both volumes and the E-field ratio, which is the ratio of ROI E-field to P99 E-field. 

## Authors
The following persons worked on this project
* Sybren Van Hoornweder
* Raf L.J. Meesen
* Kevin A. Caulfield

## How to cite
Pleaser cite the following article if you use the following script: TBP

## Requirements
The function requires [MATLAB](https://www.mathworks.com/products/matlab.html) and [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) to run. 

## How to run
To run the pipeline, the following command should be used:
outcome = EFieldExtraction(pathname, mni_coord, roi_sphere, toplot);
  pathname = path that contains .msh simulation file
  mni_coord = MNI coordinate of grey matter ROI that should be analyzed
  roi_sphere = size of sphere used for ROI analyses, typically 10mm is used
  toplot = can be either 1 (show plot) or 0 (don't show plot)

To run the pipeline, download the [code](/Code). 

## License
This software runs under a GNU General Public License v3.0.

This software uses free packages from the Internet, except MATLAB, which is a proprietary software by the MathWorks. You need a valid Matlab license to run this software.
