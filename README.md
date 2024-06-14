# OnTheResolutionOfThinBeds
 MATLAB codes for the analysis of resolution of thin beds given a 2D structurally informed prior

The scripts are available with pre-calculated look-up tables with 10000 models in each. In order to do the full analysis please obtain a license for AarhusInv forward modelling. With the codes provided here it is possible to produce lookuptables of the same size as used in the original analysis (10 million models).

# Running the Analysis
To run the analysis please open the MATLAB script named "MainCode.m" and run it. If you do not have a license for AarhusInv and/or if you have not installed it in the same folder as the scripts then you can comment out the section from line 38 to 78 before running the script. In this case only the analysis for reflection seismics will run.