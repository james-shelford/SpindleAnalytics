# SpindleAnalytics
Analysing mitotic spindle positioning

## Overview 

A data analysis project to analyse mitotic spindle positioning and other spindle parameters using fixed microscopy data. Example data and outputs are provided to allow you to run the code and reproduce the plots, or try it on your own data. 

### About the code

The ImageJ script was written to analyse a folder containing multi channel image stacks and requires the user to select the centrosomes, draw around the cell and draw a line through the metaphase plate. The code uses `ROI_csv.py` to save the positions of the selected objects to a csv file and is output to a folder chosen by the user. The R code processes the csv files, performs calculations and generates plots. 

### Data organisation 

* The parent directory `SpindleAnalytics` contains subdirectories; `Data`, where the output from ImageJ and the `lookup.csv` are located; `Output`, where the dataframes and plots are saved to the directories `Dataframe` and `Plots`, respectively; and `Scripts`, where the ImageJ and R code are located. 
* Within `Data` there are subdirectories containing the ImageJ output from separate experiments, each labelled with their unique experiment number (e.g. `JS149`). 
* Analysis was performed blind to the conditions of the experiment so `log.txt` and `lookup.csv` are required and loaded in R to add the original labels.

### Running the code

1. Add `Spindle_analytics.ijm` and `ROI_csv.py` to *Fiji/plugins*
2. Run `Spindle analytics` in Fiji, select the directory to be analysed and follow the on screen instructions
3. Run `Spindle_analytics_process.R` and select the directory containing the output from ImageJ. R will perform the calculations and generate a dataframe containing the parameters that is saved to *Output/Dataframe*
4. Run `Spindle_analytics_combine.R` to combine the dataframes from multiple experiments and produce plots that are saved to *Output/Plots*

