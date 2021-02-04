# SpindleAnalytics
Analysing mitotic spindle positioning

## Overview 

A data analysis project to analyse mitotic spindle positioning and other spindle parameters using fixed microscopy data. 
Example data and outputs are provided so that the user can run the code and reproduce the plots, or try it on your own data.

### Data organisation 

* The parent directory `SpindleAnalytics` contains subdirectories; `Data`, where the output from ImageJ, the `lookup.csv` and `log.txt` are located; `Output`, where the dataframes and plots are saved to the directories `Dataframe` and `Plots`, respectively; and `Scripts`, where the ImageJ and R code are located. 
* Within `Data` there are subdirectories containing the ImageJ output from separate experiements, each labelled with their unique experiment number (e.g. `JS149`). 
* Analysis was performed blind to the conditions of the experiment so `log.txt` and `lookup.csv` are required and loaded in R to add the original labels.

### Running the code

#### ImageJ macro

* `Spindle_analytics.ijm` requires `ROI_csv.py` to be in Fiji.app/plugins. The macro will take a multi channel image stack and require the user to select the centrosomes, draw around the cell and draw a line through the metaphase plate, these measurements will be saved to a .csv file for analysis in R. Run the script on a directory containing all of the images to be analysed and select a directory to output the data to.

#### R code

* `Spindle_analytics_process.R` will process the csv files (output from ImageJ) for each experiment directory found in `Data`, perform calculations and generate a dataframe that is saved to `Output/Dataframe`. 
* `Spindle_analytics_combine.R` will will combine the dataframes and produce plots that are saved to `Output/Plots`.
