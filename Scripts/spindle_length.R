# Script to load in .csv files containing different (x,y) coords for centrosome positions
# Make a loop to read the csv files for each ROI, scale and add the data to a dataframe along with the name of the cell. 

require(ggplot2)
require(ggbeeswarm)
require(Hmisc)
library(plyr)
library(dplyr)
library(multcomp)
library(RColorBrewer)
library(cowplot)

# List the .csv files
datadir <- rstudioapi::selectDirectory()
my_files_centrosomes <- list.files(datadir, pattern = "*centrosomes.csv",full.names = TRUE)

# Removing stuff from the file path and name
file_name <- basename(my_files_centrosomes)
file_name <- gsub("_centrosomes.csv","", file_name)

# get the experiment name
expt <- basename(datadir)

# Create a dataframe that will be used to store the results of the calculations.
headings <- c('centrosome_centrosome_dist', 'file_name')
output <- matrix(0,length(file_name),length(headings))
colnames(output) <- headings
output <- as.data.frame(output)

#######################
## function definitions

# This function requires the filenames as input and reads the csv file, output as matrix

load_the_data <- function(the_filename){
  # import data
  the_data <- read.csv(file = the_filename, header=FALSE, stringsAsFactors=FALSE)
  the_data <- data.matrix(the_data)
}


####################
# Load and process data

if(grepl("ELR", datadir) == TRUE) {
  pxpum <- 8.6207
} else {
  pxpum <- 5.4975
}

# Loop through the number of files 
for(i in 1:length(my_files_centrosomes)){
  my_filename_centrosomes <- my_files_centrosomes[i]
  # Call the function to get the data
  centrosome_data <- load_the_data(my_filename_centrosomes)
  output$file_name[i] <- file_name[i]
  a1 <- centrosome_data[1,1] - centrosome_data[2,1]
  b1 <- centrosome_data[1,2] - centrosome_data[2,2]
  centrosome_dist <- sqrt((a1 * a1) + (b1 * b1))
  # All values are in pixels, need to convert to microns.
  output$centrosome_centrosome_dist[i] <- centrosome_dist / pxpum
}

## The data need to be unblinded and categorised

# load the log.txt file
blind_log <- read.delim(paste(datadir,"log.txt", sep = "/"), header = TRUE) 

# We need a look up table

# function to find partial strings in a column and classify them
add_categories = function(x, patterns, replacements = patterns, fill = NA, ...) {
  stopifnot(length(patterns) == length(replacements))
  ans = rep_len(as.character(fill), length(x))
  empty = seq_along(x)
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]  
    empty = empty[!greps]
  }
  return(ans)
}

look_up_table <- read.table("Data/lookup.csv", header = TRUE, stringsAsFactors = F, sep = ",")

# add a new column to dataframe where categories are defined by searching original name for partial strings
blind_log$Category <- add_categories(blind_log$Original_Name,
                                     look_up_table$Search_name,
                                     look_up_table$Search_category,
                                     "NA", ignore.case = TRUE)

# Now we have a dataframe that can be used to lookup the real values

# This line looks up the correct Category from the blind_log
output$Category <- with(blind_log,
                        Category[match(output$file_name,
                                       Blinded_Name)])

################
## Plot the data

# Plotting the Categories in a specific order
output$Category <- as.factor(output$Category)
#output$Category <- factor(output$Category, levels = look_up_table$Search_category)

make_the_plot <- function(df_for_plotting, xCol, yCol, yLab){
  thePlot <- ggplot(data = df_for_plotting, aes_string(x = xCol, y = yCol, color = xCol)) +
    geom_quasirandom(alpha = 0.5, stroke=0, dodge.width = 1) + 
    stat_summary(fun.data = mean_se, geom = 'point', size = 2, aes(group = Category)) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom = 'errorbar', size = 0.8, aes(group = Category), width = 0) +
    theme(axis.text.x = element_text(face = "plain", color = 'black', size = 8, angle = 0, hjust = 0.5), axis.text.y = element_text(face = 'plain', color = 'black', size = 8)) +
    theme(axis.title.y = element_text(size = 10,face = 'plain', color = 'black')) +
    labs(y = yLab, x = NULL) + 
    theme(legend.title = element_blank()) + 
    theme(legend.position = 'none') 
  
  return(thePlot)
}

# make the plots
spindle_length_plot <- make_the_plot(output, "Category", "centrosome_centrosome_dist", "Spindle length (microns)")
ggsave(paste0("Output/Plots/",expt,"_spindle_length.pdf"), plot = spindle_length_plot)

# make directory for Data if it doesn't exist
ifelse(!dir.exists("Output"), dir.create("Output"), "Folder exists already")
ifelse(!dir.exists("Output/Dataframe"), dir.create("Output/Dataframe"), "Folder exists already")

# Add the experiment number to the df and save
output$Experiment_number <- expt
saveRDS(output, file = paste0("Output/Dataframe/length_", expt, ".rds"))
