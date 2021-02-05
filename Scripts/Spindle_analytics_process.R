# Script written by James Shelford
# Script to process .csv files containing (x,y) positions of centrosomes, cell outline and metaphase plate generated from ImageJ 'Spindle_analytics'
# The output is a dataframe saved as .rds for combining with other experiments to further process

# Working directory should be 'SpindleAnalytics'
# Lookup.csv should be in the 'Data' subdirectory
# Setup preferred directory structure in wd
ifelse(!dir.exists("Data"), dir.create("Data"), "Folder exists already")
ifelse(!dir.exists("Output"), dir.create("Output"), "Folder exists already")
ifelse(!dir.exists("Output/Dataframe"), dir.create("Output/Dataframe"), "Folder exists already")
ifelse(!dir.exists("Output/Plots"), dir.create("Output/Plots"), "Folder exists already")
ifelse(!dir.exists("Scripts"), dir.create("Scripts"), "Folder exists already")

# Select directory containing the .csv files to be processed
datadir <- rstudioapi::selectDirectory()
# List the .csv files
my_files_centrosomes <- list.files(datadir, pattern = "*centrosomes.csv",full.names = TRUE)
my_files_outlines <- gsub("_centrosomes.csv", "_outline.csv", my_files_centrosomes)
my_files_DNA <- gsub("_centrosomes.csv", "_DNA.csv", my_files_centrosomes)

# Removing stuff from the file path and name
file_name <- basename(my_files_centrosomes)
file_name <- gsub("_centrosomes.csv","", file_name)

# Extract the experiment number for use later (useful when combining experiments)
expt <- basename(datadir)

# Create a dataframe that will be used to store the results of the calculations.
headings <- c('intersect_1x', 'intersect_1y', 'intersect_2x', 'intersect_2y', 'lateral_intersect_1x', 
              'lateral_intersect_1y', 'lateral_intersect_2x', 'lateral_intersect_2y', 'centrosome_Ax','centrosome_Ay',
              'centrosome_Bx', 'centrosome_By', 'd1', 'd2', 'd1_d2_offset','metaphase_Ax', 'metaphase_Ay', 'metaphase_Bx', 
              'metaphase_By','metaphase_intersect_X', 'metaphase_intersect_Y', "D1", 'D2','D1_D2_offset', 'l1', 'l2', 'l1_l2_offset', 
              'centrosome_centrosome_dist', 'cell_length', 'cell_width', 'Centre_X', 'Centre_Y', 'offset_dist', 'spindle_width', 
              'spindle_aspect_ratio', 'cell_aspect_ratio', "spindle_angle", 'file_name')
output <- matrix(0,length(file_name),length(headings))
colnames(output) <- headings
output <- as.data.frame(output)

###------------------------------ Function Definitions ------------------------------

# This function requires the filename as input and reads the csv file, output as matrix
load_the_data <- function(the_filename){
  
  the_data <- read.csv(file = the_filename, header=FALSE, stringsAsFactors=FALSE)
  the_data <- data.matrix(the_data)
}

# Put the coords into a dataframe along with the filename
store_basic_info <- function(df, centrosome, dna, rownum) {
  
  df$file_name[rownum] <- file_name[rownum]
  df$centrosome_Ax[rownum] <- centrosome[1,1]
  df$centrosome_Ay[rownum] <- centrosome[1,2]
  df$centrosome_Bx[rownum] <- centrosome[2,1]
  df$centrosome_By[rownum] <- centrosome[2,2]
  df$metaphase_Ax[rownum] <- dna[1,1]
  df$metaphase_Ay[rownum] <- dna[1,2]
  df$metaphase_Bx[rownum] <- dna[2,1]
  df$metaphase_By[rownum] <- dna[2,2]
  
  return(df)
}

# Want to make a function that takes a matrix containing the (x,y) coords of the two centrosomes and a matrix containing (x,y) coords of outline as input. 
# Cell outline is lots of little lines (line segments), centrosome line is constant (treat as infinite line). 
# Given two points on a line we can make equation. We will have two lines so two equations, solve to find intersect. 
# check that intersection is within the line segment, if it is store it.
# There will be two points of intersection, so will need a counter to differentiate between the two. 
# Define A1, B1 and C1, and A2, B2, and C2. The first 3 will remain constant (centrosome line), while the remaining three will change as the (x,y) coords of the outline change.

find_intersections <- function(df, rownum, outline_df, object) {
  
  if (object == "centrosome") {
    x1 <- df$centrosome_Ax[rownum]
    x2 <- df$centrosome_Bx[rownum]
    y1 <- df$centrosome_Ay[rownum]
    y2 <- df$centrosome_By[rownum]
  } else {
    x1 <- df$metaphase_Ax[rownum]
    x2 <- df$metaphase_Bx[rownum]
    y1 <- df$metaphase_Ay[rownum]
    y2 <- df$metaphase_By[rownum] 
  }
  
  # Centrosome
  A1 <- y2 - y1
  B1 <- x1 - x2
  C1 <- A1 * x1 + B1 * y1
  
  # Is the point of intersection within the line segment of the outline?
  # If the counter == 4 the intersection lies within the segment. Output this coord to the dataframe. 
  
  intersect_count <- 0
  final_row <- nrow(outline_df)
  
  for (row_num in 1:final_row) { # Loop through the rows of the outline file to get pairs of coordinates to test (line segments)
    
    counter <- 0 # the count used to check if the intersect lies within the segment 
    
    if(row_num == final_row){
      x3 <- outline_df[1,1]
      x4 <- outline_df[final_row,1]
      y3 <- outline_df[1,2]
      y4 <- outline_df[final_row,2]
    } else {
      x3 <- outline_df[row_num,1]
      x4 <- outline_df[row_num+1,1]
      y3 <- outline_df[row_num,2]
      y4 <- outline_df[row_num+1,2]
    }
    
    # Outline
    A2 = y4 - y3
    B2 = x3 - x4
    C2 = A2 * x3 + B2 * y3
    
    determ <- A1 * B2 - A2 * B1
    X <- (B2 * C1 - B1 * C2) / determ
    Y <- (A1 * C2 - A2 * C1) / determ
    
    if(is.na(X) || is.na(Y)) {
      next
    }
    
    if(x3 <= X || x4 <= X){
      counter <- counter + 1
    }
    if(x3 >= X || x4 >= X){
      counter <- counter + 1
    }
    if(y3 <= Y || y4 <= Y){
      counter <- counter + 1
    }
    if(y3 >= Y || y4 >= Y){
      counter <- counter + 1
    }
    
    if (counter == 4){
      intersect_count <- intersect_count + 1
      
      # output the intersect
      
      if (intersect_count == 1) {
        if (object == "centrosome") {
          df$intersect_1x[rownum] <- X
          df$intersect_1y[rownum] <- Y
        } else {
          df$lateral_intersect_1x[rownum] <- X
          df$lateral_intersect_1y[rownum] <- Y
        }
      }
      if (intersect_count == 2){
        if (object == "centrosome") {
          df$intersect_2x[rownum] <- X
          df$intersect_2y[rownum] <- Y
          if (df$intersect_2x[rownum] != df$intersect_1x[rownum] && df$intersect_2y[rownum] != df$intersect_1y[rownum]) {
            break
          }
        } else {
          df$lateral_intersect_2x[rownum] <- X
          df$lateral_intersect_2y[rownum] <- Y
          if (df$lateral_intersect_2x[rownum] != df$lateral_intersect_1x[rownum] && df$lateral_intersect_2y[rownum] != df$lateral_intersect_1y[rownum]) {
            break
          }
        }
      }
      # the following is needed if intersect coincides with an outline coord
      # then the same point is stored for both intersects, this gives width of 0 and ratio of inf
      # so we just allow the code to look for the next intersect
      if (intersect_count == 3){
        if (object == "centrosome") {
          df$intersect_2x[rownum] <- X
          df$intersect_2y[rownum] <- Y
          break
        } else {
          df$lateral_intersect_2x[rownum] <- X
          df$lateral_intersect_2y[rownum] <- Y
          break
        }
      }
    }
  }
  
  return(df)
}

# Input df containing coords of intersections of centrosome line and cell outline
# calculate the distance between the centrosome and the cell outline. Output distance to the dataframe

d1_d2_func <- function(df, i, pixels_per_micron){
  
  # Distance between intersect 1 and both centrosome points, store the lowest number as d1
  # Distance between intersect 2 and both centrosome points, store the lowest number as d2
  
  # Pythagoras theorem a2 + b2 = c2
  
  intersect_1_x <- output$intersect_1x[i]
  intersect_1_y <- output$intersect_1y[i]
  intersect_2_x <- output$intersect_2x[i]
  intersect_2_y <- output$intersect_2y[i]
  
  Centrosome_A_x <- output$centrosome_Ax[i]
  Centrosome_A_y <- output$centrosome_Ay[i]
  Centrosome_B_x <- output$centrosome_Bx[i]
  Centrosome_B_y <- output$centrosome_By[i]
  
  a1 <- intersect_1_x - Centrosome_A_x
  b1 <- intersect_1_y - Centrosome_A_y
  a2 <- intersect_2_x - Centrosome_A_x
  b2 <- intersect_2_y - Centrosome_A_y
  
  distance_1 <- sqrt((a1*a1)+(b1*b1))
  distance_2 <- sqrt((a2*a2)+(b2*b2))
  
  if(distance_1<distance_2){
    d1 <- distance_1
  }
  else{
    d1 <- distance_2
  }
  
  a3 <- intersect_1_x - Centrosome_B_x
  b3 <- intersect_1_y - Centrosome_B_y
  a4 <- intersect_2_x - Centrosome_B_x
  b4 <- intersect_2_y - Centrosome_B_y
  
  distance_3 <- sqrt((a3*a3)+(b3*b3))
  distance_4 <- sqrt((a4*a4)+(b4*b4))
  
  if(distance_3 < distance_4){
    d2 <- distance_3
  }
  else{
    d2 <- distance_4
  }
  
  # All values are in pixels, need to convert to microns.
  d1_d2_offset <- abs(d1-d2)/pixels_per_micron
  
  output$d1_d2_offset[i] <- d1_d2_offset
  output$d1[i] <- d1
  output$d2[i] <- d2
  return(output)
}

# Function to calculate the two distances (l1 and l2) between the metaphase plate and the cell outline (laterally).
# Requires the df as input to access the two (x,y) coords of the outline intersects and the metaphase plate endpoints.

l1_l2_func <- function(df, i, pixels_per_micron){
  
  # Distance between intersect 1 and both metaphase plate endpoints, store the lowest number as l1
  # Distance between intersect 2 and both metaphase plate endpoints, store the lowest number as l2
  
  intersect_1_x <- output$lateral_intersect_1x[i]
  intersect_1_y <- output$lateral_intersect_1y[i]
  intersect_2_x <- output$lateral_intersect_2x[i]
  intersect_2_y <- output$lateral_intersect_2y[i]
  
  Metaphase_A_x <- output$metaphase_Ax[i]
  Metaphase_A_y <- output$metaphase_Ay[i]
  Metaphase_B_x <- output$metaphase_Bx[i]
  Metaphase_B_y <- output$metaphase_By[i]
  
  a1 <- intersect_1_x - Metaphase_A_x
  b1 <- intersect_1_y - Metaphase_A_y
  a2 <- intersect_2_x - Metaphase_A_x
  b2 <- intersect_2_y - Metaphase_A_y
  
  distance_1 <- sqrt((a1*a1)+(b1*b1))
  distance_2 <- sqrt((a2*a2)+(b2*b2))
  
  if(distance_1<distance_2){
    l1 <- distance_1
  }
  else{
    l1 <- distance_2
  }
  
  a3 <- intersect_1_x - Metaphase_B_x
  b3 <- intersect_1_y - Metaphase_B_y
  a4 <- intersect_2_x - Metaphase_B_x
  b4 <- intersect_2_y - Metaphase_B_y
  
  distance_3 <- sqrt((a3*a3)+(b3*b3))
  distance_4 <- sqrt((a4*a4)+(b4*b4))
  
  if(distance_3 < distance_4){
    l2 <- distance_3
  }
  else{
    l2 <- distance_4
  }
  
  # All values are in pixels, need to convert to microns.
  l1_l2_offset <- abs(l1-l2) / pixels_per_micron
  
  output$l1_l2_offset[i] <- l1_l2_offset
  output$l1[i] <- l1
  output$l2[i] <- l2
  return(output)
}

# Find the point of intersection of the centrosome line and metaphase plate. 
# This point will be considered as the centre of the spindle, and will be used in subsequent calculations 
# This point will also be used to calculate distance to cell outline. 
# Input the df containing values needed for calculation and i to output the calculation to the correct row.

metaphase_intersect <- function(df, i){
  
  # Define the parameters, (x1,y1) and (x2,y2) are the two centrosomes. (x3,y4) and (x4,y4) are the two end points of the metaphase plate. 
  
  x1 = output$centrosome_Ax[i]
  x2 = output$centrosome_Bx[i]
  y1 = output$centrosome_Ay[i]
  y2 = output$centrosome_By[i]
  
  x3 = output$metaphase_Ax[i]
  x4 = output$metaphase_Bx[i]
  y3 = output$metaphase_Ay[i]
  y4 = output$metaphase_By[i]
  
  # Centrosome
  A1 = y2 - y1
  B1 = x1 - x2
  C1 = A1*x1 + B1*y1
  
  # Outline
  A2 = y4 - y3
  B2 = x3 - x4
  C2 = A2*x3 + B2*y3
  
  determ <- A1*B2 - A2*B1
  X <- (B2*C1 - B1*C2)/determ
  Y <- (A1*C2 - A2*C1)/determ
  
  output$metaphase_intersect_X[i] <- X
  output$metaphase_intersect_Y[i] <- Y
  
  return(output)
  
}

D1_D2_func <- function(df, i, pixels_per_micron){
  
  # D1 = distance between intersect 1 and the metaphase intersect point
  # D2 = distance between intersect 2 and the metaphase intersect point
  
  # Pythagoras theorem a2 + b2 = c2
  
  intersect_1_x <- output$intersect_1x[i]
  intersect_1_y <- output$intersect_1y[i]
  intersect_2_x <- output$intersect_2x[i]
  intersect_2_y <- output$intersect_2y[i]
  
  meta_intersect_x <- output$metaphase_intersect_X[i]
  meta_intersect_y <- output$metaphase_intersect_Y[i]
  
  a1 <- intersect_1_x - meta_intersect_x
  b1 <- intersect_1_y - meta_intersect_y
  a2 <- intersect_2_x - meta_intersect_x
  b2 <- intersect_2_y - meta_intersect_y
  
  D1 <- sqrt((a1*a1)+(b1*b1))
  D2 <- sqrt((a2*a2)+(b2*b2))
  
  # All values are in pixels, need to convert to microns.
  D1_D2_offset <- abs(D1-D2) / pixels_per_micron
  
  output$D1_D2_offset[i] <- D1_D2_offset
  output$D1[i] <- D1
  output$D2[i] <- D2
  
  return(output)
}

# Function to calculate the distance between the two centrosomes

centrosome_dist_func <- function(df, i, pixels_per_micron){
  
  # Pythagoras theorem a2 + b2 = c2
  
  centrosome_A_x <- output$centrosome_Ax[i]
  centrosome_A_y <- output$centrosome_Ay[i]
  centrosome_B_x <- output$centrosome_Bx[i]
  centrosome_B_y <- output$centrosome_By[i]
  
  a1 <- centrosome_A_x - centrosome_B_x
  b1 <- centrosome_A_y - centrosome_B_y
  
  centrosome_dist <- sqrt((a1*a1)+(b1*b1))
  # All values are in pixels, need to convert to microns.
  centrosome_dist <- centrosome_dist / pixels_per_micron
  
  output$centrosome_centrosome_dist[i] <- centrosome_dist
  
  return(output)
}

# Function to calculate the distance between the two points of the cell outline where the spindle axis intersects, this is cell length

cell_length_func <- function(df, i, pixels_per_micron){
  
  # Pythagoras theorem a2 + b2 = c2
  
  intersect_1_x <- output$intersect_1x[i]
  intersect_1_y <- output$intersect_1y[i]
  intersect_2_x <- output$intersect_2x[i]
  intersect_2_y <- output$intersect_2y[i]
  
  a1 <- intersect_1_x - intersect_2_x
  b1 <- intersect_1_y - intersect_2_y
  
  cell_length <- sqrt((a1*a1)+(b1*b1))
  # All values are in pixels, need to convert to microns.
  cell_length <- cell_length / pixels_per_micron
  
  output$cell_length[i] <- cell_length
  
  return(output)
}

# Function to calculate the distance between the two points of the cell outline where the metaphase plate line intersects, this is the cell width
cell_width_func <- function(df, i, pixels_per_micron){
  
  # Pythagoras theorem a2 + b2 = c2
  
  lateral_intersect_1_x <- output$lateral_intersect_1x[i]
  lateral_intersect_1_y <- output$lateral_intersect_1y[i]
  lateral_intersect_2_x <- output$lateral_intersect_2x[i]
  lateral_intersect_2_y <- output$lateral_intersect_2y[i]
  
  a1 <- lateral_intersect_1_x - lateral_intersect_2_x
  b1 <- lateral_intersect_1_y - lateral_intersect_2_y
  
  cell_width <- sqrt((a1*a1)+(b1*b1))
  # All values are in pixels, need to convert to microns.
  cell_width <- cell_width / pixels_per_micron
  
  output$cell_width[i] <- cell_width
  
  return(output)
}

# Function to calculate the center position of the cell by calculating the mean of the (x,y) positions of the cell outline

centre_position <- function(df, i){
  
  centre_X <- mean(outline_data[,"V1"])
  centre_Y <- mean(outline_data[,"V2"])
  
  output$Centre_X[i] <- centre_X
  output$Centre_Y[i] <- centre_Y
  
  return(output)
  
}

# Function to calculate the offset distance (microns) of the center of the spindle and the center of the cell 

centre_offset <- function(df, i, pixels_per_micron){
  
  # Pythagoras theorem a2 + b2 = c2
  
  metaphase_x <- output$metaphase_intersect_X[i]
  metaphase_y <- output$metaphase_intersect_Y[i]
  centre_point_x <- output$Centre_X[i]
  centre_point_y <- output$Centre_Y[i]
  
  a1 <- metaphase_x - centre_point_x
  b1 <- metaphase_y - centre_point_y
  
  offset_dist <- sqrt((a1*a1)+(b1*b1))
  # All values are in pixels, need to convert to microns.
  offset_dist <- offset_dist / pixels_per_micron
  
  output$offset_dist[i] <- offset_dist
  
  return(output)
  
}

# Function to calculate the distance between the endpoints of the metaphase plate, this is the spindle width

spindle_width_func <- function(df, i, pixels_per_micron){
  
  # Pythagoras theorem a2 + b2 = c2
  
  meta_A_x <- output$metaphase_Ax[i]
  meta_A_y <- output$metaphase_Ay[i]
  meta_B_x <- output$metaphase_Bx[i]
  meta_B_y <- output$metaphase_By[i]
  
  a1 <- meta_A_x - meta_B_x
  b1 <- meta_A_y - meta_B_y
  
  spindle_width <- sqrt((a1*a1)+(b1*b1))
  # All values are in pixels, need to convert to microns.
  spindle_width <- spindle_width / pixels_per_micron
  
  output$spindle_width[i] <- spindle_width
  
  return(output)
}

# Function to calculate the spindle aspect ratio
# Spindle aspect ratio = spindle length / spindle width

spindle_aspect_ratio_func <- function(df, i){
  
  spindle_aspect_ratio <- output$centrosome_centrosome_dist[i]/output$spindle_width[i]
  output$spindle_aspect_ratio[i] <- spindle_aspect_ratio
  
  return(output)
  
}

# Function to calculate the cell aspect ratio 
# Cell aspect ratio = cell length / cell width

cell_aspect_ratio_func <- function(df, i){
  
  cell_aspect_ratio <- output$cell_length[i]/output$cell_width[i]
  output$cell_aspect_ratio[i] <- cell_aspect_ratio
  
  return(output)
}

# Function to calculate the acute angle between the spindle axis and the metaphase plate

angle_func <- function(df, i){
  
  # Calculate the gradient of the spindle axis and the metaphase plate  
  
  saX1 <- output$centrosome_Ax[i]
  saX2 <- output$centrosome_Bx[i]
  saY1 <- output$centrosome_Ay[i]
  saY2 <- output$centrosome_By[i]
  
  mpX1 <- output$metaphase_Ax[i]
  mpX2 <- output$metaphase_Bx[i]
  mpY1 <- output$metaphase_Ay[i]
  mpY2 <- output$metaphase_By[i]
  
  m1 <- (saY2-saY1)/(saX2-saX1)
  
  m2 <- (mpY2-mpY1)/(mpX2-mpX1)
  
  # absolute value to get acute angle  
  angle <- abs((m2-m1)/(1+(m2*m1)))
  # Inverse tan
  angle <- atan(angle)
  # Convert from radians to degrees
  angle <- angle * 180/pi
  
  output$spindle_angle[i] <- angle
  return(output)
  
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
  my_filename_outline <- my_files_outlines[i]
  my_filename_DNA <- my_files_DNA[i]
  # Call the function to get the data
  centrosome_data <- load_the_data(my_filename_centrosomes)
  outline_data <- load_the_data(my_filename_outline)
  DNA_data <- load_the_data(my_filename_DNA)
  # call the function to do the calculations
  output <- store_basic_info(output, centrosome_data, DNA_data, i)
  output <- find_intersections(output, i, outline_data, "centrosome")
  output <- find_intersections(output, i, outline_data, "dna")
  # these functions could be generalised further
  output <- d1_d2_func(output, i, pxpum)
  output <- l1_l2_func(output, i, pxpum)
  output <- metaphase_intersect(output, i)
  output <- D1_D2_func(output, i, pxpum)
  output <- centrosome_dist_func(output, i, pxpum)
  output <- centre_position(output, i)
  output <- centre_offset(output, i, pxpum)
  output <- spindle_width_func(output, i, pxpum)
  output <- spindle_aspect_ratio_func(output, i)
  output <- angle_func(output, i)
  output <- cell_length_func(output, i, pxpum)
  output <- cell_width_func(output, i, pxpum)
  output <- cell_aspect_ratio_func(output, i)
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

# Add the experiment number to the df
output$Experiment_number <- expt

# save the dateframe so it can be combined with other experiments in a new script
saveRDS(output, file = paste0("Output/Dataframe/", expt, ".rds"))

