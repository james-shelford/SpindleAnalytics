# Spindle analytics data visualisation, written by James Shelford
# Script to load and combine the dataframes generated in the spindle_analytics_process script and produce plots

# Working directory should be 'SpindleAnalytics'

# Load the required packages
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(cowplot)

# make directory for plots if it doesn't exist (it should)
ifelse(!dir.exists("Output"), dir.create("Output"), "Folder exists already")
ifelse(!dir.exists("Output/Plots"), dir.create("Output/Plots"), "Folder exists already")

# Load the dataframes for all experiments. These are rds files found in Output/Dataframe
datadir <- "Output/Dataframe"
my_files <- list.files(datadir,pattern='*.rds',full.names = TRUE)
combined_df <- bind_rows(lapply(my_files, readRDS))
combined_df$Category <- as.factor(combined_df$Category)

# How many cells per group?
summary(combined_df$Category)

# Classify the conditions using experiment numbers
combined_df <- combined_df %>% mutate(Condition = case_when(Experiment_number == 'ELR288' | Experiment_number == 'ELR293' | Experiment_number == 'ELR301' ~ 'Condition 1',
                                                            Experiment_number == 'ELR289' | Experiment_number == 'ELR302' | Experiment_number == 'ELR303' ~ 'Condition 2', 
                                                            Experiment_number == 'ELR290' | Experiment_number == 'ELR306' | Experiment_number == 'ELR307' ~ 'Condition 3', 
                                                            Experiment_number == 'JS144' | Experiment_number == 'JS146' | Experiment_number == 'JS149' ~ 'Condition 4', TRUE ~ 'Unspecified'))

combined_df$Condition <- as.factor(combined_df$Condition)

# Add treatment group 
combined_df <- combined_df %>% mutate(Treatment = case_when(Category == 'Condition 1 ctrl' | Category == 'Condition 2 ctrl' | 
                                                              Category == 'Condition 3 ctrl' | Category == 'Condition 4 ctrl' ~ 'Control',
                                                            Category == 'Condition 1 treatment' | Category == 'Condition 2 treatment' | 
                                                              Category == 'Condition 3 treatment' | Category == 'Condition 4 treatment' ~ 'Treatment', TRUE ~ 'Unspecified'))
combined_df$Treatment <- as.factor(combined_df$Treatment)

# Set the levels for plotting in a specific order
combined_df$Condition <- factor(combined_df$Condition, levels = c('Condition 1', 'Condition 2', 'Condition 3', 'Condition 4'))
combined_df$Treatment <- factor(combined_df$Treatment, levels = c('Control', 'Treatment'))

# Remove NaN
combined_df <- na.omit(combined_df)

# Function to generate the plots
# theme_cowplot removes grey background
Make_The_Plot <- function(dataframe_to_plot, measurement_to_plot, y_label){
  the_plot <- ggplot(data = dataframe_to_plot, aes_string(x = 'Condition', y = measurement_to_plot, fill = 'Treatment')) +
    theme_cowplot() +
    geom_violin(position = position_dodge(0.9)) +
    geom_boxplot(width=0.2, outlier.shape = NA, position = position_dodge(0.9), show.legend = FALSE) +
    theme(axis.text.x = element_text(face = "plain", color = 'black', size = 8, angle = 0, hjust = 0.5, vjust = 0.5), axis.text.y = element_text(face = 'plain', color = 'black', size = 8)) +
    theme(axis.title.y = element_text(size = 9,face = 'plain', color = 'black')) +
    labs(y = y_label, x = NULL) + 
    theme(legend.title = element_blank()) + 
    theme(legend.position = "bottom") +
    ylim(0, NA)
  return(the_plot)
}

# Make the plots
spindle_length <- Make_The_Plot(combined_df, 'centrosome_centrosome_dist', 'Spindle length (microns)')
spindle_width <- Make_The_Plot(combined_df, "spindle_width", "Spindle width (microns)")
spindle_aspect_ratio <- Make_The_Plot(combined_df, "spindle_aspect_ratio", "Spindle aspect ratio")
spindle_angle <- Make_The_Plot(combined_df, "spindle_angle", "Spindle angle (degrees)")
d1_d2_shift <- Make_The_Plot(combined_df, "d1_d2_offset", "Spindle shift d2-d1 (microns)")
l1_l2_shift <- Make_The_Plot(combined_df, "l1_l2_offset", "Spindle shift l2-l1 (microns)")
D1_D2_shift <- Make_The_Plot(combined_df, "D1_D2_offset", "Spindle shift D2-D1 (microns)")
offset_distance <- Make_The_Plot(combined_df, "offset_dist", "Spindle offset distance (microns)")

# combine plots
# labels = 'AUTO' generates uppercase labels
# align = 'hv' aligns the plots horizontally and vertically
combined_plot <- plot_grid(spindle_length + theme(legend.position = "none"), spindle_width + theme(legend.position = "none"), 
                           spindle_aspect_ratio + theme(legend.position = "none"), spindle_angle + theme(legend.position = "none"),
                           d1_d2_shift + theme(legend.position = "none"), D1_D2_shift + theme(legend.position = "none"),
                                 l1_l2_shift + theme(legend.position = "none"), offset_distance + theme(legend.position = "none"),
                                 labels = 'AUTO', align = 'hv', rel_widths = c(1, 1), rel_heights = c(1,1)) + theme(aspect.ratio=1)
combined_plot

# Extract a legend from one of the plots to add back
legend <- get_legend(d1_d2_shift + guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "bottom"))

# Add the legend to the bottom of the combined plot
combined_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, .05))
combined_plot

# Save the plots
ggsave("Output/Plots/combined_plot.png", plot = combined_plot, dpi = 300, width = 250, height = 250, units = 'mm')
ggsave("Output/Plots/combined_plot.pdf", plot = combined_plot, width = 250, height = 250, units = 'mm', useDingbats = FALSE)
