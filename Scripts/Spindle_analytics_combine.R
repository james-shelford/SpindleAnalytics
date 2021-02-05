# Spindle analytics data visualisation, written bu James Shelford
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

# Function to generate the plots

Make_The_Plot <- function(dataframe_to_plot, measurement_to_plot, y_label){
  the_plot <- ggplot(data = dataframe_to_plot, aes_string(x = 'Condition', y = measurement_to_plot, colour = 'Treatment')) +
    geom_quasirandom(alpha = 0.5, stroke=0, dodge.width = 1) + 
    stat_summary(fun.data = mean_se, geom = 'point', size = 2, aes(group = Treatment), position = position_dodge(1)) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom = 'errorbar', size = 0.8, aes(group = Treatment), position = position_dodge(1), width = 0) +
    theme(axis.text.x = element_text(face = "plain", color = 'black', size = 8, angle = 0, hjust = 0.5, vjust = 0.5), axis.text.y = element_text(face = 'plain', color = 'black', size = 8)) +
    theme(axis.title.y = element_text(size = 9,face = 'plain', color = 'black')) +
    labs(y = y_label, x = NULL) + 
    theme(legend.title = element_blank()) + 
    theme(legend.position = 'none') +
    ylim(0, NA)
  return(the_plot)
}

# Make the plots
spindle_length <- Make_The_Plot(combined_df, 'centrosome_centrosome_dist', 'Spindle length (microns)')
spindle_width <- Make_The_Plot(combined_df, "spindle_width", "Spindle width (microns)")
spindle_aspect_ratio <- Make_The_Plot(combined_df, "spindle_aspect_ratio", "Spindle aspect ratio (length/width)")
spindle_angle <- Make_The_Plot(combined_df, "spindle_angle", "Spindle angle (degrees)")

d1_d2_shift <- Make_The_Plot(combined_df, "d1_d2_offset", "Spindle shift d2-d1 (microns)")
l1_l2_shift <- Make_The_Plot(combined_df, "l1_l2_offset", "Spindle shift l2-l1 (microns)")
D1_D2_shift <- Make_The_Plot(combined_df, "D1_D2_offset", "Spindle shift D2-D1 (microns)")
offset_distance <- Make_The_Plot(combined_df, "offset_dist", "Spindle offset distance (microns)")

# combine plots
spindle_lengths_plot <- plot_grid(spindle_length, spindle_aspect_ratio)
# labels = 'AUTO' generates uppercase labels
# align = 'hv' aligns the plots horizontally and vertically
spindle_shifts_plot <- plot_grid(d1_d2_shift, D1_D2_shift, l1_l2_shift, offset_distance, labels = 'AUTO', align = 'hv', rel_widths = c(1, 1), rel_heights = c(1,1)) + theme(aspect.ratio=1)

# Save the plots
ggsave("Output/Plots/combined_spindle_length_plot.png", plot = spindle_lengths_plot, dpi = 300)
ggsave("Output/Plots/combined_spindle_length_plot.pdf", plot = spindle_lengths_plot, width = 250, height = 250, units = 'mm', useDingbats = FALSE)
ggsave("Output/Plots/combined_spindle_shift_plot.png", plot = spindle_shifts_plot, dpi = 300)
ggsave("Output/Plots/combined_spindle_shift_plot.pdf", plot = spindle_shifts_plot, width = 250, height = 250, units = 'mm', useDingbats = FALSE)

# Summary statistics
summary_table <- combined_df %>% group_by(Condition, Treatment) %>% summarise_at(vars(centrosome_centrosome_dist), list(mean_spindle_length = mean, standard_dev = sd))

# Multiple comparisons
# all_aov2 <- aov(centrosome_centrosome_dist ~ Treatment + Cell_Line, data = combined_df)
# summary(all_aov2)
all_aov3 <- aov(centrosome_centrosome_dist ~ Treatment * Condition, data = combined_df)
summary(all_aov3)
TukeyHSD(all_aov3, which = "Condition")
pairwise.t.test(combined_df$centrosome_centrosome_dist, combined_df$Condition,
                p.adjust.method = "BH")
# for unequal groups
library(car)
Anova(all_aov3, type = "III")
