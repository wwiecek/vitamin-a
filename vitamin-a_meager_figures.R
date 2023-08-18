# Script to generate figure(s) summarising simulation results
# Author: Hannah Balikci, based on code by Rachael Meager
# Date: August 2023

# Run this script only after running simulations using script XXX.R
# Running XXX.R will save the relevant outputs in simulations_results/
# and these are the only files that this script uses


library(reshape2)
library(ggplot2)
library(patchwork)
library(cowplot)
library(gridGraphics)

# prep table row outside of function
greeks=c(alpha='\u03b1', 
         tau='\u03c4', 
         sigma='\u03c3',
         beta='\u03b2',
         gamma='\u03b3')

# Define a function to generate a heatmap plot
generate_heatmap <- function(data_path, title) {
  load(data_path)
  
  # Organize the data into a matrix of heat values
  fe_mse_divided_by_bhm_mse <- fe_tau_error/bhm_tau_error
  row_names_vector <- c(paste0(greeks['sigma'], " : 0-5"), "5-10", "10-15", "15-20")
  col_names_vector <- c("SE : 0-5", "5-10", "10-15", "15-20")
  
  rownames(fe_mse_divided_by_bhm_mse) <- row_names_vector
  colnames(fe_mse_divided_by_bhm_mse) <- col_names_vector
  
  longData <- melt(fe_mse_divided_by_bhm_mse)

  # Create the heatmap plot
  ggplot(longData,
         aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(aes(fill = value), colour = "wheat") +
    scale_fill_gradient(low = "wheat", high = "red4") +
    coord_equal() +
    ggtitle(title) +
    labs(fill = "Ratio") +
    theme_bw()
}

# List of data paths and corresponding titles
data_paths <- c(
  "simulations_results/meta_analysis_monte_carlo_rubin_grid_if_se_versus_sigma.RData",
  "simulations_results/vitamin-A_dropbox/meta_analysis_monte_carlo_rubin_grid_of_se_versus_sigma_student_t_16_cells.RData",
  "simulations_results/meta_analysis_monte_carlo_rubin_grid_of_se_versus_sigma_location_outlier_16_cells.RData",
  "simulations_results/meta_analysis_monte_carlo_rubin_grid_of_se_versus_sigma_precision_outlier_16_cells.RData"
)
titles <- c(
  "Normal-Normal Simulations",
  "With Student t distributed effects",
  "With Location Outliers",
  "With Precision Outliers"
)

# Create figures patchwork without legends to build in larger figure
heatmap_plots <- lapply(seq_along(data_paths), function(i) {
  generate_heatmap(data_paths[i], titles[i]) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none")
})

# Create figures patchwork with legends to extract
heatmap_plots_with_legends <- lapply(seq_along(data_paths), function(i) {
  generate_heatmap(data_paths[i], titles[i]) +
    labs(x = NULL, y = NULL)
})

# Arrange in grid
grid <- plot_grid(plotlist = heatmap_plots, ncol = 2)

# Extract figure 4s legend
legend <- get_legend(heatmap_plots_with_legends[[4]])

final_plot <- plot_grid(grid, legend, ncol = 2, rel_widths = c(0.8, 0.2))

# Combine the combined plot with a common legend and common labels
final_plot <- final_plot +
  plot_annotation(
    title = "Ratio of FE MSE to BHM MSE in Normal-Normal Simulations",
    theme = theme(plot.title = element_text(hjust = 0.5))) +
  theme(legend.title = element_blank()) +
  labs(fill = "Ratio")

# Save the final plot
ggsave("heatmap_plots.png", plot = final_plot)
