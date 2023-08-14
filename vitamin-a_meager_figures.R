## Updated from monte_carlos_display_graphics_heat_map.R

### FIGURE 1: Heat map Normal-Normal Simulations ###

# load in the data with your pathname
load("/Users/hbalikci/Documents/Chicago_Harris/Life Admin/DIL/Vitamin A/vitamin-A_dropbox/meta_analysis_monte_carlo_rubin_grid_if_se_versus_sigma.RData")

# prep table row and column names
greeks=c(alpha='\u03b1', 
         tau='\u03c4', 
         sigma='\u03c3',
         beta='\u03b2',
         gamma='\u03b3')

col_names_vector <- c("SE : 0-5","5-10","10-15","15-20")
row_names_vector <- c(paste0(greeks['sigma']," : 0-5"),"5-10","10-15","15-20")

### Organise the data into a matrix of heat values ###

fe_mse_divided_by_bhm_mse <- fe_tau_error/bhm_tau_error
rownames(fe_mse_divided_by_bhm_mse) <- row_names_vector
colnames(fe_mse_divided_by_bhm_mse) <- col_names_vector

rownames(fe_tau_error) <-  row_names_vector
colnames(fe_tau_error) <- col_names_vector
rownames(bhm_tau_error) <-  row_names_vector
colnames(bhm_tau_error) <- col_names_vector

longData <- melt(fe_mse_divided_by_bhm_mse)

### define palate and make graphic ##

zp1 <- ggplot(longData,
              aes(x = Var2, y = Var1, fill = value)) +
              geom_tile(aes(fill = value), colour = "wheat") +
              scale_fill_gradient(low = "wheat", high = "red4") +
              coord_equal() +ggtitle("Ratio of FE MSE to BHM MSE in Normal-Normal Simulations") +
              theme_bw() +
              xlab("Standard errors of point estimates") +
              ylab("Standard deviation of true effects")

print(zp1)

### FIGURE 2: heat map student t ###

rm(list = ls())

# load in the data with your pathname
load("/Users/hbalikci/Documents/Chicago_Harris/Life Admin/DIL/Vitamin A/vitamin-A_dropbox/meta_analysis_monte_carlo_rubin_grid_of_se_versus_sigma_student_t_16_cells.RData")

# prep table row and column names
greeks=c(alpha='\u03b1', tau='\u03c4', sigma='\u03c3',
         beta='\u03b2',
         gamma='\u03b3')

col_names_vector <- c("SE : 0-5","5-10","10-15","15-20")
row_names_vector <- c(paste0(greeks['sigma']," : 0-5"),"5-10","10-15","15-20")

### Organise the data into a matrix of heat values ###

fe_mse_divided_by_bhm_mse_student_t <- fe_tau_error/bhm_tau_error
rownames(fe_mse_divided_by_bhm_mse_student_t) <- c("0-5", "5-10","10-15","15-20")
colnames(fe_mse_divided_by_bhm_mse_student_t) <- c("0-5","5-10","10-15","15-20")

longData <- melt(fe_mse_divided_by_bhm_mse_student_t)

### define palate and make graphic ##

zp2 <- ggplot(longData,
              aes(x = Var2, y = Var1, fill = value)) +
              geom_tile(aes(fill = value), colour = "wheat") + 
              scale_fill_gradient(low = "wheat", high = "red4") + 
              coord_equal() +
              ggtitle("Ratio of FE MSE to BHM MSE \n in Simulations with Student t distributed effects") + 
              theme_bw() +
              xlab("Standard errors of point estimates") +
              ylab("Standard deviation of true effects")

print(zp2)

### FIGURE 3: heat map with location outliers ###

rm(list = ls())

# load in the data with your pathname
load("/Users/hbalikci/Documents/Chicago_Harris/Life Admin/DIL/Vitamin A/vitamin-A_dropbox/meta_analysis_monte_carlo_rubin_grid_of_se_versus_sigma_location_outlier_16_cells.RData")

# prep table row and column names
greeks=c(alpha='\u03b1', tau='\u03c4', sigma='\u03c3',
         beta='\u03b2',
         gamma='\u03b3')

col_names_vector <- c("SE : 0-5","5-10","10-15","15-20")
row_names_vector <- c(paste0(greeks['sigma']," : 0-5"),"5-10","10-15","15-20")

### Organise the data into a matrix of heat values ###

fe_mse_divided_by_bhm_mse_location_outliers <- fe_tau_error/bhm_tau_error
rownames(fe_mse_divided_by_bhm_mse_location_outliers) <- c("0-5","5-10","10-15","15-20")
colnames(fe_mse_divided_by_bhm_mse_location_outliers) <- c("0-5","5-10","10-15","15-20")

longData <- melt(fe_mse_divided_by_bhm_mse_location_outliers)

### define palate and make graphic ##

zp3 <- ggplot(longData,
              aes(x = Var2, y = Var1, fill = value)) +
              geom_tile(aes(fill = value), colour = "wheat") +
              scale_fill_gradient(low = "wheat", high = "red4") +
              coord_equal() +
              ggtitle("Ratio of FE MSE to BHM MSE \n in Normal-Normal Simulations with Location Outliers") +
              theme_bw() +
              xlab("Standard errors of point estimates") +
              ylab("Standard deviation of true effects")
print(zp3)

### FIGURE 4: heat map with precision outliers ###

rm(list = ls())

# load in the data with your pathname
load("/Users/hbalikci/Documents/Chicago_Harris/Life Admin/DIL/Vitamin A/vitamin-A_dropbox/meta_analysis_monte_carlo_rubin_grid_of_se_versus_sigma_precision_outlier_16_cells.RData")

# prep table row and column names
greeks=c(alpha='\u03b1', tau='\u03c4', sigma='\u03c3',
         beta='\u03b2',
         gamma='\u03b3')

col_names_vector <- c("SE : 0-5","5-10","10-15","15-20")
row_names_vector <- c(paste0(greeks['sigma']," : 0-5"),"5-10","10-15","15-20")

### Organise the data into a matrix of heat values ###

fe_mse_divided_by_bhm_mse_precision_outliers <- fe_tau_error/bhm_tau_error
rownames(fe_mse_divided_by_bhm_mse_precision_outliers) <- c("0-5","5-10","10-15","15-20")
colnames(fe_mse_divided_by_bhm_mse_precision_outliers) <- c("0-5","5-10","10-15","15-20")

longData <- melt(fe_mse_divided_by_bhm_mse_precision_outliers)

### define palate and make graphic ##

zp4 <- ggplot(longData,
              aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(aes(fill = value), colour = "wheat") +
  scale_fill_gradient(low = "wheat", high = "red4") +
  coord_equal() +
  ggtitle("Ratio of FE MSE to BHM MSE \n in Normal-Normal Simulations with precision outliers") +
  theme_bw() +
  xlab("Standard errors of point estimates") +
  ylab("Standard deviation of true effects")

print(zp4)
