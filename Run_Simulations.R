install.packages("patchwork")

library(timeSeries)
library(NlinTS)
library(RTransferEntropy)
library(ggplot2)
library(tidyr)
library(patchwork)


AC.TO <- read.csv("~/Documents/McGill/MATH 410/AC.TO.csv")
VFV.TO <- read.csv("~/Documents/McGill/MATH 410/VFV.TO (1).csv")

#set.seed(124)

# Helper Functions -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to generate linear data
generate_linear_data <- function(n, coef) {
  # n is the length of the vector
  
  x <- numeric(n)
  y <- numeric(n)
  z <- numeric(n)
  
  for (i in 2:n) {
    z[i] <- 0.5 * z[i - 1] + rnorm(1)
    x[i] <- 0.7 * x[i - 1] + rnorm(1)
  }
  
  for (i in 2:n) {
    y[i] <- z[i] + coef*x[i-1]
    # Decrease the linear relationship for .2, .3
  }
  
  result_list <- list(x=x, y=y)
  return(result_list)
  
}

# Function to generate non-linear data 
generate_nonlinear_data <- function(n, coef) {
  # n is the length of the vector
  
  x <- numeric(n)
  y <- numeric(n)
  z <- numeric(n)
  
  for (i in 2:n) {
    z[i] <- 0.5 * 1/exp (0.3*z[i - 1])*rnorm(1)
    x[i] <- 0.7 * 1/exp (0.2*x[i - 1])*rnorm(1)
  }
  
  for (i in 2:n) {
    y[i] <- (1/exp(z[i])) * (1/ exp(coef*x[i-1])) #.3, .2 in front of X and Z. Make x=1
  }
  
  result_list <- list(x=x, y=y)
  return(result_list)
}

#data <- generate_nonlinear_data(1000)
#plot(data$x, data$y)
#data$x[1:5]
#data$y[1:5]
#data$x[995: 1000]
#data$y[995:1000]
#max(data$y)
#max(data$x)
#min(data$y)
#min(data$x)
#which.max(data$x)
#data$x[71]
#plot(data$x[abs(data$x)<100000])
# Polynomial, sine

# Run the Granger Causality Test
granger_causality_test <- function(n, linearity, coef) {
  # n is the length of the vector
  # linear can be equal to "Linear" or "Non Linear"

    if (linearity=="Linear"){
    data = generate_linear_data(n, coef)
    x = data$x
    y = data$y
  }
  
    else if (linearity=="Non Linear"){
    data = generate_nonlinear_data(n, coef)
    x = data$x
    y = data$y
  }
  
    else if (linearity=="Real World"){
    x = VFV.TO$Close
    y = AC.TO$Close
    
  }
  
    model = causality.test(y,x,1,TRUE)
    p_value = model$pvalue
    return(p_value)
  
}

# Run the Transfer Entropy Test
transfer_entropy_test <- function(n, linearity, coef){
  # n is the length of the vector
  # linear can be equal to "Linear" or "Non Linear"
  
  if (linearity=="Linear"){
    data = generate_linear_data(n, coef)
    x = data$x
    y = data$y
    
  }
  
  else if (linearity=="Non Linear"){
    data = generate_nonlinear_data(n, coef)
    x = data$x
    y = data$y
  }
  
  else if (linearity=="Real World"){
    x = VFV.TO$Close
    y = AC.TO$Close
    
  }
  
  model <- transfer_entropy(x,y)
  p_value <- model$coef[,4][1]
  return(p_value)
  
}

# Run experiments
run_experiment <- function(n_points, n_reps, linearity, test, coef){
  #n_points is the length of each vector
  #n_reps is the numebr of times we run the test
  #linearity is "Linear" or "Non Linear"
  #test is either "Granger Causality" or "Transfer Entropy"
  
  set.seed(124)
  store_pvalues = c()

  if (test=="Granger Causality"){
    
    for (i in 1:n_reps){
      print(i)
      p_value = granger_causality_test(n=n_points, linearity=linearity, coef)
      store_pvalues = c(store_pvalues, p_value)
      
    }
  }
  
  if (test=="Transfer Entropy"){
    
    for (i in 1:n_reps){
      
      p_value = transfer_entropy_test(n=n_points, linearity=linearity, coef)
      store_pvalues = c(store_pvalues, p_value)
      
    }
  }
  
 return(store_pvalues) 
  
}

# Tests on Simulated Data -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 1. Linear data with strong linearity
GC_Linear_strong <- run_experiment(n_points=100, n_reps=100, linearity="Linear", test="Granger Causality", coef=1)
TE_Linear_strong <- run_experiment(n_points=100, n_reps=100, linearity="Linear", test="Transfer Entropy", coef=1)

# 2. Linear data with weak linearity
GC_Linear_weak <- run_experiment(n_points=100, n_reps=100, linearity="Linear", test="Granger Causality", coef=0.2)
TE_Linear_weak <- run_experiment(n_points=100, n_reps=100, linearity="Linear", test="Transfer Entropy", coef=0.2)

# 3. Linear data with no relationship 
GC_Linear_none <- run_experiment(n_points=100, n_reps=100, linearity="Linear", test="Granger Causality", coef=0)
TE_Linear_none <- run_experiment(n_points=100, n_reps=100, linearity="Linear", test="Transfer Entropy", coef=0)

#  4. Non Linear data with strong relationship
GC_NonLinear_strong <- run_experiment(n_points=100, n_reps=100, linearity="Non Linear", test="Granger Causality", coef=1)
TE_NonLinear_strong <- run_experiment(n_points=100, n_reps=100, linearity="Non Linear", test="Transfer Entropy", coef=1)

# 5. Non linear data with weak relationship
GC_NonLinear_weak <- run_experiment(n_points=100, n_reps=100, linearity="Non Linear", test="Granger Causality", coef=0.2)
TE_NonLinear_weak <- run_experiment(n_points=100, n_reps=100, linearity="Non Linear", test="Transfer Entropy", coef=0.2)

# 6. Non linear data with no relationship
GC_NonLinear_none <- run_experiment(n_points=100, n_reps=100, linearity="Non Linear", test="Granger Causality", coef=0)
TE_NonLinear_none <- run_experiment(n_points=100, n_reps=100, linearity="Non Linear", test="Transfer Entropy", coef=0)

# 7. Real world data
TE_RealWorld <- run_experiment(n_points=FALSE, n_reps=1, linearity="Real World", test="Transfer Entropy", coef=FALSE)
GC_RealWorld <- run_experiment(n_points=FALSE, n_reps=1, linearity="Real World", test="Granger Causality", coef=FALSE)


# Save data -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

data_frame <- data.frame( GC_Linear_strong = GC_Linear_strong, TE_Linear_strong = TE_Linear_strong, # 1. 
                          GC_Linear_weak = GC_Linear_weak, TE_Linear_weak = TE_Linear_weak, # 2.
                          GC_Linear_none = GC_Linear_none, TE_Linear_none = TE_Linear_none, # 3. 
                          GC_NonLinear_strong = GC_NonLinear_strong, TE_NonLinear_strong = TE_NonLinear_strong, # 4. 
                          GC_NonLinear_weak = GC_NonLinear_weak, TE_NonLinear_weak = TE_NonLinear_weak, # 5.
                          GC_NonLinear_none = GC_NonLinear_none, TE_NonLinear_none = TE_NonLinear_none) #6. 

current_directory <- getwd()
write.csv(data_frame, file = paste0(current_directory, "/new_file.csv"), row.names = FALSE)

# Data Manipulation  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

df <- read.csv("~/Documents/McGill/ECON 468/new_file.csv")
df$ID <- seq_len(nrow(df))
df_long <- pivot_longer(df, cols = -ID, names_to = "Column", values_to = "Value")

# Distribution Plots -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

personalized_plot <- function(Strength, Linearity){

  df <- df_long %>% 
    dplyr::filter(Column == paste("GC", Linearity, Strength, sep = "_") | 
             Column == paste("TE", Linearity, Strength, sep = "_"))
  
  strength_title = ifelse(grepl("strong", Strength), "Strong", ifelse(grepl("weak", Strength), "Weak", "Absent"))
  linearity_title = ifelse(grepl("NonLinear", Linearity), "Non-Linear", ifelse(grepl("Linear", Linearity), "Linear", " "))
  title = paste(strength_title, linearity_title, "Relationship")
  
  df$Column <- ifelse(grepl("GC", df$Column), "Granger Causality", "Transfer Entropy")

  plot <- ggplot(df, aes(x = Value, color = Column)) +
    geom_density() +
    scale_color_discrete(name = " ") +  # Add a legend title
    geom_vline(xintercept = 0.05, color = "red", linetype = "dotted") +
    geom_vline(xintercept = 0.01, color = "blue", linetype = "dotted") +
    theme(legend.position = c(0.85, 0.95),  # Adjust the legend position
          legend.key.size = unit(0.3, "cm"),  # Set the size of legend keys
          legend.background = element_rect(fill = "transparent", color = NA),  # Remove legend background
          legend.key = element_rect(fill = "transparent"), # Remove legend key background
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "P-Value of Test", y = " ",  title = "")
  
  filename <- paste0(title, ".jpg")
  #ggsave(filename, plot, device = "jpg")
  return(plot)
}
  
# 1. Linear data with strong linearity
personalized_plot("strong", "Linear")
# Comments: The GC test does very well, most tests reject null AT 1%
#           The TE test does badly. 

# 2. Linear data with weak linearity
personalized_plot("weak", "Linear")
# Comments: At the 5% and 1% level GC does better than TE.

# 3. Linear data with no relationship 
personalized_plot("none", "Linear")
# Comments: Both have a high number of significant p-values, GC does worse

#  4. Non Linear data with strong relationship
personalized_plot("strong", "NonLinear")
# Comments: GC does very well. They reject every time, TE doesn't do well.

# 5. Non linear data with weak relationship
personalized_plot("weak", "NonLinear")
# Comments: GC does very well. TE does not do well.

# 6. Non linear data with no relationship
personalized_plot("none", "NonLinear")
# Comments: GC does slightly worse than TE. But not by a lot. 

# Empirical CDF plot -----------------------------------------------------------------------------------

Empirical_CDF <- function(Strength, Linearity){
  
  df <- df_long %>% 
    dplyr::filter(Column == paste("GC", Linearity, Strength, sep = "_") | 
             Column == paste("TE", Linearity, Strength, sep = "_"))
    
  strength_title = ifelse(grepl("strong", Strength), "Strong", ifelse(grepl("weak", Strength), "Weak", "Absent"))
  linearity_title = ifelse(grepl("NonLinear", Linearity), "Non-Linear", ifelse(grepl("Linear", Linearity), "Linear", " "))
  title = paste(strength_title, linearity_title, "Relationship")
  
  df$Column <- ifelse(grepl("GC", df$Column), "Granger Causality", "Transfer Entropy")
    
  plot <- ggplot(df, aes(x = Value, color = Column)) +
    stat_ecdf(geom = "step") +
    scale_color_discrete(name = " ") +  # Add a legend title
    geom_vline(xintercept = 0.05, color = "red", linetype = "dotted") +
    geom_vline(xintercept = 0.01, color = "blue", linetype = "dotted") +
    theme(legend.position = c(0.85, 0.15),  # Adjust the legend position
          legend.key.size = unit(0.3, "cm"),  # Set the size of legend keys
          legend.background = element_rect(fill = "transparent", color = NA),  # Remove legend background
          legend.key = element_rect(fill = "transparent"), # Remove legend key background
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "P-Value of Test", y = "Cumulative Distribution",  title = "")
  
  return(plot)
  
}

Empirical_CDF("weak", "Linear")

# Scatter Plot -----------------------------------------------------------------------------------

Scatter_Plot <- function(Strength, Linearity){
  
  col_1 <- paste(c("TE_", Linearity, "_", Strength), collapse = "")
  col_2 <- paste(c("GC_", Linearity, "_", Strength), collapse = "")
  
  df <- df %>% dplyr::select(col_1, col_2)
  names(df) <- c("one", "two")
  
  strength_title = ifelse(grepl("strong", Strength), "Strong", ifelse(grepl("weak", Strength), "Weak", "Absent"))
  linearity_title = ifelse(grepl("NonLinear", Linearity), "Non-Linear", ifelse(grepl("Linear", Linearity), "Linear", " "))
  title = paste(strength_title, linearity_title, "Relationship,")
  
  col_1_name <- ifelse(grepl("GC", col_1), "Granger Causality", "Transfer Entropy")
  col_2_name <- ifelse(grepl("GC", col_2), "Granger Causality", "Transfer Entropy")
  
  correlation_coefficient <- cor(df$one, df$two)
  
  plot <- ggplot(df, aes(x = one, y = two)) +
    geom_point(color = "blue") +
    labs(x = col_1_name, y = col_2_name, title = paste("r =", round(correlation_coefficient, 2)), collapse = "") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)

}

Scatter_Plot("none", "Linear")

# All Graphs -----------------------------------------------------------------------------------

all_graphs <- function( Strength, Linearity) {
  
  distribution_plot <- personalized_plot(Strength, Linearity) + theme(aspect.ratio = 1)
  empirical_cdf <- Empirical_CDF(Strength, Linearity) + theme(aspect.ratio = 1)
  scatter_plot <- Scatter_Plot(Strength, Linearity) + theme(aspect.ratio = 1)
  
  strength_title = ifelse(grepl("strong", Strength), "Strong", ifelse(grepl("weak", Strength), "Weak", "Absent"))
  linearity_title = ifelse(grepl("NonLinear", Linearity), "Non-Linear", ifelse(grepl("Linear", Linearity), "Linear", " "))
  title = paste(strength_title, linearity_title, "Relationship")
  
  combined_plot <- distribution_plot + empirical_cdf + scatter_plot + plot_layout(ncol = 3) + 
    plot_annotation(title = title, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  
  ggsave(paste(title,".jpg"), combined_plot, width = 10, height = 6, dpi = 300)
  return(combined_plot)
  
}


# 1. Linear data with strong linearity
all_graphs("strong", "Linear")

# 2. Linear data with weak linearity
all_graphs("weak", "Linear")

# 3. Linear data with no relationship 
all_graphs("none", "Linear")

#  4. Non Linear data with strong relationship
all_graphs("strong", "NonLinear")

# 5. Non linear data with weak relationship
all_graphs("weak", "NonLinear")

# 6. Non linear data with no relationship
all_graphs("none", "NonLinear")



# Tables -----------------------------------------------------------------------------------
View(df)

calculate_2x2_table <- function(p, Strength, Linearity) {
  
  # Edit column names
  col_1 <- paste(c("TE_", Linearity, "_", Strength), collapse = "")
  col_2 <- paste(c("GC_", Linearity, "_", Strength), collapse = "")
  
  print(col_1)
  df <- df %>% dplyr::select(col_1, col_2)
  names(df) <- c("TE", "GC")
  
  # Count occurrences based on conditions
  A <- sum(df$TE < p & df$GC < p)
  B <- sum(df$TE > p & df$GC < p)
  C <- sum(df$TE < p & df$GC > p)
  D <- sum(df$TE > p & df$GC > p)
  
  # Calculate proportions
  total_count <- nrow(df)
  proportion_A <- A / total_count
  proportion_B <- B / total_count
  proportion_C <- C / total_count
  proportion_D <- D / total_count
  
  # Create a 2x2 table
  result_table <- matrix(c(proportion_A, proportion_B, proportion_C, proportion_D), nrow = 2, byrow = TRUE)
  rownames(result_table) <- c("GC<p", "GC>p")
  colnames(result_table) <- c("TE<p", "TE>p")
  
  # Display the table
  print(result_table)
  
}


p_1=0.05
# 1. Linear data with strong linearity
calculate_2x2_table(p_1, "strong", "Linear")

# 2. Linear data with weak linearity
calculate_2x2_table(p_1, "weak", "Linear")

# 3. Linear data with no relationship 
calculate_2x2_table(p_1, "none", "Linear")

#  4. Non Linear data with strong relationship
calculate_2x2_table(p_1, "strong", "NonLinear")

# 5. Non linear data with weak relationship
calculate_2x2_table(p_1, "weak", "NonLinear")

# 6. Non linear data with no relationship
calculate_2x2_table(p_1, "none", "NonLinear")


p_2=0.01
# 1. Linear data with strong linearity
calculate_2x2_table(p_2, "strong", "Linear")

# 2. Linear data with weak linearity
calculate_2x2_table(p_2, "weak", "Linear")

# 3. Linear data with no relationship 
calculate_2x2_table(p_2, "none", "Linear")

#  4. Non Linear data with strong relationship
calculate_2x2_table(p_2, "strong", "NonLinear")

# 5. Non linear data with weak relationship
calculate_2x2_table(p_2, "weak", "NonLinear")

# 6. Non linear data with no relationship
calculate_2x2_table(p_2, "none", "NonLinear")

