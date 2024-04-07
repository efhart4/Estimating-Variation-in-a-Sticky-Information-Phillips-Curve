library(xtable) # table to latex
library(latex2exp) # use Latex in plot

# Load fixed data
#-------------------------------------------------------------------------------
load("/Fixed_Data.Rda")

# Access variables from the loaded data
alpha_samples_fixed <- Fixed_Data[[2]]
lambda_samples_fixed <- Fixed_Data[[3]]
hv_samples_fixed <- Fixed_Data[[4]]
regression_coef_fixed <- Fixed_Data[[5]]
percent_change_samples_fixed <- Fixed_Data[[6]]
output_gap_IV_fixed <- Fixed_Data[[7]]
marginal_likelihood_fixed <- Fixed_Data[[8]]
acceptance_rate_fixed <- Fixed_Data[[9]]


# Measure equation Parameters Fixed Model
#------------------------------
mean_alpha_fixed <- apply(alpha_samples_fixed, 1, mean)
sd_alpha_fixed <- apply(alpha_samples_fixed, 1, sd)
#print(paste("The mean of Alpha fixed",mean_alpha_fixed ))

mean_lambda_fixed <- apply(lambda_samples_fixed, 1, mean)
sd_lambda_fixed <- apply(lambda_samples_fixed, 1, sd)
#print(paste("The mean of lambda fixed",mean_lambda_fixed ))

sigma_samples_fixed  <- 1 / hv_samples_fixed
mean_sigma_error_fixed<- apply(sigma_samples_fixed, 1, mean)
sd_sigma_error_fixed <- apply(sigma_samples_fixed, 1, sd)
#print(paste("The mean of SD Error fixed",mean_sigma_error_fixed ))

mean_hv_fixed <- apply(hv_samples_fixed, 1, mean)
sd_hv_fixed <- apply(hv_samples_fixed, 1, sd)
#print(paste("The mean of lambda fixed",mean_lambda_fixed ))

# put fixed estimates into a table
Measure_table_fixed <- matrix(c(mean_alpha_fixed, mean_lambda_fixed, mean_hv_fixed, sd_alpha_fixed, sd_lambda_fixed, sd_hv_fixed), nrow = 2, ncol = 3, byrow = TRUE)
colnames(Measure_table_fixed) <- c("Alpha", "Lambda", "HV")
rownames(Measure_table_fixed) <- c("Mean", "SD")
#print(xtable(Measure_table_fixed, caption = "Measure Parameters in Fixed Model", label = "tab:mean_sd")) # nolint: line_length_linter.




# IV Parameters in Fixed Model
#------------------------------
# Calculating mean and standard deviation
mean_output_gap_IV_fixed <- apply(output_gap_IV_fixed, 1, mean)
sd_output_gap_IV_fixed <- apply(output_gap_IV_fixed, 1, sd)

# Creating the table
IV_table_fixed <- matrix(c(mean_output_gap_IV_fixed, sd_output_gap_IV_fixed), nrow = 2, ncol = 5, byrow = TRUE)
colnames(IV_table_fixed) <- c("Parameter 1", "Parameter 2", "Parameter 3", "Parameter 4", "Parameter 5")
rownames(IV_table_fixed) <- c("Mean", "SD")

# Printing the table in LaTeX format
# print(xtable(IV_table_fixed, caption = "IV Parameters in Fixed Model", label = "tab:mean_sd"))







# Load TVP data
#-------------------------------------------------------------------------------
load("~/TVP_Data.Rda")

# Access variables from the loaded data
alpha_samples_tvp <- TVP_Data[[2]]
lambda0_samples_tvp <- TVP_Data[[3]]
hv_samples_tvp <- TVP_Data[[4]]
hw_samples_tvp <- TVP_Data[[5]]
covariance_samples_tvp <- TVP_Data[[6]]
he_samples_tvp <- TVP_Data[[7]]
regression_coef_tvp <- TVP_Data[[8]]
lambda_path_samples_tvp <- TVP_Data[[9]]
percent_change_samples_tvp <- TVP_Data[[10]]
output_gap_IV_tvp <- TVP_Data[[11]]
marginal_likelihood_tvp <- TVP_Data[[12]]
acceptance_rate_TVP <- TVP_Data[[13]]

# Measure equation Parameters TVP Model
#------------------------------
mean_alpha_tvp <- apply(alpha_samples_tvp, 1, mean)
sd_alpha_tvp <- apply(alpha_samples_tvp, 1, sd)
#print(paste("The mean of Alpha TVP",mean_alpha_tvp ))

mean_lambda0_tvp <- apply(lambda0_samples_tvp, 1, mean)
sd_lambda0_tvp <- apply(lambda0_samples_tvp, 1, sd)
#print(paste("The mean of lambda0 TVP",mean_lambda0_tvp ))

mean_hv_tvp <- apply(hv_samples_tvp, 1, mean)
sd_hv_tvp <- apply(hv_samples_tvp, 1, sd)
#print(paste("The mean of lambda0 TVP",mean_lambda0_tvp ))

mean_hw_tvp <- apply(hw_samples_tvp, 1, mean)
sd_hw_tvp <- apply(hw_samples_tvp, 1, sd)
#print(paste("The mean of lambda0 TVP",mean_lambda0_tvp ))

mean_he_tvp <- apply(he_samples_tvp, 1, mean)
sd_he_tvp <- apply(he_samples_tvp, 1, sd)
#print(paste("The mean of lambda0 TVP",mean_lambda0_tvp ))

mean_covariance_tvp <- apply(covariance_samples_tvp, 1, mean)
sd_covariance_tvp <- apply(covariance_samples_tvp, 1, sd)
#print(paste("The mean of lambda0 TVP",mean_lambda0_tvp ))

# Create table
measure_table_tvp <- matrix(
  c(mean_alpha_tvp, mean_lambda0_tvp, mean_hv_tvp, mean_hw_tvp, mean_he_tvp, mean_covariance_tvp,
    sd_alpha_tvp, sd_lambda0_tvp, sd_hv_tvp, sd_hw_tvp, sd_he_tvp, sd_covariance_tvp),
  nrow = 2, ncol = 6, byrow = TRUE
)

# Assign row and column names
rownames(measure_table_tvp) <- c("Mean", "SD")
colnames(measure_table_tvp) <- c("Alpha", "Lambda0", "HV", "HW", "HE", "Covariance")

# Print table
#print(xtable(measure_table_tvp, caption = "Measure Parameters in TVP Model", label = "tab:mean_sd")) # nolint: line_length_linter.









# IV Parameters in TVP Model
#------------------------------
# Calculating mean and standard deviation
mean_output_gap_IV_tvp <- apply(output_gap_IV_tvp, 1, mean)
sd_output_gap_IV_tvp <- apply(output_gap_IV_tvp, 1, sd)

# Creating the table
IV_table_tvp <- matrix(c(mean_output_gap_IV_tvp, sd_output_gap_IV_tvp), nrow = 2, ncol = 5, byrow = TRUE)
colnames(IV_table_tvp) <- c("Parameter 1", "Parameter 2", "Parameter 3", "Parameter 4", "Parameter 5")
rownames(IV_table_tvp) <- c("Mean", "SD")

# Printing the table in LaTeX format
#print(xtable(IV_table_tvp, caption = "IV Parameters in TVP Model", label = "tab:mean_sd"))




# Graph attention
#------------------------------
n <- NCOL(lambda_path_samples_tvp)
minues_lambda <- 1 - lambda_path_samples_tvp  
lambda_filtered <- apply(minues_lambda, 1, mean, na.rm = TRUE)
per_period_std <- apply(minues_lambda, 1, sd, na.rm = TRUE)
#print(mean(lambda_filtered))
num_stds <- 1.96
per_period_std <- sd(lambda_filtered, na.rm = TRUE) / sqrt(n)
upper_bound <- lambda_filtered + num_stds * per_period_std
lower_bound <- lambda_filtered - num_stds * per_period_std

plot_matrix <- cbind(lower_bound, lambda_filtered, upper_bound)

png("/Users/modelt/Documents/Research/Estimating Shifts in Attention /project steps 2/step 9/attentionbounds2.png", width = 1200, height = 600)

matplot(Fixed_Data[[10]][1:186], plot_matrix[1:186, ],
        type = "l", xaxt = "n", xlab = "", ylab = "1 - lambda", 
        cex.main = 1.5, col = c("red","blue", "red"), lty = c(2, 1, 2))

axis(1, Fixed_Data[[10]], format(Fixed_Data[[10]], "%b %Y"), cex.axis = 1)
legend(x = "topright",legend=c(TeX(r'($1 - \bar{\lambda_t} \pm 1.96\frac{S_t}{\sqrt{n_t}$})'),TeX(r'($1 - \bar{\lambda_t}$)')),
       fill = c("red","blue"), cex=1)
dev.off()


# Analyzing graph
#------------------------------
print(cbind("Timing" = head(format(Fixed_Data[[10]], "%b %Y"), 5), "Mean" = head(lambda_filtered, 5)))
print(cbind("Timing" = tail(format(Fixed_Data[[10]], "%b %Y"), 5), "Mean" = tail(lambda_filtered, 5)))

# Find min and max values along with corresponding timing
min_lambda <- min(lambda_filtered[1:186])
max_lambda <- max(lambda_filtered[1:186])
min_timing <- Fixed_Data[[10]][which.min(lambda_filtered[1:186])]
max_timing <- Fixed_Data[[10]][which.max(lambda_filtered[1:186])]

# Print min and max values along with timing
print(paste("Minimum value of 1 - lambda:", min_lambda, "at timing:", min_timing))
print(paste("Maximum value of 1 - lambda:", max_lambda, "at timing:", max_timing))

# find the mean of the TVP path
mean_TVP = mean(lambda_filtered[1:186])
print(paste("Mean value of 1 - lambda:", mean_TVP))

# Find range of data and the first and last estimate of lambda 
# Print first and last value of lambda

print(paste("First value of 1 - lambda:", lambda_filtered[1]))
print(paste("Last value of 1 - lambda:", lambda_filtered[186]))

# Print last and first date
print(paste("First date:", format(Fixed_Data[[10]][186], "%b %Y")))
print(paste("Last date:", format(Fixed_Data[[10]][1], "%b %Y")))



# Using baseline calibration in Branch et. al (2009)
lambda1= lambda_filtered[1]
lambda2= lambda_filtered[186]
lambda_lit =  0.25

coefficient = function(i, lambda){
  rho = 0.8
  omega = 1
  alpha = 0.1
  help = (1 - lambda)^(i+1)
  
  numerator = rho^i
  denominator_addend_1 = omega * alpha^2
  denominator_addend_2 = help/ (1 - help)
  denominator = denominator_addend_1 + denominator_addend_2
  
  coeff = numerator / denominator
  coeffsqr = (coeff)^2
  return(coeffsqr)
}



VAR_price = function(lambda){

  varE = 0.1
  
  N= 1000
  sum = 0
  for (i in 1:N){
    sum = sum + coefficient(i, lambda)
  }
  var = varE * sum
  return(var)
}


VAR_output = function(lambda){
  
  varE = 1
  
  N= 100
  sum = 0
  for (i in 1:N){
    sum = sum + coefficient(i, lambda)
  }
  var = varE * sum
  return(var)
}




start_price_variance = VAR_price(lambda1)
end_price_variance = VAR_price(lambda2)
literature_price_variance = VAR_price(lambda_lit)

print(paste("Implied starting variance of price:", start_price_variance))
print(paste("Implied ending variance of price:", end_price_variance))

# percent change 
percent_change = 100 * (end_price_variance - start_price_variance)/start_price_variance
print(paste("Implied percent change in variance of price:", percent_change))

print(paste("literature variance of price:", literature_price_variance))

# Bayesian Model Comparison
#---------------------------
# Posterior Odds Ratio
ratio = marginal_likelihood_fixed/ marginal_likelihood_tvp

print(paste("The posterior odds ration is:", ratio))

# Probability of TVP Model:
tvp_p = (marginal_likelihood_tvp)/( marginal_likelihood_fixed + marginal_likelihood_tvp)
print(paste("The probability of TVP model is:", tvp_p))





# Printing Tables
#-----------------
print(xtable(Measure_table_fixed, caption = "Measure Parameters in Fixed Model", label = "tab:mean_sd")) # nolint: line_length_linter.
print(xtable(IV_table_fixed, caption = "IV Parameters in Fixed Model", label = "tab:mean_sd"))

print(xtable(measure_table_tvp, caption = "Measure Parameters in TVP Model", label = "tab:mean_sd")) # nolint: line_length_linter.
print(xtable(IV_table_tvp, caption = "IV Parameters in TVP Model", label = "tab:mean_sd"))

