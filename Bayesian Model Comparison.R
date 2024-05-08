library(MASS)  
library(Brobdingnag) # for very large numbers


# Marginal Likelihood
marginal_likelihood_TVP= function( inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z, alpha_samples_tvp, lambda0_samples_tvp, hv_samples_tvp, hw_samples_tvp, covariance_samples_tvp, he_samples_tvp, output_gap_IV_tvp){

  set.seed(24689)
  N = nrow(inflation)
  J = 5


  #Prior Parameters
    # Hyper parameters for Normal distribution describing prior for conditional mean parameters: 
    # lambda0, alpha, covariance and a delta vector for Instrumental Equation
    k = 3 + 5  
    mu = matrix(0, k, 1)
    mu[1] = 0.1 # value for alpha in theory. 
    mu[2] = 0.7 
    mu[3] = 0 # covariance between measure error and instrument error
    mu[4:k] = 0

  
  # Variance of prior 
  # I looked at a normal graph these seem reasonable. 
    V= diag(k) 
    V[1,1]= 0.01^2 # ALPHA
    V[2,2]= 1  # LAMBDA 
  

  # set variance of instrument equation high to create a diffuse prior
    for (i in 3:k){V[i,i]= 1}


  
  # Hyper parameters from Gamma distribution describing prior.E(h)= alpha*beta, variance: alpha*beta^2
    # measure disturbance
    alpha_v=100
    beta_v = 4
    
    # random walk
    alpha_w=80
    beta_w = 1
    
    # Instrumental Equation
    alpha_e = 45
    beta_e = 2
  
  
  # NUMBER OF SIMULATIONS FROM THE MH SAMPLER
    G = 2000;       #NUMBER OF POST CONVERGENCE DRAWS
    
  
  # Sampler
    num_param = k+3
    mu_prop = matrix(0, num_param, 1)
  
  
  # # Variance covariance matrix of shock to random walk proposal
    R =  diag(num_param)
  # Measure Parameters
    R[1,1] = 5e-6 # alpha
    R[2,2] = 1e-4# lambda0
    R[3,3] = 1e-13# covariance

  
  # Instrumental Parameters
    R[4,4] = 1e-7 # og_intercept
    R[5,5] = 1e-5 # og_lag1
    R[6,6] = 1e-5 # og_lag2
    R[7,7] = 1e-6 # og_lag3
    R[8,8] = 1e-5# og_inflation_expectation

  # Gamma Parameters 
    R[9,9] = 1.5 #hv # Measure Error
    R[10,10] = 1.1 #hw # Random Walk
    R[11,11] = 1.1  # he # Instrumental Error
    
    R_scale = 5 # 3 final acceptance 0.44 # 1 final acceptance 0.56 #0.15 final acceptance 0.672 #0.09  final acceptance 0.7
    R = R_scale*R
    #print(paste0("R_scale: ", R_scale))




    
    
    
    
  # Theta Tilde
    theta_tilde = matrix(0,num_param,1)  
  # # Compute posterior medians at which to compute basic marginal likelihood identity
    median_alpha_tvp <- apply(alpha_samples_tvp, 1, median)
    median_lambda0_tvp <- apply(lambda0_samples_tvp, 1, median)
    median_hv_tvp <- apply(hv_samples_tvp, 1, median)
    median_hw_tvp <- apply(hw_samples_tvp, 1, median)
    median_he_tvp <- apply(he_samples_tvp, 1, median)
    median_covariance_tvp <- apply(covariance_samples_tvp, 1, median)
    median_output_gap_IV_tvp <- apply(output_gap_IV_tvp, 1, median)  
  # # Measure Parameters
    theta_tilde[1] = median_alpha_tvp
    theta_tilde[2] = median_lambda0_tvp
    theta_tilde[3] = median_covariance_tvp
  # # Instrumental Parameters
    theta_tilde[4:8] = median_output_gap_IV_tvp
  # # Gamma Parameters  
    theta_tilde[9] = median_hv_tvp #hv # Measure Error
    theta_tilde[10] = median_hw_tvp #hw # Random Walk
    theta_tilde[11] = median_he_tvp # he # Instrumental Error

  # # Evaluate likelihood function and prior at theta_tilde
    particle_list_tilde= IV_Particle_Filter(theta_tilde,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex)
    log_lik_theta_tilde = particle_list_tilde[[1]]
    log_prior_theta_tilde = Log_MVN_pdf_kernel(theta_tilde[1:8],mu,V)+ Log_GAMMA_pdf_kernel(theta_tilde[9],alpha_v,beta_v) + Log_GAMMA_pdf_kernel(theta_tilde[10],alpha_w,beta_w)+ Log_GAMMA_pdf_kernel(theta_tilde[11],alpha_e,beta_e)
    
    
  # # Evaluate likelihood function and prior at theta_tilde (for final calculation of marginal likelihood)
    lik_val_theta_tilde = exp(as.brob(log_lik_theta_tilde))
    # you have an issue here, lik_val_theta_tilde is infinity
    prior_val_theta_tilde = exp( Log_MVN_pdf(theta_tilde[1:8],mu,V)+ Log_GAMMA_pdf(theta_tilde[9],alpha_v,beta_v) + Log_GAMMA_pdf(theta_tilde[10],alpha_w,beta_w)+ Log_GAMMA_pdf(theta_tilde[11],alpha_e,beta_e) )
    
    
    
    # Experimental addition
    log_prior_val_theta_tilde = Log_MVN_pdf(theta_tilde[1:8],mu,V)+ Log_GAMMA_pdf(theta_tilde[9],alpha_v,beta_v) + Log_GAMMA_pdf(theta_tilde[10],alpha_w,beta_w)+ Log_GAMMA_pdf(theta_tilde[11],alpha_e,beta_e) 
    
    
    
    
    if(is.na(log_lik_theta_tilde))
    {print("data has zero probability at theta tilde")
      MedianData = list()
      MedianData[[1]]= FALSE
      return(MedianData)}
  
    
  # Theta g
  theta_g = matrix(0,num_param,1)  

  #Storage 
  func1_save = 0
  func2_save = 0
    
    for (itr in 1:G) {
      # # Measure Parameters
      theta_g[1] = alpha_samples_tvp[itr]
      theta_g[2] = lambda0_samples_tvp[itr]
      theta_g[3] = covariance_samples_tvp[itr]
      # # Instrumental Parameters
      theta_g[4:8] = output_gap_IV_tvp[,itr]
      # # Gamma Parameters  
      theta_g[9] = hv_samples_tvp[itr] 
      theta_g[10] = hw_samples_tvp[itr]
      theta_g[11] = he_samples_tvp[itr] 
      
      # # Evaluate likelihood function and prior at theta_g (our current sample)
      particle_list_g= IV_Particle_Filter(theta_g,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex)
      log_lik_theta_g = particle_list_g[[1]]
      log_prior_theta_g = Log_MVN_pdf_kernel(theta_g[1:8],mu,V)+ Log_GAMMA_pdf_kernel(theta_g[9],alpha_v,beta_v) + Log_GAMMA_pdf_kernel(theta_g[10],alpha_w,beta_w)+ Log_GAMMA_pdf_kernel(theta_g[11],alpha_e,beta_e)
      
      
      # compare the sample to the median of all samples
      log_acceptance_prob = (log_lik_theta_tilde + log_prior_theta_tilde) - (log_prior_theta_g + log_lik_theta_g)
      acceptance_prob = exp(log_acceptance_prob)

      acceptance_prob = min(acceptance_prob, 1)
      
      

      # this is taking the difference between the median and the current sample, then finding the likelihood of this difference given the R matrix used for the random walk in MH sampler. 
      propose_theta_tilde_g = exp(Log_MVN_pdf(abs(theta_tilde - theta_g), mu_prop, R)) 
      func1_save = acceptance_prob * propose_theta_tilde_g + func1_save
      
      
      # Part 2
      #------------------------------------------------
      theta_star = theta_tilde + mvrnorm(n = 1, mu_prop,R )
    
      particle_list= IV_Particle_Filter(theta_star,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex)
      log_lik_theta_star = particle_list[[1]]

      

      # If this is the case, the probability of accepting the proposed parameters should be 0. 
      if(is.na(log_lik_theta_star)){
        acceptance_prob2 =0
        # print("theta_star has 0 chance of acceptance")
      }else{
        # # likelihood of the parameters
        log_prior_theta_star  = Log_MVN_pdf_kernel(theta_star[1:8],mu,V) + Log_GAMMA_pdf_kernel(theta_star[9],alpha_v,beta_v)+ Log_GAMMA_pdf_kernel(theta_star[10],alpha_w,beta_w)+ Log_GAMMA_pdf_kernel(theta_star[11],alpha_e,beta_e)
        
        # # collect the pieces and find the final probability of acceptance. 
        log_acceptance_prob2 = (log_lik_theta_star+log_prior_theta_star)-(log_lik_theta_tilde + log_prior_theta_tilde)
        acceptance_prob2 = exp(log_acceptance_prob2)
        acceptance_prob2 = min(acceptance_prob2,1)
      }# end else 
      
      func2_save = func2_save + acceptance_prob2;
      
      if(itr%%1000 ==0)
      {
        print(paste("iteration is: ", itr))
      }
      
    
    }# endloop                                                                                       
    
    func1_save = func1_save/G
    
    func2_save = func2_save/G
    
    post_val_theta_tilde = func1_save/func2_save
    
      
    marg_lik_cj = (lik_val_theta_tilde*prior_val_theta_tilde) / post_val_theta_tilde
      
      
      
    # OUTPUT SAMPLES
    SummaryData = list()
    SummaryData[[1]]= TRUE
    SummaryData[[2]]= marg_lik_cj
    return(SummaryData)
    
}# end main function



# Marginal Likelihood
marginal_likelihood_Fixed= function( inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z, alpha_samples_fixed, lambda0_samples_fixed, hv_samples_fixed, covariance_samples_fixed, he_samples_fixed, output_gap_IV_fixed){
  
  set.seed(24689)
  N = nrow(inflation)
  J = 5
  
  # reduce data to that which is used in the TVP estimation from J+1 to N-1
  inflation = inflation[(J+1): (N-1)]
  output_gap = output_gap[(J+1):( N-1)]
  delta_natural_output = delta_natural_output[(J+1): (N-1)]
  inflation_ex = inflation_ex[(J+1):(N-1),]
  output_change_ex = output_change_ex[(J+1):(N-1),]
  Z = Z[(J+1):(N-1),]
  
  
  #Prior Parameters
  # Hyper parameters for Normal distribution describing prior for conditional mean parameters: 
  # lambda0, alpha, covariance and a delta vector for Instrumental Equation
  k = 3 + 5  
  mu = matrix(0, k, 1)
  mu[1] = 0.1 # value for alpha in theory. 
  mu[2] = 0.7 
  mu[3] = 0 # covariance between measure error and instrument error
  mu[4:k] = 0
  
  
  # Variance of prior 
  # I looked at a normal graph these seem reasonable. 
  V= diag(k) 
  V[1,1]= 0.01^2 # ALPHA
  V[2,2]= 1  # LAMBDA 
  
  
  # set variance of instrument equation high to create a diffuse prior
  for (i in 3:k){V[i,i]= 1}
  
  
  
  # Hyper parameters from Gamma distribution describing prior.E(h)= alpha*beta, variance: alpha*beta^2
  # measure disturbance
  alpha_v=100
  beta_v = 4
  
  # random walk
  #alpha_w=80
  #beta_w = 1
  
  # Instrumental Equation
  alpha_e = 45
  beta_e = 2
  
  
  # NUMBER OF SIMULATIONS FROM THE MH SAMPLER
  G = 2000;       #NUMBER OF POST CONVERGENCE DRAWS
  
  
  # Sampler
  num_param = k+2
  mu_prop = matrix(0, num_param, 1)
  
  
  # # Variance covariance matrix of shock to random walk proposal
  R =  diag(num_param)
  # Measure Parameters
  R[1,1] = 5e-6 # alpha
  R[2,2] = 1e-4# lambda0
  R[3,3] = 1e-13# covariance
  
  
  # Instrumental Parameters
  R[4,4] = 1e-7 # og_intercept
  R[5,5] = 1e-5 # og_lag1
  R[6,6] = 1e-5 # og_lag2
  R[7,7] = 1e-6 # og_lag3
  R[8,8] = 1e-5# og_inflation_expectation
  
  # Gamma Parameters 
  R[9,9] = 1.5 #hv # Measure Error
  #R[10,10] = 1.1 #hw # Random Walk
  R[10,10] = 1.1  # he # Instrumental Error
  
  R_scale = 10
  R = R_scale*R
  #print(paste0("R_scale: ", R_scale))
  
  
  
  
  
  
  
  
  # Theta Tilde
  theta_tilde = matrix(0,num_param,1)  
  # # Compute posterior medians at which to compute basic marginal likelihood identity
  median_alpha_fixed <- apply(alpha_samples_fixed, 1, median)
  median_lambda0_fixed <- apply(lambda0_samples_fixed, 1, median)
  median_hv_fixed <- apply(hv_samples_fixed, 1, median)
  #median_hw_fixed <- apply(hw_samples_fixed, 1, median)
  median_he_fixed <- apply(he_samples_fixed, 1, median)
  median_covariance_fixed <- apply(covariance_samples_fixed, 1, median)
  median_output_gap_IV_fixed <- apply(output_gap_IV_fixed, 1, median)  
  # # Measure Parameters
  theta_tilde[1] = median_alpha_fixed
  theta_tilde[2] = median_lambda0_fixed
  theta_tilde[3] = median_covariance_fixed
  # # Instrumental Parameters
  theta_tilde[4:8] = median_output_gap_IV_fixed
  # # Gamma Parameters  
  theta_tilde[9] = median_hv_fixed #hv # Measure Error
  # theta_tilde[10] = median_hw_fixed #hw # Random Walk
  theta_tilde[10] = median_he_fixed # he # Instrumental Error
  
  
  
  # # Evaluate likelihood function and prior at theta_tilde
  #particle_list_tilde= IV_Particle_Filter(theta_tilde,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex)
  log_lik_theta_tilde = likelihood_fixed_lambda(theta_tilde,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z)
  log_prior_theta_tilde = Log_MVN_pdf_kernel(theta_tilde[1:8],mu,V)+ Log_GAMMA_pdf_kernel(theta_tilde[9],alpha_v,beta_v) + Log_GAMMA_pdf_kernel(theta_tilde[10],alpha_e,beta_e)
  
  
  # # Evaluate likelihood function and prior at theta_tilde (for final calculation of marginal likelihood)
  lik_val_theta_tilde = exp(as.brob(log_lik_theta_tilde))
  # you have an issue here, lik_val_theta_tilde is infinity
  prior_val_theta_tilde = exp( Log_MVN_pdf(theta_tilde[1:8],mu,V)+ Log_GAMMA_pdf(theta_tilde[9],alpha_v,beta_v) + Log_GAMMA_pdf(theta_tilde[10],alpha_e,beta_e) )
  
  
  
  # Experimental addition
  log_prior_val_theta_tilde = Log_MVN_pdf(theta_tilde[1:8],mu,V)+ Log_GAMMA_pdf(theta_tilde[9],alpha_v,beta_v) + Log_GAMMA_pdf(theta_tilde[10],alpha_e,beta_e) 
  
  
  
  
  if(is.na(log_lik_theta_tilde))
  {print("data has zero probability at theta tilde")
    MedianData = list()
    MedianData[[1]]= FALSE
    return(MedianData)}
  
  
  # Theta g
  theta_g = matrix(0,num_param,1)  
  
  #Storage 
  func1_save = 0
  func2_save = 0
  
  for (itr in 1:G) {
    # # Measure Parameters
    theta_g[1] = alpha_samples_fixed[itr]
    theta_g[2] = lambda0_samples_fixed[itr]
    theta_g[3] = covariance_samples_fixed[itr]
    # # Instrumental Parameters
    theta_g[4:8] = output_gap_IV_fixed[,itr]
    # # Gamma Parameters  
    theta_g[9] = hv_samples_fixed[itr] 
    #theta_g[10] = hw_samples_fixed[itr]
    theta_g[10] = he_samples_fixed[itr] 
    
    # # Evaluate likelihood function and prior at theta_g (our current sample)
    # particle_list_g= IV_Particle_Filter(theta_g,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex)
    log_lik_theta_g = likelihood_fixed_lambda(theta_g,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z)
    log_prior_theta_g = Log_MVN_pdf_kernel(theta_g[1:8],mu,V)+ Log_GAMMA_pdf_kernel(theta_g[9],alpha_v,beta_v) + Log_GAMMA_pdf_kernel(theta_g[10],alpha_e,beta_e)
    
    
    # compare the sample to the median of all samples
    log_acceptance_prob = (log_lik_theta_tilde + log_prior_theta_tilde) - (log_prior_theta_g + log_lik_theta_g)
    acceptance_prob = exp(log_acceptance_prob)
    
    acceptance_prob = min(acceptance_prob, 1)
    
    
    
    # this is taking the difference between the median and the current sample, then finding the likelihood of this difference given the R matrix used for the random walk in MH sampler. 
    propose_theta_tilde_g = exp(Log_MVN_pdf(abs(theta_tilde - theta_g), mu_prop, R)) 
    func1_save = acceptance_prob * propose_theta_tilde_g + func1_save
    
    
    # Part 2
    #------------------------------------------------
    theta_star = theta_tilde + mvrnorm(n = 1, mu_prop,R )
    
    #particle_list= IV_Particle_Filter(theta_star,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex)
    log_lik_theta_star =likelihood_fixed_lambda(theta_star,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z)
    
    
    
    # If this is the case, the probability of accepting the proposed parameters should be 0. 
    if(is.na(log_lik_theta_star)){
      acceptance_prob2 =0
      # print("theta_star has 0 chance of acceptance")
    }else{
      # # likelihood of the parameters
      log_prior_theta_star  = Log_MVN_pdf_kernel(theta_star[1:8],mu,V) + Log_GAMMA_pdf_kernel(theta_star[9],alpha_v,beta_v)+  Log_GAMMA_pdf_kernel(theta_star[10],alpha_e,beta_e)
      
      # # collect the pieces and find the final probability of acceptance. 
      log_acceptance_prob2 = (log_lik_theta_star+log_prior_theta_star)-(log_lik_theta_tilde + log_prior_theta_tilde)
      acceptance_prob2 = exp(log_acceptance_prob2)
      acceptance_prob2 = min(acceptance_prob2,1)
    }# end else 
    
    func2_save = func2_save + acceptance_prob2;
    
    if(itr%%1000 ==0)
    {
      print(paste("iteration is: ", itr))
    }
    
    
  }# endloop                                                                                       
  
  func1_save = func1_save/G
  
  func2_save = func2_save/G
  
  post_val_theta_tilde = func1_save/func2_save
  
  
  marg_lik_cj = (lik_val_theta_tilde*prior_val_theta_tilde) / post_val_theta_tilde
  
  
  
  # OUTPUT SAMPLES
  SummaryData = list()
  SummaryData[[1]]= TRUE
  SummaryData[[2]]= marg_lik_cj
  return(SummaryData)
  
}# end main function




likelihood_fixed_lambda = function(theta_in,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z){
  
  Sample = NROW(inflation)
  # Theta in 
  alpha_in   = theta_in[1]  
  lambda_in = theta_in[2] 
  covariance_in = theta_in[3] 
  delta_in = theta_in[4:8]
  
  # Measure Error   
  hv_in = theta_in[9]  
  sig_v_in = hv_in^(-1)
  
  
  # Instrumental Error
  he_in  = theta_in[10]
  sig_e_in = he_in^(-1)
  
  #Variance-Covariance Matrix between error terms 
  Sigma_error_in = c(sig_v_in^2,covariance_in, covariance_in,sig_e_in^2 )
  dim(Sigma_error_in) = c(2,2)
  my_eigen = eigen(Sigma_error_in, symmetric = TRUE)
  det_Sigma_error_in = det(Sigma_error_in)
  solve_Sigma_error_in = solve(Sigma_error_in)
  mu = c(0,0)
  
  # Variance Covariance Matrix must be PSD. PSD if both eigen values are greater than zero. 
  if(my_eigen$values[1]<=0 || my_eigen$values[2] <=0){
    # print("Sigma_error_in")
    # print(Sigma_error_in)
    # print("Eigen Value")
    # print(my_eigen$values)
    print("Not PSD Matrix" )
    return(NA)
  }
  
  
  
  # Expectation of optimal price weighted by the percent of firms using that info set
  beliefs_less1 = (1-lambda_in)*( inflation_ex[,1]+ alpha_in*(output_change_ex[,1] -delta_natural_output) )
  beliefs_less2 = (1-lambda_in)*(lambda_in^1)*( inflation_ex[,2]+ alpha_in*(output_change_ex[,2] -delta_natural_output) )
  beliefs_less3 = (1-lambda_in)*(lambda_in^2)*( inflation_ex[,3]+ alpha_in*(output_change_ex[,3] -delta_natural_output) )
  beliefs_less4 = (1-lambda_in)*(lambda_in^3)*( inflation_ex[,4]+ alpha_in*(output_change_ex[,4] -delta_natural_output) )
  beliefs_less5 = (1-lambda_in)*(lambda_in^4)*( inflation_ex[,5]+ alpha_in*(output_change_ex[,5] -delta_natural_output) )
  
  # vector of exogenous variables in instrument equation
  
  # Simultaneous Equation
  output_gap_prediction = delta_in %*% t(Z) 
  output_gap_prediction = t(output_gap_prediction)
  
  # Measurement Equation
  prediction_vec =  ((1-lambda_in)/lambda_in)*alpha_in* output_gap + beliefs_less1 +  beliefs_less2 + beliefs_less3 +beliefs_less4+beliefs_less5
  
  # Measure 
  error_vec = inflation- prediction_vec   
  epsilon_vec = output_gap-output_gap_prediction
  
  
  # Combine error_vec and replicated instrument_error_t into a matrix
  X <- cbind(error_vec, epsilon_vec)
  
  # Compute the PDF for each row in the matrix
  # pdf_values <- apply(X, 1, function(x) dmvnorm(x, mean = c(0,0), sigma = Sigma_error_in))
  
  const <- 1 / (sqrt((2 * pi) ^ length(mu) * det_Sigma_error_in))
  exponent <- exp(-0.5 * rowSums((X - mu) %*% solve_Sigma_error_in * (X - mu)))
  pdf_values <- const * exponent
  
  log_pdf_values <- log(pdf_values)
  log_likelihood = sum(log_pdf_values)
  
  return(log_likelihood)    
}



IV_Particle_Filter = function(theta_in, inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex){
  
  
  Sample = NROW(inflation)
  # Theta in 
  alpha_in   = theta_in[1]  
  lambda0_in = theta_in[2] 
  covariance_in = theta_in[3] 
  delta_in = theta_in[4:8]
  
  # Measure Error   
  hv_in = theta_in[9]  
  sig_v_in = hv_in^(-1)
  
  # Stochastic Process of Random Walk
  hw_in = theta_in[10]  
  sig_w_in = hw_in^(-1)
  
  # Instrumental Error
  he_in  = theta_in[11]
  sig_e_in = he_in^(-1)
  
  #Variance-Covariance Matrix between error terms 
  Sigma_error_in = c(sig_v_in^2,covariance_in, covariance_in,sig_e_in^2 )
  dim(Sigma_error_in) = c(2,2)
  my_eigen = eigen(Sigma_error_in, symmetric = TRUE)
  det_Sigma_error_in = det(Sigma_error_in)
  solve_Sigma_error_in = solve(Sigma_error_in)
  mu = c(0,0)
  
  # Variance Covariance Matrix must be PSD. PSD if both eigen values are greater than zero. 
  if(my_eigen$values[1]<=0 || my_eigen$values[2] <=0){
    # print("Sigma_error_in")
    # print(Sigma_error_in)
    # print("Eigen Value")
    # print(my_eigen$values)
    print("Not PSD Matrix" )
    return(NA)
  }
  
  
  
  # Set parameters of particle filter and initialize storage
  num_particles = 500;
  J = 5
  
  
  # Matrices for Filtered Estimates
  lambda_est_matrix = matrix(0, nrow = Sample, ncol =num_particles)
  lambda_est_matrix[1:J,]=  lambda0_in
  lambda_filtered = matrix(0, nrow = Sample, ncol =1)
  lambda_filtered[1:J]=  lambda0_in
  
  # Initialize log likelihood value
  log_lik = 0;
  
  
  
  for (k in (J+1):(Sample-1)){
    # Generate initial  sample of particles
    # State Equation
    lambda_proposed = mean(lambda_est_matrix[k-1,]) + rnorm(num_particles, mean =0, sd =sig_w_in)
    
    # Measurement Equation
    # weights for different information sets
    w1 = (1-lambda_proposed)
    w2 = lambda_proposed*(1-lambda_est_matrix[k-1,])
    w3 = lambda_proposed*lambda_est_matrix[k-1,]*(1-lambda_est_matrix[k-2,])
    w4 = lambda_proposed*lambda_est_matrix[k-2,]*lambda_est_matrix[k-1,]*(1-lambda_est_matrix[k-3,])
    w5 = lambda_proposed*lambda_est_matrix[k-3,]*lambda_est_matrix[k-2,]*lambda_est_matrix[k-1,]*(1-lambda_est_matrix[k-4,])
    
    # Expectation of optimal price weighted by the percent of firms using that info set
    beliefs_less1 = w1*( inflation_ex[k,1]+ alpha_in*(output_change_ex[k,1] -delta_natural_output[k]) )
    beliefs_less2 = w2*( inflation_ex[k,2]+ alpha_in*(output_change_ex[k,2] -delta_natural_output[k]) )
    beliefs_less3 = w3*( inflation_ex[k,3]+ alpha_in*(output_change_ex[k,3] -delta_natural_output[k]) )
    beliefs_less4 = w4*( inflation_ex[k,4]+ alpha_in*(output_change_ex[k,4] -delta_natural_output[k]) )
    beliefs_less5 = w5*( inflation_ex[k,5]+ alpha_in*(output_change_ex[k,5] -delta_natural_output[k]) )
    
    # vector of exogenous variables in instrument equation
    Zt = c(1,output_gap[k-1],output_gap[k-2],output_gap[k-3],inflation_ex[k+1,1] )
    
    
    # Measurement Equation
    prediction_vec =  ((1-lambda_proposed)/lambda_proposed)*alpha_in* output_gap[k]   +
      beliefs_less1 +  beliefs_less2 + beliefs_less3 +beliefs_less4+beliefs_less5
    
    
    # Simultaneous Equation
    output_gap_prediction = delta_in %*% Zt 
    
    
    # Measure 
    error_vec = inflation[k]- prediction_vec   
    epsilon = output_gap[k]-output_gap_prediction
    
    
    
    # Replicate instrument_error epsilon to match the length of error_vec
    instrument_error_replicated <- rep(epsilon, length(error_vec))
    
    # Combine error_vec and replicated instrument_error_t into a matrix
    X <- cbind(error_vec, instrument_error_replicated)
    
    # Compute the PDF for each row in the matrix
    # likelihood_over_particles <- apply(X, 1, function(x) dmvnorm(x, mean = c(0,0), sigma = Sigma_error_in))
    const <- 1 / (sqrt((2 * pi) ^ length(mu) * det_Sigma_error_in))
    exponent <- exp(-0.5 * rowSums((X - mu) %*% solve_Sigma_error_in * (X - mu)))
    likelihood_over_particles <- const * exponent
    
    
    # Set NA equal to 0 probability
    likelihood_over_particles[is.na(likelihood_over_particles)]=0
    
    if(sum(likelihood_over_particles) ==0)
    {
      print(" sum(likelihood_over_particles) ==0")
      return(NA)
    }
    
    probability_over_particles = likelihood_over_particles/sum(likelihood_over_particles)
    
    
    # Chosen States and Likelihood function
    # Draw from state indices using a weighted probability vector
    chosen_states = sample(1:num_particles,num_particles,replace = TRUE, prob = probability_over_particles)
    
    # save the chosen states
    lambda_chosen = lambda_proposed[chosen_states]
    lambda_est_matrix[k,] = lambda_chosen
    
    # Find mean likelihood over chosen states
    log_pdf_itr = log(mean(likelihood_over_particles[chosen_states]));
    
    # Find total likelihood over time series
    log_lik = log_lik+log_pdf_itr;
    
    # Save mean lambda each period
    lambda_filtered[k] = mean(lambda_chosen) 
    
  }# end for loop particle filter
  
  
  
  # OUTPUT SAMPLES
  particle_list = list()
  particle_list[[1]]= log_lik
  particle_list[[2]]= lambda_filtered
  
  
  return(particle_list)
}# End Particle Filter


  

  
  

  
# Statistics 
#------------------------------------------------------------------------------------ 
Log_MVN_pdf_kernel= function(theta_in, mu_in, V_in){
  # Compute log of kernel of mvn pdf 
  log_pdf_val = -0.5*t(theta_in-mu_in) %*% solve(V_in) %*% (theta_in-mu_in)
  return(log_pdf_val)
}

Log_MVN_pdf = function(theta_in, mu_in, V_in){
  # Compute log of mvn pdf 
  D = length(mu_in)
  det_V = det(V_in)
  log_pdf_val = -0.5 * (D * log(2 * pi) + log(det_V) + t(theta_in - mu_in) %*% solve(V_in) %*% (theta_in - mu_in))
  return(log_pdf_val)
}


Log_GAMMA_pdf_kernel= function(h_in, alpha_in, beta_in){
  # Compute log of kernel of gamma pdf 
  log_pdf_val = (alpha_in-1)*log(h_in) - (h_in/beta_in)
  return(log_pdf_val)
}

Log_GAMMA_pdf = function(h_in, alpha_in, beta_in){
  # Compute log of gamma pdf 
  log_pdf_val = (alpha_in - 1) * log(h_in) - h_in / beta_in - alpha_in * log(beta_in) - lgamma(alpha_in)
  return(log_pdf_val)
}



#Data
# First load the real data, then we load in the estiamtes of the MH sampler
load("MyData.Rda")


N = nrow(MyData)

inflation_ex = matrix(0, N, 5)
inflation_ex[,1]= as.matrix(MyData$expectations_inflation_less1)
inflation_ex[,2]= as.matrix(MyData$expectations_inflation_less2)
inflation_ex[,3]= as.matrix(MyData$expectations_inflation_less3)
inflation_ex[,4]= as.matrix(MyData$expectations_inflation_less4)
inflation_ex[,5]= as.matrix(MyData$expectations_inflation_less5)

output_change_ex = matrix(0, N, 5)
output_change_ex[,1]= as.matrix(MyData$expectations_change_output_less1)
output_change_ex[,2]= as.matrix(MyData$expectations_change_output_less2)
output_change_ex[,3]= as.matrix(MyData$expectations_change_output_less3)
output_change_ex[,4]= as.matrix(MyData$expectations_change_output_less4)
output_change_ex[,5]= as.matrix(MyData$expectations_change_output_less5)


inflation = as.matrix(MyData$inflation)
delta_natural_output = as.matrix(MyData$delta_natural_output)
output_gap =as.matrix(MyData$output_gap)

timing =MyData$t

Z = matrix(0, N, 5)
Z[,1]= as.matrix(MyData$constant)
Z[,2]= as.matrix(MyData$Output_Gap_lag1)
Z[,3]= as.matrix(MyData$Output_Gap_lag2)
Z[,4]= as.matrix(MyData$Output_Gap_lag3)
Z[,5]= as.matrix(MyData$inflation_instrument)




# Load fixed data
#-------------------------------------------------------------------------------
load("Fixed_Data.Rda")

# Access variables from the loaded data
alpha_samples_fixed <- Fixed_Data[[2]]
lambda_samples_fixed <- Fixed_Data[[3]]
hv_samples_fixed <- Fixed_Data[[4]]
covariance_samples_fixed <- Fixed_Data[[5]]
he_samples_fixed  <- Fixed_Data[[6]]
output_gap_IV_fixed <- Fixed_Data[[7]]



# Load TVP data
#-------------------------------------------------------------------------------
load("TVP_Data.Rda")

# Access variables from the loaded data
alpha_samples_tvp <- TVP_Data[[2]]
lambda0_samples_tvp <- TVP_Data[[3]]
hv_samples_tvp <- TVP_Data[[4]]
hw_samples_tvp <- TVP_Data[[5]]
covariance_samples_tvp <- TVP_Data[[6]]
he_samples_tvp <- TVP_Data[[7]]
output_gap_IV_tvp <- TVP_Data[[11]]





Fixed_marginal_likelihood = marginal_likelihood_Fixed( inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z, alpha_samples_fixed, lambda_samples_fixed, hv_samples_fixed, covariance_samples_fixed, he_samples_fixed, output_gap_IV_fixed)
save(Fixed_marginal_likelihood,file="Fixed_marginal_likelihood.Rda") 


TVP_marginal_likelihood =  marginal_likelihood_TVP( inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z, alpha_samples_tvp, lambda0_samples_tvp, hv_samples_tvp, hw_samples_tvp, covariance_samples_tvp, he_samples_tvp, output_gap_IV_tvp)
save(TVP_marginal_likelihood,file="TVP_marginal_likelihood.Rda") 


load("Fixed_marginal_likelihood.Rda")
load("TVP_marginal_likelihood.Rda")

VP_marginal_likelihood[[2]]/ (Fixed_marginal_likelihood[[2]] +TVP_marginal_likelihood[[2]] )



Fixed_marginal_likelihood[[2]]/ (Fixed_marginal_likelihood[[2]] +TVP_marginal_likelihood[[2]] )

 TVP_marginal_likelihood[[2]]/ (Fixed_marginal_likelihood[[2]])