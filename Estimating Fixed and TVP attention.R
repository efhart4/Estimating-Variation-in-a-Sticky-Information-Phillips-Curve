library(MASS)  



# Sampler  
Metropolis_Hastings_Sampler= function( inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z){
  
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
  G0 = 1000;       #NUMBER OF BURN IN DRAWS
  G = 2000;       #NUMBER OF POST CONVERGENCE DRAWS
  total_draws = G0+G  #TOTAL NUMBER OF DRAWS
    
  
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
  
  
  
  
    # Initial Parameters: linear instrument equation
    # # run regression, pull out pieces we need
      Instrumental_regression = lm.fit (Z, output_gap)
      Instrumental_variance = var(Instrumental_regression$residuals)
      coefficients = Instrumental_regression$coefficients
  
      
      
      
    # Initial Parameters: linear instrumental equation   
    # inflation_ex[k,1]+ alpha_in*(output_change_ex[k,1] -delta_natural_output[k]
      yy = inflation
      XX = cbind(inflation_ex, output_change_ex, delta_natural_output)
      Linear_measure_regression = lm.fit (XX[(J+1):(N-1),], yy[(J+1):(N-1)])
      Linear_measure_variance = var(Linear_measure_regression$residuals)
      
      
      
    # Initial Parameters: theta_g
    theta_g = matrix(0,num_param,1)  
    # Normal Parameters
    # # Measure Parameters
      theta_g[1] = 0.1 # alpha
      theta_g[2] = 0.7 # lambda0
      theta_g[3] = Linear_measure_variance^(0.5)*Instrumental_variance^(0.5)*0.5# covariance
    
    # # Instrumental Parameters
      theta_g[4:8] = coefficients
     
    # Gamma Parameters  
    theta_g[9] = Linear_measure_variance^(-0.5)#hv # Measure Error
    theta_g[10] = 100 #hw # Random Walk
    theta_g[11] = Instrumental_variance^(-0.5) # he # Instrumental Error
    #print(theta_g)
  
    # Likelihood of Initial parameters: theta_g
    log_prior_theta_g     = Log_MVN_PDF(theta_g[1:8],mu,V)+ Log_GAMMA_PDF(theta_g[9],alpha_v,beta_v) + Log_GAMMA_PDF(theta_g[10],alpha_w,beta_w)+ Log_GAMMA_PDF(theta_g[11],alpha_e,beta_e)
    #print(log_prior_theta_g)
    
    
    particle_list_g= IV_Particle_Filter(theta_g,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex)
    log_lik_theta_g = particle_list_g[[1]]
    lambda_sample_path_g = particle_list_g[[2]]
    
    if(is.na(log_lik_theta_g))
    {print("theta g has zero probability")
      MedianData = list()
      MedianData[[1]]= FALSE
      return(MedianData)}
    
  
  
  
  
    # storage for samples
    alpha_samples = matrix(0, 1, G)
    lambda0_samples = matrix(0, 1, G)
    hv_samples = matrix(0, 1, G)
    hw_samples = matrix(0, 1, G)
    he_samples = matrix(0, 1, G)
    covariance_samples = matrix(0, 1, G)
    per_sample_weighted_likelihood = matrix(0, 1, G)
    percent_change_samples = matrix(0, num_param, G)
    output_gap_samples = matrix(0,5,G)
    lambda_path_samples = matrix(0, nrow = N, ncol =G)
   
    
    # Counters
    accept=0
    itr=1
  
    
    while(itr<=total_draws){
      
      # using multivariate normal distribution to create a proposal for theta g+1, theta star 
      theta_star = theta_g +mvrnorm(n = 1, mu_prop,R )
      
      # # Print statements to see the variation between draws of the theta vector. 
      if(itr>G0){ percent_change_samples[,itr-G0] = ((theta_star-theta_g)/theta_g)*100}
  
        
      # construct the acceptance probability for theta_star
        # # likelihood across the data
        particle_list= IV_Particle_Filter(theta_star,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex)
        log_lik_theta_star = particle_list[[1]]
        lambda_sample_path_star = particle_list[[2]]
  
        # Particle_lik() will return NA if the proposed parameters produce a vector of states that all have 0 probability
        # If this is the case, the probability of accepting the proposed parameters should be 0. 
        if(is.na(log_lik_theta_star)){
          acceptance_prob =0
          # print("theta_star has 0 chance of acceptance")
        }else{
          # # likelihood of the parameters
          log_prior_theta_star  = Log_MVN_PDF(theta_star[1:8],mu,V) + Log_GAMMA_PDF(theta_star[9],alpha_v,beta_v)+ Log_GAMMA_PDF(theta_star[10],alpha_w,beta_w)+ Log_GAMMA_PDF(theta_star[11],alpha_e,beta_e)
  
          # # collect the pieces and find the final probability of acceptance. 
          log_acceptance_prob = (log_lik_theta_star+log_prior_theta_star)-(log_lik_theta_g+log_prior_theta_g)
          acceptance_prob = exp(log_acceptance_prob)
          acceptance_prob = min(acceptance_prob,1)
          }# end else 
      
      
      
      # accept or reject proposal
      u = runif(1)
      if(u<= acceptance_prob)
        {
          theta_g = theta_star
          lambda_sample_path_g = lambda_sample_path_star
          log_lik_theta_g = log_lik_theta_star
          log_prior_theta_g = log_prior_theta_star
          accept = accept+1
        }

      
      #Print out acceptance rate every 1000th iteration
      if(itr%%1000 ==0){
        print(itr)
        print("Acceptance Rate: ")
        print(accept/itr)
        #print(theta_star)
      }
      
  
      
      
      
      if(itr>G0)
      {
        # Collect samples
        alpha_samples[1,(itr-G0) ] = theta_g[1]
        lambda0_samples[1,(itr-G0) ] = theta_g[2]
        covariance_samples[1,(itr-G0) ] = theta_g[3]
        output_gap_samples[,(itr-G0) ] =theta_g[4:8] 
        
        hv_samples[1,(itr-G0) ] = theta_g[9]
        hw_samples[1,(itr-G0) ] = theta_g[10]
        he_samples[1,(itr-G0) ] = theta_g[11]  
        
        per_sample_weighted_likelihood[1,(itr-G0) ] = log_lik_theta_g+ log_prior_theta_g
        
        lambda_path_samples[,(itr-G0) ]=  lambda_sample_path_g
  
      }
      
      itr = itr+1
    }# endloop                                                                                       
    
    exp_per_sample_weighted_likelihood = exp(per_sample_weighted_likelihood)
    marginal_likelihood = mean(exp_per_sample_weighted_likelihood)
    marginal_likelihood2 = mean(per_sample_weighted_likelihood)
    print(marginal_likelihood)
    print(marginal_likelihood2)
      
      
      
    # OUTPUT SAMPLES
    SummaryData = list()
    SummaryData[[1]]= TRUE
    SummaryData[[2]]= alpha_samples
    SummaryData[[3]]= lambda0_samples
    SummaryData[[4]]= hv_samples
    SummaryData[[5]]= hw_samples
    SummaryData[[6]]= covariance_samples
    SummaryData[[7]]= he_samples
    SummaryData[[8]]= Instrumental_regression$coefficients
    SummaryData[[9]]= lambda_path_samples
    SummaryData[[10]]= percent_change_samples
    SummaryData[[11]]=output_gap_samples
    SummaryData[[12]]= marginal_likelihood
    SummaryData[[13]]= (accept/itr)
    return(SummaryData)
    
}# end main function

  
IV_Particle_Filter = function(theta_in, pi,output_gap,delta_natural_output,inflation_ex,output_change_ex){
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
    
    # Variance Covariance Matrix must be PSD. PSD if both eigen values are greater than zero. 
    if(my_eigen$values[1]<=0 || my_eigen$values[2] <=0){
      # print("Sigma_error_in")
      # print(Sigma_error_in)
      # print("Eigen Value")
      # print(my_eigen$values)
      print("Not PSD Matrix" )
      return(NA)
    }
    
    Sigma_inv = solve(Sigma_error_in)
    Sigma_det = det(Sigma_error_in)
  

  
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
    

    # Evaluate the probability between vector of potential measure errors and instrument error
    likelihood_over_particles = Multivariate_norm_particle_filter(error_vec,epsilon, Sigma_inv,Sigma_det)
    
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
Log_MVN_PDF= function(theta_in, mu_in, V_in){
  # Compute log of kernal of mvn pdf 
  log_pdf_val = -0.5*t(theta_in-mu_in) %*% solve(V_in) %*% (theta_in-mu_in)
  return(log_pdf_val)
}


Log_GAMMA_PDF= function(h_in, alpha_in, beta_in){
  # Compute log of kernal of gamma pdf 
  log_pdf_val = (alpha_in-1)*log(h_in) - (h_in/beta_in)
  return(log_pdf_val)
}



Multivariate_norm_particle_filter = function(error_vec,instrument_error_t, Sigma_inv,Sigma_det){
  particle_count = NROW(error_vec)
  likelihood_over_particles = matrix(0, particle_count,1)
  
  for (i in 1:particle_count){
    likelihood_over_particles[i] = (2*pi)^(-1) * Sigma_det^(-0.5) * exp(-0.5*(c(error_vec[i],instrument_error_t) %*% Sigma_inv %*% c(error_vec[i],instrument_error_t)))
  }
  return(likelihood_over_particles)
}
  
  
  
  
# Sampler  
Metropolis_Hastings_Sampler_Fixed= function( inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z){
  
  set.seed(24689)
  N = nrow(inflation)
  J = 5
  
  
  #Prior Parameters
  # Hyper parameters for Normal distribution describing prior for conditional mean parameters: 
  # lambda0, alpha, covariance and a delta vector for Instrumental Equation
  k = 2 + 5  
  mu = matrix(0, k, 1)
  mu[1] = 0.1 # value for alpha in theory. 
  mu[2] = 0.7 # estimate from other research
  mu[3:k] = 0
  
  
  # Variance of prior 
  # I looked at a normal graph these seem reasonable. 
  V= diag(k) 
  V[1,1]= 0.01^2 # ALPHA
  V[2,2]= 1  # LAMBDA 
  
  
  # set variance of instrument equation high to create a diffuse prior
  for (i in 3:k){V[i,i]= 1}
  
  
  
  # Hyper parameters from Gamma distribution describing prior.E(h)= alpha*beta, variance: alpha*beta^2
  # compound disturbance
  alpha_v=100
  beta_v = 4
  
  

  
  
  # NUMBER OF SIMULATIONS FROM THE MH SAMPLER
  G0 = 1000;       #NUMBER OF BURN IN DRAWS
  G = 2000;       #NUMBER OF POST CONVERGENCE DRAWS
  total_draws = G0+G  #TOTAL NUMBER OF DRAWS
  
  
  # Sampler
  num_param = k+1
  mu_prop = matrix(0, num_param, 1)
  
  
  # # Variance covariance matrix of shock to random walk proposal
  R =  diag(num_param)
  # Measure Parameters
  R[1,1] = 5e-6 # alpha
  R[2,2] = 1e-4# lambda0
  
  # Instrumental Parameters
  R[3,3] = 1e-7 # og_intercept
  R[4,4] = 1e-5 # og_lag1
  R[5,5] = 1e-5 # og_lag2
  R[6,6] = 1e-6 # og_lag3
  R[7,7] = 1e-5# og_inflation_expectation
  
  # Gamma Parameters 
  R[8,8] = 1.5 #hv # compound Error
  
  R_scale = 100 #0.09
  R = R_scale*R
  #print(paste0("R_scale: ", R_scale))
  
  
  
  
  # Initial Parameters: linear instrument equation
  # # run regression, pull out pieces we need
  Instrumental_regression = lm.fit (Z, output_gap)
  Instrumental_variance = var(Instrumental_regression$residuals)
  coefficients = Instrumental_regression$coefficients
  
  
  
  
  # Initial Parameters: linear instrumental equation   
  # inflation_ex[k,1]+ alpha_in*(output_change_ex[k,1] -delta_natural_output[k]
  yy = inflation
  XX = cbind(inflation_ex, output_change_ex, delta_natural_output)
  Linear_measure_regression = lm.fit (XX[(J+1):(N-1),], yy[(J+1):(N-1)])
  Linear_measure_variance = var(Linear_measure_regression$residuals)
  
  
  
  # Initial Parameters: theta_g
  theta_g = matrix(0,num_param,1)  
  # Normal Parameters
  # # Measure Parameters
  theta_g[1] = 0.1 # alpha
  theta_g[2] = 0.7 # lambda
  
  # # Instrumental Parameters
  theta_g[3:k] = coefficients
  
  # Gamma Parameters  
  theta_g[8] = Linear_measure_variance^(-0.5)#hv # Measure compound

  
  
  # Likelihood of Initial parameters: theta_g
  log_prior_theta_g     = Log_MVN_PDF(theta_g[1:7],mu,V)+ Log_GAMMA_PDF(theta_g[8],alpha_v,beta_v) 
  log_lik_theta_g = likelihood_fixed_lambda(theta_g,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z)
  

  # storage for samples
  alpha_samples = matrix(0, 1, G)
  lambda0_samples = matrix(0, 1, G)
  hv_samples = matrix(0, 1, G)
  per_sample_weighted_likelihood = matrix(0, 1, G)
  percent_change_samples = matrix(0, num_param, G)
  output_gap_samples = matrix(0,5,G)

  
  # Counters
  accept=0
  itr=1
  
  
  while(itr<=total_draws){
    
    # using multivariate normal distribution to create a proposal for theta g+1, theta star 
    theta_star = theta_g +mvrnorm(n = 1, mu_prop,R )
    
    # # Print statements to see the variation between draws of the theta vector. 
    if(itr>G0){ percent_change_samples[,itr-G0] = ((theta_star-theta_g)/theta_g)*100}
    
    
    # construct the acceptance probability for theta_star
    log_lik_theta_star = likelihood_fixed_lambda(theta_star,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z)
    if(is.na(log_lik_theta_star)){
      acceptance_prob =0
      # print("theta_star has 0 chance of acceptance")
    }else{
      # # likelihood of the parameters
      log_prior_theta_star  = Log_MVN_PDF(theta_star[1:7],mu,V) + Log_GAMMA_PDF(theta_star[8],alpha_v,beta_v)
      
      # # collect the pieces and find the final probability of acceptance. 
      log_acceptance_prob = (log_lik_theta_star+log_prior_theta_star)-(log_lik_theta_g+log_prior_theta_g)
      acceptance_prob = exp(log_acceptance_prob)
      acceptance_prob = min(acceptance_prob,1)
    }# end else 
    
    
    
    # accept or reject proposal
    u = runif(1)
    if(u<= acceptance_prob)
    {
      theta_g = theta_star
      log_lik_theta_g = log_lik_theta_star
      log_prior_theta_g = log_prior_theta_star
      accept = accept+1
    }
    # else, keep theta_g the same
    
    
    #Print out acceptance rate every 1000th iteration
    if(itr%%1000 ==0){
      print(itr)
      print("Acceptance Rate: ")
      print(accept/itr)
      #print(theta_star)
    }
    
    
    
    
    
    if(itr>G0)
    {
      # Collect theta_g as a sample
      alpha_samples[1,(itr-G0) ] = theta_g[1]
      lambda0_samples[1,(itr-G0) ] = theta_g[2]
      output_gap_samples[,(itr-G0) ] =theta_g[3:k] 
      hv_samples[1,(itr-G0) ] = theta_g[8]    
      per_sample_weighted_likelihood[1,(itr-G0) ] = log_lik_theta_g+ log_prior_theta_g # log_lik_theta_g+log_prior_theta_g
      #print(log_lik_theta_g+ log_prior_theta_g)
      #print(exp(log_lik_theta_g+ log_prior_theta_g))
    }
    
    itr = itr+1
  }# endloop                                                                                       
  
  exp_per_sample_weighted_likelihood = exp(per_sample_weighted_likelihood)
  marginal_likelihood = mean(exp_per_sample_weighted_likelihood)
  marginal_likelihood2 = mean(per_sample_weighted_likelihood)
  print(marginal_likelihood)
  print(marginal_likelihood2)
  
  
  # OUTPUT SAMPLES
  SummaryData = list()
  SummaryData[[1]]= TRUE
  SummaryData[[2]]= alpha_samples
  SummaryData[[3]]= lambda0_samples
  SummaryData[[4]]= hv_samples 
  SummaryData[[5]]= Instrumental_regression$coefficients
  SummaryData[[6]]= percent_change_samples
  SummaryData[[7]]=output_gap_samples
  SummaryData[[8]]= marginal_likelihood
  SummaryData[[9]]= (accept/itr)
  return(SummaryData)
  
}# end main function
  
  
  
  
likelihood_fixed_lambda = function(theta_in,inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z){

  Sample = NROW(inflation)
  # Theta in 
  alpha_in   = theta_in[1]  
  lambda_in = theta_in[2] 
  delta_in = theta_in[3:7]
  
  # Measure Error   
  hv_in = theta_in[8]  
  #sig_v_in = hv_in^(-1)

  J = 5
  
  # Initialize log likelihood value
  log_lik = 0;


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
  prediction_vec =  ((1-lambda_in)/lambda_in)*alpha_in* output_gap_prediction + beliefs_less1 +  beliefs_less2 + beliefs_less3 +beliefs_less4+beliefs_less5

  # Measure 
  error_vec = inflation- prediction_vec   
  epsilon = output_gap-output_gap_prediction


  # Likelihood
  # Error from using theta in
  #resid = Y - (gamma1*(X1^gamma3) + gamma2*(X2^gamma3))^(1/gamma3)
  
  # likehood of the outcome conditioned on theta and independent variables
  #log_lik_val = -(N/2)*log(2*pi) + (N/2)*log(h_in) -(h_in/2)*(t(resid) %*% resid)
  
  #log_lik_val = -log_lik_val
  
  log_lik_val = -(Sample/2)*log(2*pi) + (Sample/2)*log(hv_in) -(hv_in/2)*(t(error_vec) %*% error_vec)
  
  return(log_lik_val)    
}








#Data
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


# # Extracting Estimates

Fixed_Data = Metropolis_Hastings_Sampler_Fixed( inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z)
Fixed_Data[[10]] = timing
#save(Fixed_Data,file="Fixed_Data.Rda") 





TVP_Data =  Metropolis_Hastings_Sampler( inflation,output_gap,delta_natural_output,inflation_ex,output_change_ex, Z)
#save(TVP_Data,file="TVP_Data.Rda") 

