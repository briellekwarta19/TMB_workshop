FIRST_SEED <- 300
# Prepare the data & parameters, then estimate for different numbers of hidden states
set.seed(FIRST_SEED)
# Parameters and covariates --------------------------
m <- M_LIST_SIMU3
if (m == 2) {
  true_gamma <- matrix(c(0.95, 0.05,
                         0.15, 0.85),
                       byrow = TRUE,
                       nrow = m,
                       ncol = m)
} else if (m == 3) {
  true_gamma <- matrix(c(0.950, 0.025, 0.025,
                         0.050, 0.900, 0.050,
                         0.075, 0.075, 0.850),
                       byrow = TRUE,
                       nrow = m,
                       ncol = m)
} else if (m == 4) {
  true_gamma <- matrix(c(0.850, 0.050, 0.050, 0.05,
                         0.050, 0.850, 0.050, 0.05,
                         0.050, 0.100, 0.800, 0.05,
                         0.034, 0.033, 0.033, 0.90),
                       byrow = TRUE,
                       nrow = m,
                       ncol = m)
}
true_lambda <- seq(1,
                   13,
                   length.out = m)
true_delta <- stat.dist(true_gamma)

simu3_data <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                 mod = list(m = m,
                                                            lambda = true_lambda,
                                                            gamma = true_gamma,
                                                            delta = true_delta),
                                                 testing_params = list(m = m,
                                                                       lambda = true_lambda,
                                                                       gamma = true_gamma,
                                                                       delta = true_delta))$data


# Parameters & covariates for TMB ------------------
working_true_params <- pois.HMM.pn2pw(m = m,
                                      lambda = true_lambda,
                                      gamma = true_gamma,
                                      delta = true_delta)
TMB_data <- list(x = simu3_data,
                 m = m)

# Estimation ------------------------------------
dm <- DM.estimate(x = simu3_data,
                  m = m,
                  lambda0 = true_lambda,
                  gamma0 = true_gamma,
                  delta0 = true_delta)
tmb <- TMB.estimate(TMB_data = TMB_data,
                    parameters = working_true_params)
tmb_g <- TMB.estimate(TMB_data = TMB_data,
                      parameters = working_true_params,
                      gradient = TRUE)
tmb_h <- TMB.estimate(TMB_data = TMB_data,
                      parameters = working_true_params,
                      hessian = TRUE)
tmb_gh <- TMB.estimate(TMB_data = TMB_data,
                       parameters = working_true_params,
                       gradient = TRUE,
                       hessian = TRUE,
                       std_error = TRUE)

# If one doesn't converge successfully, stop
if (dm$convergence != 0) {
  stop(paste("dm didn't converge properly, simu3 dataset, m =",
             m))
}
if (tmb$convergence != 0) {
  stop(paste("tmb didn't converge properly, simu3 dataset, m =",
             m))
}
if (tmb_g$convergence != 0) {
  stop(paste("tmb_g didn't converge properly, simu3 dataset, m =",
             m))
}
if (tmb_h$convergence != 0) {
  stop(paste("tmb_h didn't converge properly, simu3 dataset, m =",
             m))
}
if (tmb_gh$convergence != 0) {
  stop(paste("tmb_gh didn't converge properly, simu3 dataset, m =",
             m))
}

# Creating variables for the CIs -----------------
params_names_latex <- paste0(rep("$\\lambda_{",
                                 m),
                             1:m,
                             "}$")
for (col in 1:m) {
  # Get row and column indices for gamma instead of the default
  # columnwise index: the default indices are 1:m for the 1st column,
  # then (m + 1):(2 * m) for the 2nd, etc...
  params_names_latex <- c(params_names_latex,
                          paste0(sapply(X = 1:m,
                                        FUN = function(row) {paste0("$\\gamma_{",
                                                                    row,
                                                                    col,
                                                                    "}$")})))
}
params_names_latex <- c(params_names_latex,
                        paste0(rep("$\\delta_{",
                                   m),
                               1:m,
                               "}$"))
len_par <- length(params_names_latex)
indices <- (length(conf_int_simu3$m) + 1):(length(conf_int_simu3$m) + len_par)
conf_int_simu3[indices, "m"] <- m
conf_int_simu3[indices, "Parameter"] <- params_names_latex
# Reminder, PARAMS_NAMES contains c("lambda", "gamma", "delta")
conf_int_simu3[indices, "Estimate"] <- unlist(tmb_gh[PARAMS_NAMES])
conf_int_simu3[indices, "True.value"] <- as.numeric(c(true_lambda,
                                                      true_gamma,
                                                      true_delta))

param_tmb_gh <- pois.HMM.pn2pw(m = m,
                               lambda = tmb_gh$lambda,
                               gamma = tmb_gh$gamma)

if (m == 1) {
  w_params_names <- c("tlambda1")
} else {
  w_params_names <- c(paste0(rep("tlambda",
                                 m),
                             1:m),
                      paste0(rep("tgamma",
                                 m ^ 2 - m),
                             1:(m ^ 2 - m)))
}
len_w_par <- length(w_params_names)

lambda_indices <- 1:m
gamma_indices <- max(lambda_indices) + 1:(m ^ 2)
delta_indices <- max(gamma_indices) + 1:m
tgamma_indices <- (m + 1):(m ^ 2)

# Benchmarks ------------
set.seed(FIRST_SEED + 1)
if (BENCHMARK_SAMPLES != 0) {
  for (idx_counter in 1:BENCHMARK_SAMPLES) {
    # Generate data that can be estimated by TMB_GH
    # and is tested on the slightly off parameters from the beginning of this file
    # The goal is to have a dataset that poses no estimation problem
    benchmark_model <- pois.HMM.generate.estimable.sample(
      ns = DATA_SIZE_SIMU3,
      mod = list(m = m,
                 lambda = tmb_gh$lambda,
                 gamma = tmb_gh$gamma,
                 delta = tmb_gh$delta),
      testing_params = list(m = m,
                            lambda = true_lambda,
                            gamma = true_gamma,
                            delta = true_delta)
    )
    benchmark_data <- benchmark_model$data
    # Benchmark all different combinations of gradient and hessians with DM ----------------
    # Parameters & covariates for DM and TMB
    TMB_benchmark_data <- list(x = benchmark_data,
                               m = m)
    obj_benchmark <- MakeADFun(TMB_benchmark_data,
                               working_true_params,
                               DLL = "poi_hmm",
                               silent = TRUE)
    parvect_benchmark_TMB <- pois.HMM.pn2pw(m = m,
                                            lambda = true_lambda,
                                            gamma = true_gamma,
                                            delta = true_delta)
    # nlminb needs a vector, not a list
    parvect_benchmark_DM <- unlist(parvect_benchmark_TMB)
    # Estimation benchmark
    temp <- microbenchmark("DM" = nlminb(parvect_benchmark_DM,
                                         pois.HMM.mllk,
                                         x_alias = benchmark_data,
                                         m_alias = m)$convergence==0,
                           "TMB" = nlminb(obj_benchmark$par,
                                          obj_benchmark$fn)$convergence==0,
                           "TMB_G" = nlminb(obj_benchmark$par,
                                            obj_benchmark$fn,
                                            gradient = obj_benchmark$gr)$convergence==0,
                           "TMB_H" = nlminb(obj_benchmark$par,
                                            obj_benchmark$fn,
                                            hessian = obj_benchmark$he)$convergence==0,
                           "TMB_GH" = nlminb(obj_benchmark$par,
                                             obj_benchmark$fn,
                                             gradient = obj_benchmark$gr,
                                             hessian = obj_benchmark$he)$convergence==0,
                           times = 1,
                           check = "equal",
                           setup = obj_benchmark <<- MakeADFun(TMB_benchmark_data,
                                                               working_true_params,
                                                               DLL = "poi_hmm",
                                                               silent = TRUE))
    times <- temp$time / 10^9
    timeDM <- times[temp$expr == "DM"]
    timeTMB <- times[temp$expr == "TMB"]
    timeTMB_G <- times[temp$expr == "TMB_G"]
    timeTMB_H <- times[temp$expr == "TMB_H"]
    timeTMB_GH <- times[temp$expr == "TMB_GH"]
    
    iterDM <- nlminb(parvect_benchmark_DM,
                     pois.HMM.mllk,
                     x_alias = benchmark_data,
                     m_alias = m)$iterations
    iterTMB <- nlminb(obj_benchmark$par,
                      obj_benchmark$fn)$iterations
    iterTMB_G <- nlminb(obj_benchmark$par,
                        obj_benchmark$fn,
                        gradient = obj_benchmark$gr)$iterations
    iterTMB_H <- nlminb(obj_benchmark$par,
                        obj_benchmark$fn,
                        hessian = obj_benchmark$he)$iterations
    iterTMB_GH <- nlminb(obj_benchmark$par,
                         obj_benchmark$fn,
                         gradient = obj_benchmark$gr,
                         hessian = obj_benchmark$he)$iterations
    
    estim_benchmarks_df_simu3 <- rbind(estim_benchmarks_df_simu3,
                                       data.frame(time = c(timeDM,
                                                           timeTMB,
                                                           timeTMB_G,
                                                           timeTMB_H,
                                                           timeTMB_GH),
                                                  m = rep(m,
                                                          length(PROCEDURES)),
                                                  procedure = PROCEDURES,
                                                  iterations = c(iterDM,
                                                                 iterTMB,
                                                                 iterTMB_G,
                                                                 iterTMB_H,
                                                                 iterTMB_GH),
                                                  dataset_number = rep(idx_counter,
                                                                       length(PROCEDURES))))
    
  }
}

# Profiling the likelihood --------------------------
registerDoParallel(cores = CORES)
working_conf_int <- foreach (idx_param = 1:len_w_par,
                             .packages = "TMB",
                             .inorder = TRUE,
                             .combine = rbind) %dopar% {
                               # TMB::compile("code/poi_hmm.cpp")
                               dyn.load(dynlib("code/poi_hmm"))
                               profile <- tmbprofile(obj = tmb_gh$obj,
                                                     name = idx_param,
                                                     trace = FALSE)
                               
                               ci <- tryCatch({
                                 confint(profile)
                               },
                               error = function(e){
                                 return(rep(NA,
                                            2))
                               })
                               
                               dyn.unload(dynlib("code/poi_hmm"))
                               return(ci)
                             }
stopImplicitCluster()
working_conf_int <- as.data.frame(working_conf_int)
rownames(working_conf_int) <- w_params_names

# Transform the working parameters into natural ones
# Lambda (m values)
conf_int_simu3$Profile.L[lambda_indices] <- exp(working_conf_int$lower[lambda_indices])
# Gamma (m^2-m working parameters, m^2 natural ones)
if (!anyNA(working_conf_int$lower[tgamma_indices])) {
  natural_gamma <- as.numeric(gamma.w2n(m,
                                        working_conf_int$lower[tgamma_indices]))
  conf_int_simu3$Profile.L[gamma_indices] <- natural_gamma
}
# Lambda (m values)
conf_int_simu3$Profile.U[lambda_indices] <- exp(working_conf_int$upper[lambda_indices])
# Gamma (m^2-m working parameters, m^2 natural ones)
if (!anyNA(working_conf_int$upper[tgamma_indices])) {
  natural_gamma <- as.numeric(gamma.w2n(m,
                                        working_conf_int$upper[tgamma_indices]))
  conf_int_simu3$Profile.U[gamma_indices] <- natural_gamma
}

# Bootstrap ---------------------------
set.seed(FIRST_SEED + 2)
if (BOOTSTRAP_SAMPLES != 0) {
  registerDoParallel(cores = CORES)
  bootstrap_simu3 <- foreach (idx_sample = 1:BOOTSTRAP_SAMPLES,
                              .packages = "TMB",
                              .combine = rbind) %dopar% {
                                # TMB::compile("code/poi_hmm.cpp")
                                dyn.load(dynlib("code/poi_hmm"))
                                temp <- pois.HMM.generate.estimable.sample(
                                  ns = DATA_SIZE_SIMU3,
                                  mod = list(m = m,
                                             lambda = tmb_gh$lambda,
                                             gamma = tmb_gh$gamma,
                                             delta = tmb_gh$delta),
                                  testing_params = list(m = m,
                                                        lambda = true_lambda,
                                                        gamma = true_gamma,
                                                        delta = true_delta)
                                )
                                # The values from gamma are taken columnwise
                                natural_parameters <- temp$natural_parameters
                                natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
                                dyn.unload(dynlib("code/poi_hmm"))
                                return(natural_parameters)
                              }
  stopImplicitCluster()
  if (BOOTSTRAP_SAMPLES == 1) {
    names(bootstrap_simu3) <- params_names_latex
  } else {
    colnames(bootstrap_simu3) <- params_names_latex
  }
  q <- apply(bootstrap_simu3,
             2,
             quantile.colwise)
  conf_int_simu3$Bootstrap.L <- q[1, ]
  conf_int_simu3$Bootstrap.U <- q[2, ]
}

# TMB confidence intervals --------------
# Manually cap values at their natural bound
# lambda must be strictly above 0
lambda_L <- pmax(0.0001,
                 tmb_gh$lambda - q95_norm * tmb_gh$lambda_std_error)
# gamma must be 0 or more
gamma_L <- pmax(0,
                tmb_gh$gamma - q95_norm * tmb_gh$gamma_std_error)
# delta must be 0 or above
delta_L <- pmax(0,
                tmb_gh$delta - q95_norm * tmb_gh$delta_std_error)
conf_int_simu3$TMB.L <- c(lambda_L,
                          gamma_L,
                          delta_L)
# no upper bound on lambda
# gamma must be 1 or less
gamma_U <- pmin(1,
                tmb_gh$gamma + q95_norm * tmb_gh$gamma_std_error)
# delta must be 1 or less
delta_U <- pmin(1,
                tmb_gh$delta + q95_norm * tmb_gh$delta_std_error)
conf_int_simu3$TMB.U <- c(tmb_gh$lambda + q95_norm * tmb_gh$lambda_std_error,
                          gamma_U,
                          delta_U)
# Coverage probabilities of the 3 CI methods -----------------
set.seed(FIRST_SEED + 3)
parameter_names <- paste0(rep("lambda",
                              m),
                          1:m)
for (col in 1:m) {
  # Get row and column indices for gamma instead of the default
  # columnwise index: the default indices are 1:m for the 1st column,
  # then (m + 1):(2 * m) for the 2nd, etc...
  parameter_names <- c(parameter_names,
                       paste0(sapply(X = 1:m,
                                     FUN = function(row) {paste0("gamma",
                                                                 row, ", ",
                                                                 col)})))
}
parameter_names <- c(parameter_names,
                     paste0(rep("delta",
                                m),
                            1:m))

# Record the times where the profile CI of a parameter successfully contains the parameter's true value
coverage_count_profile_simu3 <-
  coverage_count_bootstrap <-
  coverage_count_tmb <-
  data.frame(parameter = parameter_names,
             count = 0,
             ratio = 0,
             missing = 0)
coverage_count_bootstrap <- coverage_count_tmb <- data.frame(parameter = parameter_names,
                                                             count = 0,
                                                             ratio = 0,
                                                             missing = 0)
idx_coverage <- 0
while (idx_coverage < COVERAGE_SAMPLES) {
  idx_coverage <- idx_coverage + 1
  # Generate a data sample where nlminb converges
  # Loop as long as there is an issue with nlminb
  # Estimate a model
  # Unlike with the other datasets, we know the true parameters of this one
  coverage_model <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                       mod = list(m = m,
                                                                  lambda = true_lambda,
                                                                  gamma = true_gamma,
                                                                  delta = true_delta),
                                                       testing_params = list(m = m,
                                                                             lambda = true_lambda,
                                                                             gamma = true_gamma,
                                                                             delta = true_delta),
                                                       std_error = TRUE)
  
  # Save the occurrences of failures to generate a sample for which parameters can be estimated
  for (reason in c("state_number", "TMB_null", "TMB_converge", "TMB_G_null",
                   "TMB_G_converge", "TMB_H_null", "TMB_H_converge", "TMG_GH_null",
                   "TMG_GH_converge", "NA_value")) {
    coverage_skips_simu3[coverage_skips_simu3$m == m, reason] <-
      coverage_skips_simu3[coverage_skips_simu3$m == m, reason] +
      coverage_model$failure[reason]
  }
  
  # # Confidence interval profiling -------------------------------------
  # Less than half of the coverage samples have converging profile CIs, so we do not include
  # them in the paper. Therefore, there is no need to run this part (lines 453-535)
  # If you want to verify this, feel free to uncomment lines 453-535 and run the code below
  # then look at the variable coverage_count_profile_simu3
  # The column "missing" indicates how many times a parameter has no profile CI.
  # Here is the content of that variable that was generated
  # 
  # > coverage_count_profile_simu3
  # parameter count     ratio missing
  # 1    lambda1   957 0.9570000       0
  # 2    lambda2   939 0.9390000       0
  # 3    lambda3   937 0.9370000       0
  # 4    lambda4   953 0.9530000       0
  # 5  gamma1, 1   854 1.0000000     146
  # 6  gamma2, 1   653 0.8932969     269
  # 7  gamma3, 1   374 0.9234568     595
  # 8  gamma4, 1   463 0.9042969     488
  # 9  gamma1, 2   804 0.9414520     146
  # 10 gamma2, 2   731 1.0000000     269
  # 11 gamma3, 2   364 0.8987654     595
  # 12 gamma4, 2   502 0.9804688     488
  # 13 gamma1, 3   823 0.9637002     146
  # 14 gamma2, 3   705 0.9644323     269
  # 15 gamma3, 3   405 1.0000000     595
  # 16 gamma4, 3   510 0.9960938     488
  # 17 gamma1, 4   781 0.9145199     146
  # 18 gamma2, 4   694 0.9493844     269
  # 19 gamma3, 4   385 0.9506173     595
  # 20 gamma4, 4   512 1.0000000     488
  # 21    delta1     0        NA       0
  # 22    delta2     0        NA       0
  # 23    delta3     0        NA       0
  # 24    delta4     0        NA       0
  # 
  # Profile CIs for the third row of the TPM did not converge 595 times out of 1000.
  # For information, the column "count" counts the times where the profile CI successfully
  # contained the true parameter value.
  # The column "ratio" is "count" / (1000 - "missing").
  # 
  # registerDoParallel(cores = CORES)
  # working_conf_int <- foreach (idx_param = 1:len_w_par,
  #                              .packages = "TMB",
  #                              .inorder = TRUE,
  #                              .combine = rbind) %dopar% {
  #                                # TMB::compile("code/poi_hmm.cpp")
  #                                dyn.load(dynlib("code/poi_hmm"))
  #                                profile <- tmbprofile(obj = coverage_model$mod$obj,
  #                                                      name = idx_param,
  #                                                      trace = FALSE)
  #                                
  #                                ci <- tryCatch({
  #                                  confint(profile)
  #                                },
  #                                error = function(e){
  #                                  return(rep(NA,
  #                                             2))
  #                                })
  #                                dyn.unload(dynlib("code/poi_hmm"))
  #                                return(ci)
  #                              }
  # stopImplicitCluster()
  # working_conf_int <- as.data.frame(working_conf_int)
  # rownames(working_conf_int) <- w_params_names
  # 
  # lambda_profile_lower <- exp(working_conf_int$lower[lambda_indices])
  # lambda_profile_upper <- exp(working_conf_int$upper[lambda_indices])
  # gamma_profile_lower <- as.numeric(gamma.w2n(m,
  #                                             working_conf_int$lower[tgamma_indices]))
  # gamma_profile_upper <- as.numeric(gamma.w2n(m,
  #                                             working_conf_int$upper[tgamma_indices]))
  # 
  ### MAIN DIFFERENCE FROM THE OTHER DATASETS !!!
  ###
  ### Unlike in the other datasets, we have empirically found that profile CIs diverge when
  ### lambda=seq(1, 13, length.out = m)
  ### Therefore, requiring profile CIs to converge will cause an infinite loop since we chose
  ### lambda=seq(1, 13, length.out = m) to showcase that profile CIs sometimes diverge.
  # 
  # # If profiling doesn't yield results for all parameters, record it (because we don't guarantee
  # # profile CIs in this simulation)
  # estimates_coverage <- c(lambda_profile_lower,
  #                         lambda_profile_upper,
  #                         gamma_profile_lower,
  #                         gamma_profile_upper)
  # 
  # test_null <- sapply(X = estimates_coverage, FUN = is.null)
  # test_finite <- sapply(X = estimates_coverage, FUN = is.finite)
  # # If some CI bounds are NULL or missing (NA) or infinite (Inf), record it
  # if (any(test_null == TRUE) | any(test_finite == FALSE)) {
  #   test_null_lower <- sapply(X = c(lambda_profile_lower,
  #                                   gamma_profile_lower), FUN = is.null)
  #   test_null_upper <- sapply(X = c(lambda_profile_upper,
  #                                   gamma_profile_upper), FUN = is.null)
  #   test_finite_lower <- sapply(X = c(lambda_profile_lower,
  #                                     gamma_profile_lower), FUN = is.finite)
  #   test_finite_upper <- sapply(X = c(lambda_profile_upper,
  #                                     gamma_profile_upper), FUN = is.finite)
  #   
  #   indices_lower <- which(test_null_lower | !test_finite_lower)
  #   indices_upper <- which(test_null_upper | !test_finite_upper)
  #   indices <- union(indices_lower, indices_upper)
  #   coverage_count_profile_simu3[indices, "missing"] <- coverage_count_profile_simu3[indices, "missing"] + 1
  # }
  # 
  # # If the true value for lambda is in the CI, then increase the count
  # real_lambda_profile_lower <- pmin(lambda_profile_lower,
  #                                   lambda_profile_upper)
  # real_lambda_profile_upper <- pmax(lambda_profile_lower,
  #                                   lambda_profile_upper)
  # indices <- which(true_lambda >= real_lambda_profile_lower &
  #                    true_lambda <= real_lambda_profile_upper)
  # coverage_count_profile_simu3[indices, "count"] <- coverage_count_profile_simu3[indices, "count"] + 1
  # 
  # # Same for gamma
  # real_gamma_profile_lower <- pmin(gamma_profile_lower,
  #                                  gamma_profile_upper)
  # real_gamma_profile_upper <- pmax(gamma_profile_lower,
  #                                  gamma_profile_upper)
  # indices <- which(as.vector(true_gamma) >= real_gamma_profile_lower &
  #                    as.vector(true_gamma) <= real_gamma_profile_upper)
  # indices <- indices + m
  # coverage_count_profile_simu3[indices, "count"] <- coverage_count_profile_simu3[indices, "count"] + 1
  
  # Confidence interval bootstrap -----------------------------------
  if (BOOTSTRAP_SAMPLES != 0) {
    registerDoParallel(cores = CORES)
    bootstrap_simu3 <- foreach (idx_sample = 1:BOOTSTRAP_SAMPLES,
                                .packages = "TMB",
                                .combine = rbind) %dopar% {
                                  # TMB::compile("code/poi_hmm.cpp")
                                  dyn.load(dynlib("code/poi_hmm"))
                                  temp <- pois.HMM.generate.estimable.sample(
                                    ns = DATA_SIZE_SIMU3,
                                    mod = list(m = m,
                                               lambda = coverage_model$natural_parameters$lambda,
                                               gamma = coverage_model$natural_parameters$gamma,
                                               delta = coverage_model$natural_parameters$delta),
                                    testing_params = list(m = m,
                                                          lambda = true_lambda,
                                                          gamma = true_gamma,
                                                          delta = true_delta))
                                  # The values from gamma are taken columnwise
                                  natural_parameters <- temp$natural_parameters
                                  natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
                                  dyn.unload(dynlib("code/poi_hmm"))
                                  return(natural_parameters)
                                }
    stopImplicitCluster()
    q <- apply(bootstrap_simu3,
               2,
               quantile.colwise)
    indices <- which(as.vector(true_lambda) >= q[1, lambda_indices] &
                       as.vector(true_lambda) <= q[2, lambda_indices])
    coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
    indices <- which(as.vector(true_gamma) >= q[1, gamma_indices] &
                       as.vector(true_gamma) <= q[2, gamma_indices]) + m
    coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
    indices <- which(as.vector(true_delta) >= q[1, delta_indices] &
                       as.vector(true_delta) <= q[2, delta_indices]) + m + m ^ 2
    coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
  }
  
  # Confidence interval TMB -----------------------------------------
  nat_par <- coverage_model$natural_parameters
  lambda_tmb_lower <- nat_par$lambda - q95_norm * nat_par$lambda_std_error
  lambda_tmb_upper <- nat_par$lambda + q95_norm * nat_par$lambda_std_error
  gamma_tmb_lower <- nat_par$gamma - q95_norm * nat_par$gamma_std_error
  gamma_tmb_upper <- nat_par$gamma + q95_norm * nat_par$gamma_std_error
  delta_tmb_lower <- nat_par$delta - q95_norm * nat_par$delta_std_error
  delta_tmb_upper <- nat_par$delta + q95_norm * nat_par$delta_std_error
  
  indices <- which(as.vector(true_lambda) >= lambda_tmb_lower &
                     as.vector(true_lambda) <= lambda_tmb_upper)
  coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
  indices <- which(as.vector(true_gamma) >= gamma_tmb_lower &
                     as.vector(true_gamma) <= gamma_tmb_upper)
  indices <- indices + m
  coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
  indices <- which(as.vector(true_delta) >= delta_tmb_lower &
                     as.vector(true_delta) <= delta_tmb_upper)
  indices <- indices + m + m ^ 2
  coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
}
temp <- coverage_count_profile_simu3[lambda_indices, ]
coverage_count_profile_simu3[lambda_indices, "ratio"] <- temp[, "count"] /
  (COVERAGE_SAMPLES - temp[, "missing"])
temp <- coverage_count_profile_simu3[gamma_indices, ]
coverage_count_profile_simu3[gamma_indices, "ratio"] <- temp[, "count"] /
  (COVERAGE_SAMPLES - temp[, "missing"])
# delta is not a parameter for us, so it has no profile CI
coverage_count_profile_simu3[delta_indices, "ratio"] <- NA

coverage_count_tmb$ratio <- coverage_count_tmb$count / COVERAGE_SAMPLES

coverage_count_bootstrap$ratio <- coverage_count_bootstrap$count / COVERAGE_SAMPLES

# Most of the time, profile CIs for the stationary distribution is NA.
conf_int_simu3[1:(m ^ 2 + 2 * m), "Coverage.Profile"] <- coverage_count_profile_simu3$ratio * 100
conf_int_simu3[1:(m ^ 2 + 2 * m), "Coverage.Bootstrap"] <- coverage_count_bootstrap$ratio * 100
conf_int_simu3[1:(m ^ 2 + 2 * m), "Coverage.TMB"] <- coverage_count_tmb$ratio * 100

# Fixes -------------------------
# Fix label switching in conf_int_simu3
ordered_params <- pois.HMM.label.order(m = m,
                                       lambda = true_lambda,
                                       gamma = true_gamma,
                                       delta = true_delta)

new_lambda_indices <- ordered_params$ordered_lambda_indices
new_gamma_indices <- ordered_params$ordered_gamma_vector_indices
new_delta_indices <- ordered_params$ordered_delta_indices

conf_int_simu3[lambda_indices, - 2] <- conf_int_simu3[lambda_indices, - 2][new_lambda_indices, ]
conf_int_simu3[gamma_indices, - 2] <- conf_int_simu3[gamma_indices, - 2][new_gamma_indices, ]
conf_int_simu3[delta_indices, - 2] <- conf_int_simu3[delta_indices, - 2][new_delta_indices, ]

# Reorder the TPM row-wise instead of column-wise
# Lexicographical parameter sort for gamma (sort on the parameter name)
new_gamma_indices_truncated_table <- order(conf_int_simu3[gamma_indices, "Parameter"])
# Replace rows by sorted rows
conf_int_simu3[gamma_indices, ] <- conf_int_simu3[gamma_indices, ][new_gamma_indices_truncated_table, ]

# The profile CIs may not be sorted, so we sort them manually
for (i in 1:length(conf_int_simu3[, 1])) {
  row <- conf_int_simu3[i, c("Profile.L", "Profile.U")]
  conf_int_simu3[i, c("Profile.L", "Profile.U")] <- cbind(min(row), max(row))
}
conf_int_simu3$m <- as.integer(conf_int_simu3$m)

estim_benchmarks_df_simu3$m <- factor(estim_benchmarks_df_simu3$m,
                                      levels = M_LIST_SIMU3)
