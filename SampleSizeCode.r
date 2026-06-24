# Written 06.20.26 
# By A.Williams

# Code for Aim2: Power calculations for the fixed effects model for 
# stepped wedge designs with a small number of clusters

# Code to calculate the power and required sample size for a
# stepped wedge design for the fixed effects model for continuous outcomes

# Model 1 (fixed effects model, with time effects): 
# Outcome for cluster i in period j for individual k:
# Y_{ijk} = sum_{j = 1}^{T-1} \beta_j I(step = j) + 
#                  sum_{i = 1}^{I} \alpha_j I(cluster = j) + 
#                  \theta_1*X_{ij} + e_{ijk}
# Assume individual level error e_{ijk} ~ N(0, sigma_e^2)

# Model 2 (fixed effects model, without time effects): 
# Outcome for cluster i in period j for individual k:
# Y_{ijk} = sum_{i = 1}^{I} \alpha_j I(cluster = j) + 
#                  \theta_1*X_{ij} + e_{ijk}
# Assume individual level error e_{ijk} ~ N(0, sigma_e^2)


########################################################################
# Power of fixed effects model for SWD, time effects included in model (Model 1)
########################################################################
power_model1_general <- function(num_I, num_T, num_N, ICC, Var_Y, design, theta, alpha) {
  # Assumptions:
  
  # Analysis model: Model 1
  
  # num_I     : number of clusters
  # num_T     : number of periods
  # num_N     : number of observations sampled per cluster-period
  # ICC       : Intraclass correlation coefficient
  # Var_Y     : total variance, sigma_e^2 = Var_Y*(1-ICC)
  # design    : Matrix specifying the stepped wedge treatment allocation scheme. 
  #                   Rows correspond to clusters and columns correspond to time periods.
  #                   Entries are coded as 0 for control and 1 for intervention.
  # theta     : intervention effect
  # alpha     : type 1 error rate
  
  stopifnot(
    num_I > 0,
    num_T > 0,
    num_N > 0,
    ICC >= 0,
    ICC < 1,
    Var_Y > 0,
    alpha > 0,
    alpha < 1
  )
  
  z_alpha <- qnorm(1-alpha/2)
  
  U  <- sum(design)              # SUM_i (SUM_j X_{ij})
  W  <- sum(colSums(design)^2)   # SUM_j (SUM_i X_ij^2)
  V  <- sum(rowSums(design)^2)   # SUM_i (SUM_j X_ij^2)
  
  # Calculate design constants
  var_theta <- Var_Y*(1-ICC)*num_I*num_T/(num_N*(num_I*num_T*U - num_T*W + U^2-num_I*V))
  power     <- pnorm(abs(theta)/sqrt(var_theta) -z_alpha) + pnorm(-abs(theta)/sqrt(var_theta) -z_alpha)
  
  return(power)
}


########################################################################
# Required sample size of fixed effects model for SWD, time effects
# included in model (Model 1)
########################################################################
sample_size_model1_general <- function(num_I, num_T, ICC, Var_Y, design, theta, alpha, power) {
  # Returns the per cluster period sample size required to achieve the desired power
  
  # Assumptions:
  
  # Analysis model: Model 1
  
  # num_I     : number of clusters
  # num_T     : number of periods
  # ICC       : Intraclass correlation coefficient
  # Var_Y     : total variance, sigma_e^2 = Var_Y*(1-ICC)
  # design    : Matrix specifying the stepped wedge treatment allocation scheme. 
  #                   Rows correspond to clusters and columns correspond to time periods.
  #                   Entries are coded as 0 for control and 1 for intervention.
  # theta     : intervention effect
  # alpha     : type 1 error rate
  # power     : desired power, as fraction i.e. 0.80
  
  stopifnot(
    num_I > 0,
    num_T > 0,
    ICC >= 0,
    ICC < 1,
    Var_Y > 0,
    alpha > 0,
    alpha < 1,
    power > 0,
    power < 1
  )
  
  z_alpha <- qnorm(1-alpha/2)
  z_beta  <- qnorm(power)
  
  U  <- sum(design)              # SUM_i (SUM_j X_{ij})
  W  <- sum(colSums(design)^2)   # SUM_j (SUM_i X_ij^2)
  V  <- sum(rowSums(design)^2)   # SUM_i (SUM_j X_ij^2)
  
  # Calculate design constants
  num_N <- Var_Y*(1-ICC)*num_I*num_T/(num_I*num_T*U - num_T*W + U^2-num_I*V)*((z_alpha + z_beta)/theta)^2
  
  return(ceiling(num_N))
}

# Example 1:
# Real world examples based on PROMPT Study
# Model 1, no transition period included, assuming 5 clusters and 6 time periods
design_example <- matrix(c(0,1,1,1,1,1,
                   0,0,1,1,1,1,
                   0,0,0,1,1,1,
                   0,0,0,0,1,1,
                   0,0,0,0,0,1), 
                 nrow = 5, ncol = 6, byrow =TRUE)

power_model1_general(num_I = 5, num_T = 6, num_N = 84,
                    ICC = 0.03, Var_Y = 81, design = design_example, 
                    theta = 1.78, alpha = 0.05)

sample_size_model1_general(num_I = 5, num_T = 6,
                     ICC = 0.03, Var_Y = 81, design = design_example, 
                     theta = 1.78, alpha = 0.05, power = 0.8)


########################################################################
# Power of fixed effects model, time effects included, c clusters randomized
# per time period starting with period 2 (standard SWD)
########################################################################
power_model1_standard <- function(num_T, num_N, c, ICC, Var_Y, theta, alpha) {
  # Assumptions:
  
  # Analysis model: Model 1
  # A fixed number of clusters, c, cross over per time step (starting with period 2)
  
  # num_T     : number of periods
  # num_N     : sampled per cluster-period
  # c         : clusters crossed over per time period beginning in period 2
  # ICC       : Intraclass correlation coefficient
  # Var_Y     : total variance, sigma_e^2 = Var_Y*(1-ICC)
  # theta     : intervention effect
  # alpha     : type 1 error rate
  
  stopifnot(
    num_I > 0,
    num_T > 0,
    num_N > 0,
    ICC >= 0,
    ICC < 1,
    Var_Y > 0,
    alpha > 0,
    alpha < 1
  )
  
  z_alpha <- qnorm(1-alpha/2)
  
  # Calculate design constants
  var_theta <- Var_Y*(1-ICC)*12/(num_N*c*(num_T+1)*(num_T -2))
  power     <- pnorm(abs(theta)/sqrt(var_theta) -z_alpha) + pnorm(-abs(theta)/sqrt(var_theta) -z_alpha)
  
  return(power)
}

########################################################################
# Required sample size of fixed effects model, time effects included, c clusters randomized
# per time period starting with period 2 (standard SWD)
########################################################################
sample_size_model1_standard <- function(num_T, c, ICC, Var_Y, theta, alpha, power) {
  # Returns the per cluster period sample size required to achieve the desired power
  
  # Assumptions:
  # Analysis model: Model 1
  # A fixed number of clusters, c, cross over per time period (starting with period 2)
  
  # num_T     : number of periods
  # c         : clusters crossed over per time period beginning in period 2
  # ICC       : Intraclass correlation coefficient
  # Var_Y     : total variance, sigma_e^2 = Var_Y*(1-ICC)
  # theta     : intervention effect
  # alpha     : type 1 error rate
  # power     : desired power, as fraction i.e. 0.80
  
  stopifnot(
    c >= 1,
    num_T > 0,
    ICC >= 0,
    ICC < 1,
    Var_Y > 0,
    alpha > 0,
    alpha < 1,
    power > 0,
    power < 1
  )
  
  z_alpha <- qnorm(1-alpha/2)
  z_beta  <- qnorm(power)
  
  # Calculate design constants
  N <- Var_Y*(1-ICC)*12/(c*(num_T+1)*(num_T -2))*((z_alpha + z_beta)/theta)^2
  
  return(ceiling(N))
}

# Example 2:
# Real world examples based on PROMPT Study
# Model 1, no transition period included, assuming 5 clusters and 6 time periods, c = 1
power_model1_standard( num_T = 6, c= 1, num_N = 84,
                     ICC = 0.03, Var_Y = 81, 
                     theta = 1.78, alpha = 0.05)

sample_size_model1_standard( num_T = 6, c= 1,
                       ICC = 0.03, Var_Y = 81, 
                       theta = 1.78, alpha = 0.05, power = 0.8)



########################################################################
# Power of fixed effects model for SWD, without time effects (Model 2)
########################################################################
power_model2_general <- function(num_I, num_T, num_N, ICC, Var_Y, design, theta, alpha) {
  # Assumptions:
  
  # Analysis model: Model 2
  
  # num_I     : number of clusters
  # num_T     : number of periods
  # num_N     : number of observations sampled per cluster-period
  # ICC       : Intraclass correlation coefficient
  # Var_Y     : total variance, sigma_e^2 = Var_Y*(1-ICC)
  # design    : Matrix specifying the stepped wedge treatment allocation scheme. 
  #                   Rows correspond to clusters and columns correspond to time periods.
  #                   Entries are coded as 0 for control and 1 for intervention.
  # theta     : intervention effect
  # alpha     : type 1 error rate
  
  stopifnot(
    num_I > 0,
    num_T > 0,
    num_N > 0,
    ICC >= 0,
    ICC < 1,
    Var_Y > 0,
    alpha > 0,
    alpha < 1
  )
  
  z_alpha <- qnorm(1-alpha/2)
  
  U  <- sum(design)              # SUM_i (SUM_j X_{ij})
  W  <- sum(colSums(design)^2)   # SUM_j (SUM_i X_ij^2)
  V  <- sum(rowSums(design)^2)   # SUM_i (SUM_j X_ij^2)
  
  # Calculate design constants
  var_theta <- Var_Y*(1-ICC)*num_T/(num_N*(num_T*U - V))
  power     <- pnorm(abs(theta)/sqrt(var_theta) -z_alpha) + pnorm(-abs(theta)/sqrt(var_theta) -z_alpha)
  
  return(power)
}


########################################################################
# Required sample size of fixed effects model for SWD, without time effects (Model 2)
########################################################################
sample_size_model2_general <- function(num_I, num_T, ICC, Var_Y, design, theta, alpha, power) {
  # Returns the per cluster period sample size required to achieve the desired power
  
  # Assumptions:
  
  # Analysis model: Model 2
  
  # num_I     : number of clusters
  # num_T     : number of periods
  # ICC       : Intraclass correlation coefficient
  # Var_Y     : total variance, sigma_e^2 = Var_Y*(1-ICC)
  # design    : Matrix specifying the stepped wedge treatment allocation scheme. 
  #                   Rows correspond to clusters and columns correspond to time periods.
  #                   Entries are coded as 0 for control and 1 for intervention.
  # theta     : intervention effect
  # alpha     : type 1 error rate
  # power     : desired power, as fraction i.e. 0.80
  
  stopifnot(
    num_I > 0,
    num_T > 0,
    num_N > 0,
    ICC >= 0,
    ICC < 1,
    Var_Y > 0,
    alpha > 0,
    alpha < 1,
    power > 0,
    power < 1
  )
  
  z_alpha <- qnorm(1-alpha/2)
  z_beta  <- qnorm(power)
  
  U  <- sum(design)              # SUM_i (SUM_j X_{ij})
  W  <- sum(colSums(design)^2)   # SUM_j (SUM_i X_ij^2)
  V  <- sum(rowSums(design)^2)   # SUM_i (SUM_j X_ij^2)
  
  # Calculate design constants
  N <- Var_Y*(1-ICC)*num_T/(num_T*U - V)*((z_alpha + z_beta)/theta)^2
  
  return(ceiling(N))
}

# Example 3:
# Real world examples based on PROMPT Study
# Model 2, no transition period included, assuming 5 clusters and 6 time periods
design_example <- matrix(c(0,1,1,1,1,1,
                           0,0,1,1,1,1,
                           0,0,0,1,1,1,
                           0,0,0,0,1,1,
                           0,0,0,0,0,1), 
                         nrow = 5, ncol = 6, byrow =TRUE)

power_model2_general(num_I = 5, num_T = 6, num_N = 34,
                     ICC = 0.03, Var_Y = 81, design = design_example, 
                     theta = 1.78, alpha = 0.05)

sample_size_model2_general(num_I = 5, num_T = 6,
                 ICC = 0.03, Var_Y = 81, design = design_example, 
                 theta = 1.78, alpha = 0.05, power = 0.8)


########################################################################
# Power of fixed effects model for SWD, with transition period, 
# with time effects (Model 1)
########################################################################
power_model1_transition <- function( num_I, num_N, ICC, Var_Y, theta, alpha ) {
  # Assumptions:
  
  # Analysis model: Model 1
  # SWD with 1 cluster crossed over per time period beginning in period 3
  # and a transition period between control and treatment conditions
  # I > 2, I = T - 2
  
  # num_I     : number of clusters
  # num_N     : number of observations sampled per cluster-period
  # ICC       : Intraclass correlation coefficient
  # Var_Y     : total variance, sigma_e^2 = Var_Y*(1-ICC)
  # theta     : intervention effect
  # alpha     : type 1 error rate
  
  stopifnot(
    num_I > 2,
    num_N > 0,
    ICC >= 0,
    ICC < 1,
    Var_Y > 0,
    alpha > 0,
    alpha < 1
  )
  
  
  z_alpha <- qnorm(1-alpha/2)
  num_T <- num_I + 2
  
  # X[i, j]   : 1 if cluster i is in treatment at time step j, else 0
  X <- outer(seq_len(num_I), seq_len(num_T), function(i, j) +(j > i + 1))
  
  total_sum <- sum(X)
  col_sums  <- colSums(X)   # length num_T, col_sums[j] = SUM_i X_ij
  row_sums  <- rowSums(X)   # length num_I, row_sums[i] = SUM_j X_ij
  
  # interior_cols columns: indices 2..(T-1)
  interior_cols <- 2:(num_T - 1)
  
  interior_cols_total <- sum(col_sums[interior_cols])
  # sum of X over all interior_cols columns except d;
  sum_excluding_column <- function(d) interior_cols_total - col_sums[d]
  
  # denominator term: (T - 2)(T^2 - 4T + 2)
  denom <- (num_T - 2) * (num_T^2 - 4*num_T + 2)
  
  #calculation:
  #  c1 * T1 + c2 * T2
  c1 <- 1/(num_T - 3) + (num_T - 3)*(num_T - 1)/denom
  T1 <- sum(col_sums[interior_cols]^2)   # Sum_j (Sum_i X_ij)^2
  c2 <- (2*num_T^2 - 10*num_T + 11)/denom
  T2 <- sum(row_sums^2)             # Sum_i (Sum_j X_ij)^2

  # c3 * (T3 - 2*T4)
  c3 <- (num_T^2 - 5*num_T + 5)/denom
  T3 <- sum(row_sums * (total_sum - row_sums))       # SUM_q row_q * SUM_{i!=q} row_i
  T4 <- sum(row_sums * col_sums[seq_len(num_I) + 1])    # SUM_s row_s * col_{s+1}
  
  # c4 * (T5 + T6)
  c4 <- -(num_T - 3)*(num_T - 1)/denom
  T5 <- sum(col_sums[interior_cols] * (total_sum - row_sums[interior_cols - 1]))  # SUM_q col_q * SUM_{i!=q-1} row_i
  excl_col <- interior_cols_total - col_sums[seq_len(num_I) + 1]
  T6 <- sum(row_sums * excl_col)                     # SUM_s row_s * SUM_{j!=s+1} col_j
  
  # c5 * T7
  c5 <- (num_T^3 - 7*num_T^2 + 14*num_T - 7) / (denom * (num_T - 3))
  T7 <- sum(col_sums[seq_len(num_T - 1)] * (interior_cols_total - col_sums[seq_len(num_T - 1)]))
  
  m <- c1 * T1 + c2 * T2 +  c3 * (T3 - 2*T4) + c4 * (T5 + T6) + c5 * T7

  
  var_theta <-   Var_Y*(1-ICC)/ num_N/(total_sum - m)
  power     <- pnorm(abs(theta)/sqrt(var_theta) -z_alpha) + pnorm(-abs(theta)/sqrt(var_theta) -z_alpha)
  
  return(power)
}



########################################################################
# Required sample size of the fixed effects model for SWD, with transition period, 
# with time effects (Model 1)
########################################################################
sample_size_model1_transition <- function( num_I, ICC, Var_Y, theta, alpha, power ) {
  # Returns: (per cluster-period sample size, and total sample size)
  # Assumptions:
  
  # Analysis model: Model 1
  # SWD with 1 cluster crossed over per time period beginning in period 3
  # and a transition period between control and treatment conditions
  # I > 2, I = T - 2

  # num_I     : number of clusters
  # ICC       : Intraclass correlation coefficient
  # Var_Y     : total variance, sigma_e^2 = Var_Y*(1-ICC)
  # theta     : intervention effect
  # alpha     : type 1 error rate
  # power     : desired power, as fraction i.e. 0.80
  
  stopifnot(
    num_I > 2,
    ICC >= 0,
    ICC < 1,
    Var_Y > 0,
    alpha > 0,
    alpha < 1,
    power > 0,
    power < 1
  )
  
  z_alpha <- qnorm(1-alpha/2)
  z_beta  <- qnorm(power)
  
  num_T <- num_I + 2
  
  # X[i, j]   : 1 if cluster i is in treatment at time step j, else 0
  X <- outer(seq_len(num_I), seq_len(num_T), function(i, j) +(j > i + 1))
  
  total_sum <- sum(X)
  col_sums  <- colSums(X)   # length num_T, col_sums[j] = SUM_i X_ij
  row_sums  <- rowSums(X)   # length num_I, row_sums[i] = SUM_j X_ij
  
  # interior_cols columns: indices 2..(T-1)
  interior_cols <- 2:(num_T - 1)
  
  interior_cols_total <- sum(col_sums[interior_cols])
  # sum of X over all interior_cols columns except d;
  sum_excluding_column <- function(d) interior_cols_total - col_sums[d]
  
  # denominator term: (T - 2)(T^2 - 4T + 2)
  denom <- (num_T - 2) * (num_T^2 - 4*num_T + 2)
  
  #calculation:
  #  c1 * T1 + c2 * T2
  c1 <- 1/(num_T - 3) + (num_T - 3)*(num_T - 1)/denom
  T1 <- sum(col_sums[interior_cols]^2)   # Sum_j (Sum_i X_ij)^2
  c2 <- (2*num_T^2 - 10*num_T + 11)/denom
  T2 <- sum(row_sums^2)             # Sum_i (Sum_j X_ij)^2
  
  # c3 * (T3 - 2*T4)
  c3 <- (num_T^2 - 5*num_T + 5)/denom
  T3 <- sum(row_sums * (total_sum - row_sums))       # SUM_q row_q * SUM_{i!=q} row_i
  T4 <- sum(row_sums * col_sums[seq_len(num_I) + 1])    # SUM_s row_s * col_{s+1}
  
  # c4 * (T5 + T6)
  c4 <- -(num_T - 3)*(num_T - 1)/denom
  T5 <- sum(col_sums[interior_cols] * (total_sum - row_sums[interior_cols - 1]))  # SUM_q col_q * SUM_{i!=q-1} row_i
  excl_col <- interior_cols_total - col_sums[seq_len(num_I) + 1]
  T6 <- sum(row_sums * excl_col)                     # SUM_s row_s * SUM_{j!=s+1} col_j
  
  # c5 * T7
  c5 <- (num_T^3 - 7*num_T^2 + 14*num_T - 7) / (denom * (num_T - 3))
  T7 <- sum(col_sums[seq_len(num_T - 1)] * (interior_cols_total - col_sums[seq_len(num_T - 1)]))
  
  m <- c1 * T1 + c2 * T2 +  c3 * (T3 - 2*T4) + c4 * (T5 + T6) + c5 * T7
  

  num_N <-   Var_Y*(1-ICC)/(total_sum - m)*((z_alpha + z_beta)/theta)^2
  
  return(c(ceiling(num_N), ceiling(num_N)*(num_T-1)*num_I))
}


# Example 4:
# Real world examples based on PROMPT Study
# Model 1, transition period included, assuming 5 clusters and 7 time periods
power_model1_transition( num_I = 5, num_N = 128, 
                         ICC = 0.03, Var_Y = 81,
                         theta = 1.78, alpha  = 0.05)
sample_size_model1_transition( num_I = 5,  
                         ICC = 0.03, Var_Y = 81,
                         theta = 1.78, alpha  = 0.05, power = 0.8)
