#' Function to compute the equilibrium residential and commercial floorspace and equilibrium price.
#'
#' @param K Nx1 array. Land supply.
#' @param w Nx1 array. Wages in each location.
#' @param L_j Nx1 array. Number of workers in each location.
#' @param y_bar Nx1 array. Average income in each location.
#' @param L_i Nx1 array. Number of residents in each location.
#' @param beta Float. Output elasticity with respect to labor.
#' @param alpha Float. Expenditure share on consumption (1-alpha on housing).
#' @param mu Float. Floorspace production elasticity: F = M^mu * K^(1−mu).
#' @param FS Nx1 array. Total floorspace.
#'
#' @return A list containing the equilibrium values for the floorspace market. 
#' It returns the geometric mean floorspace price (Q_mean), the normalized floorspace price (Q_norm),
#' the commercial floorspace (FS_f), the residential floorspace (FS_r) and the density of development (varphi).
#'   
density_development_noprices = function(K,
                               w,
                               L_j,
                               y_bar,
                               L_i,
                               beta,
                               alpha,
                               mu,
                               FS) {
  
  # Income component
  X = array_operator(y_bar, L_i, '*')
  
  # Total expenditure on floorspace: residential + commercial
  D = (1 - alpha) * X +
    ((1 - beta)/beta ) * array_operator(w, L_j, '*')
  
  # Development density
  varphi = array_operator(FS, K^(1 - mu), '/')
  
  #Floorspace prices
  Q = array_operator(D, FS, '/')
  Q_mean = exp(mean(log(Q)))
  Q_norm = Q
  
  FS_f   = array_operator(((1 - beta)/beta ) * array_operator(w, L_j, '*'), Q, '/')
  FS_r   = array_operator((1 - alpha) * X, Q, '/')
  
  return(list(
    Q_mean = Q_mean,
    Q_norm = Q_norm,     
    varphi = varphi,
    FS_f   = FS_f,
    FS_r   = FS_r,
    FS     = FS
  ))
}