
#' Function to solve counterfactuals.
#'
#' @param N Integer - Number of locations.
#' @param L_i Nx1 array - Number of residents in each location
#' @param L_j Nx1 array - Number of workers in each location
#' @param K Nx1 array - Land supply
#' @param t_ij NxN matrix - Travel times across locations
#' @param a Nx1 array - Total Factor Productivity in each location
#' @param b Nx1 array - Vector of amenities in each location
#' @param varphi Nx1 array - Density of development
#' @param w_eq Nx1 array - Initial vector of wages
#' @param u_eq Nx1 array - Initial vector of welfare
#' @param U_eq 1x1 array - Initial  welfare
#' @param Q_eq Nx1 array - Initial price for floorspace
#' @param ttheta_eq Nx1 array - Share of floorspace used commercially 
#' @param alpha Float - Exp. share in consumption, 1-alpha exp. share in housing
#' @param beta Float - Output elasticity with respect to labor
#' @param theta Float - Commuting and migration elasticity.
#' @param mu Float - Floorspace prod function: output elasticity wrt capital
#' @param delta Float - Decay parameter agglomeration force
#' @param lambda Float - agglomeration externality
#' @param rho Float - decay parameter for amenities
#' @param eta Float - amenity externality
#' @param epsilon Float - Parameter that transforms travel times to commuting costs
#' @param zeta Float - convergence parameter
#' @param tol Float - tolerance factor
#' @param varrho Float - migration elasticty
#' @param sh_city Share of population living in the city
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=1000.
#' 
#' @return A list containing new equilibrium values for the variables of the model. 
#' The reported variables are wages (w), residents market access (W_i), endogenous productivities (A), 
#' endogenous amenities (B), population (Li), employment (Lj), share of commercial floorspace (ttheta),
#' real income (u), welfare (U), firms market access (FCMA), firms output (Y). 

#' @export
#'

solveModel = function(N,
                      L_i,
                      L_j,
                      K,
                      t_ij,
                      a,
                      b,
                      varphi,
                      w_eq,
                      u_eq,
                      U_eq, 
                      Q_eq,
                      ttheta_eq,
                      alpha=0.7,
                      beta=0.7,
                      theta=7,
                      mu=0.3,
                      delta=0.3585,
                      lambda=0.01,
                      rho=0.9094,
                      eta=0.15,
                      epsilon=0.01,
                      zeta=0.95,
                      varrho=3,
                      sh_city = 0.1,
                      tol=10^-10,
                      maxiter=1000){
  
  # Formatting of input data
  if(is.data.frame(L_i)){
    L_i = array(unlist(L_i), dim(L_i))
  } else if(is.null(dim(L_i))){
    L_i = array(L_i, dim=c(N,1))
  }
  
  if(is.data.frame(L_j)){
    L_j = array(unlist(L_j), dim(L_j))
  } else if(is.null(dim(L_j))){
    L_j = array(L_j, dim=c(N,1))
  }
  if(is.data.frame(K)){
    K = array(unlist(K), dim(K))  
  } else if(is.null(dim(K))){
    K = array(K, dim=c(N,1))
  }
  
  t_ij = array(unlist(t_ij), dim(t_ij))  
  
  if(is.null(dim(a))){
    a = array(a, dim=c(N,1))
  }
  if(is.null(dim(b))){
    b = array(b, dim=c(N,1))
  }
  if(is.null(dim(varphi))){
    varphi = array(varphi, dim=c(N,1))
  }
  if(is.null(dim(w_eq))){
    w_eq = array(w_eq, dim=c(N,1))
  }
  if(is.null(dim(u_eq))){
    u_eq = array(u_eq, dim=c(N,1))
  }
  
  if(is.null(dim(U_eq))){
    U_eq = array(U_eq, dim=c(1,1))
  }
  
  if(is.null(dim(Q_eq))){
    Q_eq = array(Q_eq, dim=c(N,1))
  }
  if(is.null(dim(ttheta_eq))){
    ttheta_eq = array(ttheta_eq, dim=c(N,1))
  }
  
  # Normalize L_i to have the same size as L_j
  L_i=L_i*sum(L_j)/sum(L_i)
  L_bar=sum(L_i)
  L_bar_init = L_bar
  lambda_i = array_operator(L_i, L_bar, '/')
  
  D = commuting_matrix(t_ij=t_ij, epsilon = epsilon)
  tau = D$tau
  L_i = array(unlist(L_i),dim(L_i))
  L_j = array(unlist(L_j),dim(L_j))
  K = array(unlist(K), dim(K))
  
  # Settings
  outerdiff = Inf;
  w = w_eq;
  u = u_eq;
  U_init = U_eq;
  Q = Q_eq;
  Vo_varrho  = ((U_eq^varrho)*(1-sh_city))/(sh_city)
  ttheta = ttheta_eq;
  iter = 0;
  zeta_init = zeta;
  
  cat("Solving model...\n")
  
  while(outerdiff>tol & iter < maxiter){
    # 1) Labor supply equation
    w_tr = aperm(array(w, dim=c(N,1)), c(2,1));
    rep_w_tr = kronecker(w_tr^theta, array(1, dim=c(N, 1)));
    # Constructing employment shares
    w_tr_tau = array_operator(w_tr^theta, tau^(-theta), '*');
    lambda_ij_i = array_operator(w_tr_tau, sumDims2(w_tr_tau,2), '/');
    
    W_i = (sumDims2(w_tr_tau,2))^(1/theta);
    
    # Labor is equal to probabilities * total number of residents * proportion of workers in each sector.
    L_ij = array_operator(L_i, lambda_ij_i, '*')
    L_j = sumDims2(L_ij, 1)
    
    
    # 2 average income
    av_income = av_income_simple(lambda_ij_i=lambda_ij_i,w_tr = w_tr)
    ybar = av_income$y_bar
    
    # 3 Total floorspace
    FS = array_operator(varphi,K^(1-mu),"*")
    
    # 4 Agglomeration externalities
    L_j_dens = (array_operator(L_j, K, '/'));
    L_j_dens_per = aperm(array(L_j_dens, dim=c(N,1)), c(2,1));
    L_j_dens_rep = kronecker(L_j_dens_per, array(1, dim=c(N, 1)));
    Upsilon = sumDims2(array_operator(exp(-delta*t_ij), L_j_dens_rep, '*'), 2);    
    A = array_operator(a, Upsilon^lambda, '*')
    
    # 5 Amenities
    L_i_dens = (array_operator(L_i, K, '/'));
    L_i_dens_per = aperm(array(L_i_dens, dim=c(N,1)), c(2,1));
    L_i_dens_rep = kronecker(L_i_dens_per, array(1, dim=c(N, 1)));
    Omega = sumDims2(array_operator(exp(-rho*t_ij), L_i_dens_rep, '*'), 2);
    B = array_operator(b, Omega^(eta),'*')
    
    # 6 Residents, probabilities, and welfare
    u =  array_operator(array_operator(W_i, Q^(1-alpha), '/'), B, '*')
    U = sum(u^theta)
    lambda_i_upd = (u^theta)/U
    U = U^(1/theta)
    
    # 7 Total output by location
    FS_f = array_operator(ttheta,array_operator(varphi, K^(1-mu), '*'), '*')
    Y = array_operator(A, array_operator(L_j^beta, FS_f^(1-beta), '*'), '*')
    Q_upd1 = (1-beta)*array_operator(Y,FS_f, '/')
    w_upd = beta*array_operator(Y, L_j, '/')
    
    # 8 Housing prices
    FS_r = array_operator((1-ttheta), array_operator(varphi, K^(1-mu), '*'), '*')
    X = array_operator(ybar, L_i, '*')
    Q_upd2 = (1-alpha)*array_operator(X, FS_r, '/')
    Q_upd = Q_upd1*(a>0) + Q_upd2*(a==0 & b>0)
    
    # 9 Share of commercial floorspace
    #LP = array_operator(Q_upd1, array_operator(varphi, K^(1-mu), '*'), '*')
    #ttheta_upd = (1-beta)*array_operator(Y, LP, '/')
    #ttheta_upd = (b==0)+ttheta_upd*(b>0)
    FS_f_upd    = (1-beta)*array_operator(Y,Q_upd, '/')
    FS_r_upd    = (1-alpha)*array_operator(X,Q_upd, '/')
    ttheta_upd  = array_operator(FS_f_upd, FS_f_upd+FS_r_upd,'/')

    # 10 Calculating the main differences
    z_w = array_operator(w, w_upd, '-')
    z_L = array_operator(lambda_i, lambda_i_upd, '-')
    z_Q = array_operator(Q, Q_upd, '-')
    z_theta = array_operator(ttheta, ttheta_upd, '-')
    #z_U = array_operator(U, U_init, '-')
    
    
    if (varrho>0){
      sh_city_upd = (U^(varrho))/((U^(varrho))+Vo_varrho)
      Lhat        = sh_city_upd/sh_city
      L_bar_upd   = Lhat*L_bar_init
      z_L_bar     = L_bar-L_bar_upd
      L_bar       = zeta*L_bar_upd + (1-zeta)*L_bar
    }
    
    if (varrho==0){
      sh_city_upd = sh_city
      L_bar       = L_bar
      Lhat        = 1
      z_L_bar     = 0
    }
    
    #outerdiff = max(c(max(abs(z_w)), max(abs(z_L)), max(abs(z_Q)), max(abs(z_theta)), max(abs(z_U)), abs(z_L_bar)))
    outerdiff = max(c(max(abs(z_w)), max(abs(z_L)), max(abs(z_Q)),  max(abs(z_theta)), abs(z_L_bar)))
    
    #Lbar_upd  = L_bar*(1 + 0.25*(array_operator(array_operator(U, U_init, '-'), U_init, '/'))); 
    
    iter = iter+1
    
    # 11 New vector of variables
    lambda_i = zeta*lambda_i + (1-zeta)*lambda_i_upd
    Q = zeta*Q + (1-zeta)*Q_upd
    w = zeta*w + (1-zeta)*w_upd
    ttheta = zeta*ttheta + (1-zeta)*ttheta_upd
    
    #L_bar = zeta*Lbar_upd + (1-zeta)*L_bar
    L_i = array_operator(L_bar, lambda_i, '*')
    
    
    #12 Firms market access
    WA = array_operator(L_i, tau^(-theta), '*');
    FCMA = array_operator(WA, W_i, '/');
    FCMA = (sumDims2(FCMA,2))^(1/theta);   
    
    
    if(iter %% 25 == 0){
      cat(paste0("Iteration: ", iter, ", error: ", round(outerdiff, 10), ".\n"))
    }
  }
  if(outerdiff<=tol){
    cat(paste0("Converged after ", iter, " iterations. Error=", round(outerdiff, 10), ".\n"))
  } else{
    cat(paste0("Reached maximum number of iterations (", iter, "). Error=", round(outerdiff, 10), ".\n"))
  }
  
  return(list(w=w, W_i=W_i, B=B, A=A, Q=Q, L_i=L_i, L_j=L_j,
              ybar=ybar, lambda_i=lambda_i, ttheta=ttheta, u=u, U=U, L_bar=L_bar, FCMA=FCMA, Y=Y))
}