library(rdist)
library(MASS)

#### Helper Functions

find_second_smallest <- function (m_arr) {
  smallest_ind <- 0
  second_ind <- 0
  smallest <- Inf
  second <- Inf
  
  for (i in 1:length(m_arr)) {
    if (m_arr[i] < smallest) {
      second_ind <- smallest_ind
      second <- smallest
      smallest_ind <- i
      smallest <- m_arr[i]
    } else if (m_arr[i] < second) {
      second_ind <- i
      second <- m_arr[i]
    } 
  } 
  return(second_ind)
}

find_smallest_non_zero <- function (m_arr) {
  smallest_ind <- 0
  smallest <- Inf
  for (i in 1:length(m_arr)) {
    if (m_arr[i] < smallest) {
      if(m_arr[i] > 0) {
        smallest <- m_arr[i]
        smallest_ind <- i
      }
    }
  } 
  return(smallest_ind)
}

row_min_inds <- function( m_mat, print_flag = FALSE ) {
  output <- rep(0, dim(m_mat)[1])
  for (i in 1:dim(m_mat)[1]) {
    output[i] <- find_smallest_non_zero(m_mat[i,])
  }
  if (print_flag) {
    counts <- sort(table( output ), decreasing=T )
    print(counts)
  }
  return(output)
}

select_from_ind <- function( m_mat, m_inds ) {
  output <- rep(0, dim(m_mat)[1])
  for (i in 1:dim(m_mat)[1] ) {
    output[i] <- m_mat[i, m_inds[i]]
  }
  return(output)
}


#### Testing Functions

AsymVarTest <- function( m_X, m_Y, d_X, norm_Y, m_V, m_mu_g, m_p_t, m_thresh, m_m, m_g_dot, m_g_star) {
  n <- dim(m_X)[1]
  dist_mat <- matrix(0, nrow = m_m , ncol = n)
  dif_mat <- matrix(0, nrow = m_m, ncol = n)
  
  for (i in 1:m_m) {
    g_i <- sample( m_mu_g ,1)
    ind <- sample(1:n, 1)
    X_i <- m_X[ind,]
    Y_i <- m_Y[ind,]
    g_cdot_X_i <- m_g_dot(g_i, X_i)
    g_star_Y_i <- m_g_star(g_i, Y_i)
    dist_mat[i,] <- d_X( g_cdot_X_i, m_X )
    dif_mat[i,] <- norm_Y( g_star_Y_i, m_Y)
  }
  
  inds <- row_min_inds(dist_mat)
  dists <- select_from_ind(dist_mat, inds)
  difs <- select_from_ind(dif_mat, inds)
  
  p_val <- pbinom(sum(I( m_V(dists) < difs - m_thresh) ) - 1, 
                  sum(I(difs>0)), m_p_t(m_thresh), lower.tail = FALSE )
  
  return(p_val)
}

PermVarTest <- function( m_X, m_Y, d_X, norm_Y, m_V_Cal, m_q, m_mu_g, m_m, m_B, m_g_dot, m_g_star ) {
  n <- dim(m_X)[1] 
  ratios <- rep(0, m_B)
  for (j in 1:m_B) {
    dist_mat <- matrix(0, nrow = m_m , ncol = n)
    dif_mat <- matrix(0, nrow = m_m, ncol = n)
    for (k in 1:m_m) {
      g_i <- sample(m_mu_g, 1)
      ind <- sample(1:n, 1)
      X_i <- m_X[ind,]
      Y_i <- m_Y[ind,]
      g_dot_X_i <- m_g_dot( g_i, X_i )
      g_star_Y_i <- m_g_star( g_i , Y_i )
      dist_mat[k,] <- d_X( g_dot_X_i, m_X )
      dif_mat[k,] <- norm_Y(g_star_Y_i, m_Y)
    }
    
    rats <- dif_mat[which(dist_mat != 0)] / dist_mat[which(dist_mat != 0)]
    ratios[j] <- quantile(rats, probs = m_q, type = 1, names = FALSE)  }
  
  dist_mat_orig <- matrix(0, nrow = m_m , ncol = n - 1)
  dif_mat_orig <- matrix(0, nrow = m_m, ncol = n - 1)
  
  for (k in 1:m) {
    ind <- sample(1:n, 1)
    X_i <- m_X[ind,]
    Y_i <- m_Y[ind,]
    dist_mat_orig[k,] <- d_X( X_i, m_X[-ind, ] )   # Important to exclude the original data point so we don't just see 0/0
    dif_mat_orig[k,] <- norm_Y( Y_i, m_Y[-ind,] )
  }
  
  orig_rats <- dif_mat_orig[which(dist_mat_orig != 0)] / dist_mat_orig[which(dist_mat_orig != 0)] ## FIX FOR NEAREST NEIGHBOURS
  orig_ratio <- quantile(orig_rats, probs = m_q, type = 1, names = FALSE)
  
  p_val <- mean( I(ratios <= orig_ratio) )
  return(p_val)
}

#### Simulations

## Set-up Pars

g_dot <- function( g, X ) {
  X_new <- X
  X_new[1] <- (X[1] * cos(g)) - (X[2] * sin(g))
  X_new[2] <- (X[1] * sin(g)) + (X[2] * cos(g))
  return( X_new )
}

g_star <- function( g, X ) {
  return( g_dot( 2 * g, X ) )
}

g_bullet <- function( g, X ) {
  return( X )
} 

d_X <- function ( x_i, X ) {
  return( sqrt( colSums( (x_i - t(X) )^2 ) ) )
}

norm_Y <- function ( y_i, Y ) {
  return( sqrt( colSums( (y_i - t(Y) )^2 ) ) )
}

a_V <- function( d ) {
  return( d )
}

a_V_e <- function( d ) {
  return( exp(-1) * d )
}

a_V_cal <- function( d ) {
  return( d )
}

f_d <- function( X ) {
  norm <- sqrt( X[,1]^2 )
  return( exp( -norm ) )
}

a_p_t <- function ( m_thres ) {
  return( 2 * exp( - (m_thres / 0.05) ^2 / 4) / ( (m_thres / 0.05) * sqrt(2 * pi )) )
}

#### f_d Sims
set.seed(1010)
num_sims = 100
alpha <- 0.05
ns <- c(20,30,40,50,60,70,80,90,100,120,150,200,250,300)
d <- 4
mu_X <- rep(0, d)
Sigma_X = diag(rep(2,d))

sigma <- 0.05
mu_eps <- rep(0, d)
Sigma_eps = diag(rep(sigma^2,d))

G <- c( pi/2, pi, 3*pi/2 )

a_q <- 0.95
a_B <- 100

rejections_fd <- matrix(0, ncol = length(ns), nrow = 4)
pb <- txtProgressBar(min = 0, max = sum(ns) * num_sims, style = 3, width = 50, char = "=") 
counter <- 0

for (k in 1:length(ns)) {
  n <- ns[k]
  m <- ns[k]
  p_vals_fd_H0 <- rep(0, num_sims)
  p_vals_fd_H1 <- rep(0, num_sims)  
  p_vals_fd_H0_P <- rep(0, num_sims)
  p_vals_fd_H1_P <- rep(0, num_sims)
  
  for (i in 1:num_sims) {
    X <- mvrnorm(n, mu_X, Sigma_X)
    Y <- as.matrix( f_d(X) + rnorm(n, 0, sigma) )
    p_vals_fd_H0[i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, G, a_p_t, 2*sigma, m, g_star, g_bullet)
    p_vals_fd_H1[i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, G, a_p_t, 2*sigma, m, g_dot, g_bullet)
    p_vals_fd_H0_P[i] <- PermVarTest(X, Y, d_X, norm_Y, a_V_cal, a_q, G, m, a_B, g_star, g_bullet)
    p_vals_fd_H1_P[i] <- PermVarTest(X, Y, d_X, norm_Y, a_V_cal, a_q, G, m, a_B, g_dot, g_bullet)
    counter <- counter + n
    setTxtProgressBar(pb, counter)
  }
  
  rejections_fd[1,k] <- sum(I(p_vals_fd_H0 < alpha)) / num_sims
  rejections_fd[2,k] <- sum(I(p_vals_fd_H1 < alpha)) / num_sims
  rejections_fd[3,k] <- sum(I(p_vals_fd_H0_P < alpha)) / num_sims
  rejections_fd[4,k] <- sum(I(p_vals_fd_H1_P < alpha)) / num_sims
  
  print(paste("--Rejections for n=", n, "--"))
  print(paste( rejections_fd[1,k], rejections_fd[2,k], rejections_fd[3,k], rejections_fd[4,k]) )
}



## Plots for f_d

plot(ns, rejections_fd[4,], type = "p", pch = 17, ylab  = "", xlab = "n", ylim = c(0,1), log = "x")
points(ns, rejections_fd[3,], pch = 6)
points(ns, rejections_fd[2,], pch = 15)
points(ns, rejections_fd[1,], pch = 22)

lines(ns, rep(0.05, 14), lty = 3)
lines(ns, rejections_fd[4,], lty = 2)
lines(ns, rejections_fd[2,], lty = 4)




#### Plots for Fig 3

set.seed(1011)
n <- 20
d <- 2
mu_X <- rep(0, d)
Sigma_X = diag(rep(2,d))

sigma <- 0.05
mu_eps <- rep(0, d)
Sigma_eps = diag(rep(sigma^2,d))

G <- c( pi/2, pi, 3*pi/2 )

X <- mvrnorm(n, mu_X, Sigma_X)
Y_2 <- as.matrix( f_d(X) + rnorm(n, 0, sigma) )

X_d <- pdist(X)
Y_d <- pdist(Y_2)


G_dot_X <- matrix( nrow = 3 * n, ncol = d )
for (j in 1:3) {
  for (i in 1:n) {
    G_dot_X[(j-1)*n + i,] <- g_dot( G[j], X[i,] )
  }
}
G_star_X <- matrix( nrow = 3 * n, ncol = d )
for (j in 1:3) {
  for (i in 1:n) {
    G_star_X[(j-1)*n + i,] <- g_star( G[j], X[i,] )
  }
}
G_Y <- rbind( Y_2, Y_2, Y_2 )

G_star_X_d <- pdist(G_star_X)
G_dot_X_d <- pdist(G_dot_X)
G_Y_d <- pdist(G_Y)

plot(G_dot_X_d, G_Y_d, type = "p", pch = 4, ylab = "", xlab = "")
points(X_d, Y_d, pch = 4, col = "red")
lines(c(0,8), c(0,8*exp(-1)), col ="red")

plot(G_star_X_d, G_Y_d, type = "p", pch = 4, ylab = "", xlab = "")
points(X_d, Y_d, pch = 4, col = "red")
lines(c(0,8), c(0,8*exp(-1)), col ="red")



#### Effect of changing L

set.seed(1010)
num_sims = 100
alpha <- 0.05
ns <- c(20,30,40,50,60,70,80,90,100,120,150,200,250,300)
d <- 2
mu_X <- rep(0, d)
Sigma_X = diag(rep(2,d))

sigma <- 0.05
mu_eps <- rep(0, d)
Sigma_eps = diag(rep(sigma^2,d))

G <- c( pi/2, pi, 3*pi/2 )

Ls <- c(exp(-3), exp(-2), exp(-1.2), exp(-1),  0.5, 1, 2)

rejections_fd_H0 <- matrix(0, ncol = length(ns), nrow = length(Ls))
rejections_fd_H1 <- matrix(0, ncol = length(ns), nrow = length(Ls))


for (j in 1:length(Ls)) {
  a_V_L <- function (d) { return( Ls[j] * d ) }
  pb <- txtProgressBar(min = 0, max = sum(ns) * num_sims, style = 3, width = 50, char = "=") 
  counter <- 0
    for (k in 1:length(ns)) {
    n <- ns[k]
    m <- ns[k]
    p_vals <- rep(0, num_sims)
    for (i in 1:num_sims) {
      X <- mvrnorm(n, mu_X, Sigma_X)
      Y <- as.matrix( f_d(X) + rnorm(n, 0, sigma) )
      p_vals_fd_H0[i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_L, G, a_p_t, 2*sigma, m, g_star, g_bullet)
      p_vals_fd_H1[i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_L, G, a_p_t, 2*sigma, m, g_dot, g_bullet)
      counter <- counter + n
      setTxtProgressBar(pb, counter)
    }
    rejections_fd_H0[j,k] <- sum(I(p_vals_fd_H0 < alpha)) / num_sims
    rejections_fd_H1[j,k] <- sum(I(p_vals_fd_H1 < alpha)) / num_sims
    }
  print(paste("--------", Ls[j], "--------"))
  print(rejections_fd_H0[j,])
  print(rejections_fd_H1[j,])
  
}

## Estimated Power Plots
plot(ns, rejections_fd_H1[4,], type = "p", pch = 17, ylab  = "", xlab = "n", ylim = c(0,1), log = "x")
points(ns, rejections_fd_H1[5,], pch = 6)
points(ns, rejections_fd_H1[6,], pch = 15)
points(ns, rejections_fd_H1[7,], pch = 22)

lines(ns, rejections_fd_H1[4,], lty = 1)
lines(ns, rejections_fd_H1[5,], lty = 2)
lines(ns, rejections_fd_H1[6,], lty = 3)
lines(ns, rejections_fd_H1[7,], lty = 4)

## Empirical Size Plots

plot(ns, rejections_fd_H0[4,], type = "p", pch = 17, ylab  = "", xlab = "n", ylim = c(0,1), log = "x")
points(ns, rejections_fd_H0[3,], pch = 6)
points(ns, rejections_fd_H0[2,], pch = 15)
points(ns, rejections_fd_H0[1,], pch = 22)

lines(ns, rejections_fd_H0[4,], lty = 1)
lines(ns, rejections_fd_H0[3,], lty = 2)
lines(ns, rejections_fd_H0[2,], lty = 3)
lines(ns, rejections_fd_H0[1,], lty = 4)
