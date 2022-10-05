library(rdist)
library(keras)

set.seed(1010)

#### Data
data.mnist <- dataset_mnist()
train_x<-data.mnist$train$x
train_y<-data.mnist$train$y
test_x<-data.mnist$test$x
test_y<-data.mnist$test$y
train_x <- array(as.numeric(train_x), dim = c(dim(train_x)[[1]], 784))
test_x <- array(as.numeric(test_x), dim = c(dim(test_x)[[1]], 784))
train_x <- train_x / 255
test_x <- test_x / 255
train_y<-to_categorical(train_y,10)
test_y<-to_categorical(test_y,10)


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



#### Estimating L

est_L <- function (digit_ind, sample_num = 0) {
  digit.inds <- which(train_y[,digit_ind] == 1, arr.ind = TRUE)
  non.digit.inds <- which(train_y[,digit_ind] == 0, arr.ind = TRUE)
  digit.x <- train_x[digit.inds,]
  non.digit.x <- train_x[sample(non.digit.inds, length(digit.inds)), ]
  if (sample_num == 0) {
    digits.dists.out <- cdist(digit.x, non.digit.x)
    return( 1 / min( digits.dists.out ) )
  } else {
    inds <- sample( 1:length(digit.inds), sample_num )
    digits.dists.out <- cdist(digit.x[inds,], non.digit.x[inds,])
    return( 1 / min( digits.dists.out ) )
  }
  return(0)
}

## These functions take approximately 7 minutes to run each. 

# system.time( zeros_L <-  est_L( 1, 5000 ) )
# print(zeros_L)
# system.time( ones_L <-  est_L( 2, 5000 ) )
# print(ones_L)
system.time( twos_L <- est_L( 3, 5000 ) )
print(twos_L) ## Calculated as 0.2694746
system.time( threes_L <-  est_L( 4, 5000 ) )
print(threes_L) ## Calculated as 0.2650288
system.time( fours_L <-  est_L( 5, 5000 ) )
print(fours_L) ## Calculated as 0.2914684
system.time( fives_L <-  est_L( 6, 5000 ) )
print(fives_L) ## Calculated as 0.2370678
system.time( sixs_L <-  est_L( 7, 5000 )  )
print(sixs_L) ## Calculated as 0.2569751
system.time( sevens_L <-  est_L( 8, 5000 ) )
print(sevens_L) ## Calculated as 0.3384539
# system.time( eights_L <-  est_L( 9, 5000 ) )
# print(eights_L)
system.time( nines_L <-  est_L( 10, 5000 ) )
print(nines_L) ## Calculated as 0.3534019


#### G = D_4
rotate = function(original, degree){
  # rotate 28x28 greyscale images by the given degrees clockwise
  rotated = original
  img = as.matrix(original[,-1])
  dim(img) <- c(1,784)
  ind = matrix(FALSE,40,40); ind[7:34,7:34] = TRUE; ind = as.vector(ind)
  tmp = matrix(0, dim(img)[1], 40*40); tmp[, ind] = img
  i = rep(7:34, 28); j = rep(7:34, each=28)
  a = degree/180*pi
  i_rotated = (i-20.5)*cos(a) + (j-20.5)*sin(a) + 20.5
  j_rotated = -(i-20.5)*sin(a) + (j-20.5)*cos(a) + 20.5
  i0 = floor(i_rotated); i1 = ceiling(i_rotated); lambda_i = i_rotated-i0
  j0 = floor(j_rotated); j1 = ceiling(j_rotated); lambda_j = j_rotated-j0
  rotated[,-1] = tmp[, i0+(j0-1)*40] * (1-lambda_i) * (1-lambda_j) + 
    tmp[, i0+(j1-1)*40] * (1-lambda_i) * lambda_j + 
    tmp[, i1+(j0-1)*40] * lambda_i * (1-lambda_j) + 
    tmp[, i1+(j1-1)*40] * lambda_i * lambda_j
  return(rotated)
}

rotate_set = function(original, degree){
  # rotate 28x28 greyscale images by the given degrees clockwise
  rotated = original
  img = as.matrix(original[,-1])
  ind = matrix(FALSE,40,40); ind[7:34,7:34] = TRUE; ind = as.vector(ind)
  tmp = matrix(0, dim(img)[1], 40*40); tmp[, ind] = img
  i = rep(7:34, 28); j = rep(7:34, each=28)
  a = degree/180*pi
  i_rotated = (i-20.5)*cos(a) + (j-20.5)*sin(a) + 20.5
  j_rotated = -(i-20.5)*sin(a) + (j-20.5)*cos(a) + 20.5
  i0 = floor(i_rotated); i1 = ceiling(i_rotated); lambda_i = i_rotated-i0
  j0 = floor(j_rotated); j1 = ceiling(j_rotated); lambda_j = j_rotated-j0
  rotated[,-1] = tmp[, i0+(j0-1)*40] * (1-lambda_i) * (1-lambda_j) + 
    tmp[, i0+(j1-1)*40] * (1-lambda_i) * lambda_j + 
    tmp[, i1+(j0-1)*40] * lambda_i * (1-lambda_j) + 
    tmp[, i1+(j1-1)*40] * lambda_i * lambda_j
  return(rotated)
}

reflect_v <- function(original) {
  # Reflect a 28x28 greyscale image through the vertical axis
  img = as.matrix(original[,-1], nrow = 28, ncol = 28)
  dim(img) <- c(28,28)
  reflected <- img[28:1,]
  dim(reflected) <- c(1,784)
  output <- cbind(as.matrix(original[1,1]), reflected)
  dim(output) <- c(1,785)
  return(output)
}

reflect_h <- function(original) {
  # Reflect a 28x28 greyscale image through the horizontal axis
  img = as.matrix(original[,-1], nrow = 28, ncol = 28)
  dim(img) <- c(28,28)
  reflected <- img[,28:1]
  dim(reflected) <- c(1,784)
  output <- cbind(as.matrix(original[1,1]), reflected)
  dim(output) <- c(1,785)
  return(output)
}

id  <- function(x) { return( x ) }
g_1 <- function(x) { return( rotate(x, 270) ) }  # Rotate 90 AC
g_2 <- function(x) { return( rotate(x, 180) ) } # Rotate 180 
g_3 <- function(x) { return( rotate(x, 90) ) }  # Rotate 270
g_4 <- function(x) { return( reflect_v(x)  ) }  # Reflect horizontally
g_5 <- function(x) { return( reflect_h(x) ) }  # Reflect vertically
g_6 <- function(x) { return( g_1(g_4(x)) ) }     # Reflect Down Diagonal
g_7 <- function(x) { return( g_3(g_4(x)) ) }     # Reflect Up Diagonal

D_4 <- c( id, g_1, g_2, g_3, g_4, g_5, g_6, g_7 )


#### Testing Functions

run_test <- function( X, Y, G, g_inds, t, p_t, V, m ) {
  n <- dim(X)[1]
  
  dist_mat <- matrix(0, nrow = m , ncol = n)
  dif_mat <- matrix(0, nrow = m, ncol = n)
  
  pb <- txtProgressBar(min = 0, max = m, style = 3, width = 50, char = "=")
  
  for (k in 1:m) {
    i <- sample(1:n, 1)
    X_i <- X[i,] 
    g_ind <- g_inds[sample.int(length(g_inds), size = 1)]
    X_prime <- as.matrix(G[[g_ind]](cbind(1, X_i))[,-1])
    dim(X_prime) <- c(1,784)
    dist_mat[k,] <- cdist(X_prime, X, metric = "euclidean", p = 2)
    dif_mat[k,] <- abs( Y[i] - Y )
    setTxtProgressBar(pb, k)
  }
  
  close(pb)
  
  inds <- row_min_inds(dist_mat)
  dists <- select_from_ind(dist_mat, inds)
  difs <- select_from_ind(dif_mat, inds)
  
  print(max(difs / dists))
  return( p_val(dists, difs, V, t, p_t) )
}

p_val <- function (dists, difs, V, t, p_t) {
  print(paste("N_t = ", sum(I( V(dists) < difs - t) ), " out of ", length(difs) ) )
  return( pbinom(sum(I( V(dists) < difs - t) ) - 1, 
                 length(difs), p_t, lower.tail = FALSE ) )
}


#### Running Tests

## Test Params
m <- 1000
t_mnist <- 0
p_t_mnist <- 0

test_digit <- function (digit_ind, G, g_inds) {
  digit.inds <- which(train_y[,digit_ind] == 1, arr.ind = TRUE)
  digit.x <- train_x[digit.inds,]
  dig.m <- floor(length(digit.inds)/2)
  X_Y <- matrix(0, nrow = dig.m, ncol = 785)
  digit.to.rot <- cbind(rep(1,dig.m), digit.x[1:dig.m,] )
  for (i in 1:dig.m) {
    X_Y[i,] <- reflect_h(t(as.matrix(digit.to.rot[i,])))
  }
  X <- rbind(X_Y[,-1], digit.x[(dig.m+1):length(digit.inds),])
  Y <- rep(1, dim(X)[1])
  Y[1:dig.m] <- 0    
  system.time(p <- run_test(X, Y, G, g_inds, t_mnist, p_t_mnist, V_mnist, m))
  return(p)
}

## NB: The digit.ind for n is n+1
set.seed(1010)

V_mnist <- function(d) {
  # Conservative Constraint from est_L
  return (d * twos_L )
}
test_digit(3, D_4, c(2,6))
test_digit(3, D_4, c(2))
test_digit(3, D_4, c(6))

V_mnist <- function(d) {
  # Conservative Constraint from est_L
  return (d * threes_L )
}
test_digit(4, D_4, c(2,6))
test_digit(4, D_4, c(2))
test_digit(4, D_4, c(6))

V_mnist <- function(d) {
  # Conservative Constraint from est_L
  return (d * fours_L )
}
test_digit(5, D_4, c(2,6))
test_digit(5, D_4, c(2))
test_digit(5, D_4, c(6))

V_mnist <- function(d) {
  # Conservative Constraint from est_L
  return (d * fives_L )
}
test_digit(6, D_4, c(2,6))
test_digit(6, D_4, c(2))
test_digit(6, D_4, c(6))

V_mnist <- function(d) {
  # Conservative Constraint from est_L
  return (d * sixs_L )
}
test_digit(7, D_4, c(2,6))
test_digit(7, D_4, c(2))
test_digit(7, D_4, c(6))

V_mnist <- function(d) {
  # Conservative Constraint from est_L
  return (d * sevens_L )
}
test_digit(8, D_4, c(2,6))
test_digit(8, D_4, c(2))
test_digit(8, D_4, c(6))

V_mnist <- function(d) {
  # Conservative Constraint from est_L
  return (d * nines_L )
}
test_digit(10, D_4, c(2,6))
test_digit(10, D_4, c(2))
test_digit(10, D_4, c(6))

#### Time Check
system.time(test_digit(4,D_4,c(2,6)))





