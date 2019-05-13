#library(MCMCpack)
#library(transport)
library(MixSim)
library(MASS)
library(winference)
library(FNN)
library(Rcpp)
library(RcppArmadillo)
library(inline)

# Weighting function
compute_weight <- function(x, q, eps) {
    return(exp(-abs(x)**q/eps))
}

# Rejection ABC with a tolerance level
ABC_output <- function(distance, quantile_level) {
    threshold <- quantile(distance, probs = quantile_level)
    (distance < threshold)
}

# Generate a Gaussien mixture data set and parameter values from priors
BGM_dataset_simu <- function(true_data, true_covar, nb_iters) {
    # Initialize parameter values
    rand_Pi <- array(NA, c(1,2))
    rand_Mu <- array(NA, c(2,2))
    simu_data <- array(NA, c(nb_iters,nrow(true_data),ncol(true_data)))
    stored_Pi <- array(NA, c(nb_iters,2))
    stored_Mu <- array(NA, c(nb_iters,2,2))
    
    for (k in 1:nb_iters) {
        # rand_Pi <- rdirichlet(1, c(1, 1))
        rand_Pi[1] <- runif(1, min = 0, max = 1)
        rand_Pi[2] <- 1.0 - rand_Pi[1]
        rand_Mu[1,] <- runif(2, min = -1, max = 1)
        rand_Mu[2,] <- runif(2, min = -1, max = 1)
        simu_data[k,,] <- simdataset(nrow(true_data), rand_Pi, rand_Mu, 
                                    true_covar)[[1]]
        stored_Pi[k,] <- rand_Pi
        stored_Mu[k,1,] <- rand_Mu[1,]
        stored_Mu[k,2,] <- rand_Mu[2,]
    }
    
    list(Data = simu_data, Pi = stored_Pi, Mu = stored_Mu)
}

# Generate a MA(2) data set and parameter values from priors
MA2_dataset_simu <- function(true_data, nb_iters) {
    # Initialize parameter values
    rand_Theta <- array(NA, 2)
    simu_data <- array(NA, c(nb_iters,nrow(true_data),ncol(true_data)))
    stored_Theta <- array(NA, c(nb_iters,2))
    
    for (k in 1:nb_iters) {
        rand_Theta[1] <- runif(1, min = -2, max = 2)
        rand_Theta[2] <- runif(1, min = -1, max = 1)
        for (i in 1:nrow(true_data)) {
            simu_data[k,i,1:ncol(true_data)] <- arima.sim(n=ncol(true_data), 
                                        list(ma=rand_Theta), 
                                        rand.gen = function(n) rt(n, df = 5))
        }
        stored_Theta[k,] <- rand_Theta
    }
    
    list(Data = simu_data, Theta = stored_Theta)
}

# Generate a BBM data set and parameter values from priors
BBM_dataset_simu <- function(n_true, nb_iters) {
    # Initialize parameter values
    rand_Theta <- array(NA, 5)
    simu_data <- array(NA, c(nb_iters,n_true,2))
    stored_Theta <- array(NA, c(nb_iters,5))
    
    for (k in 1:nb_iters) {
        rand_Theta <- runif(5, min = 0, max = 5)
        U1 <- rgamma(n_true,rand_Theta[1],1)
        U2 <- rgamma(n_true,rand_Theta[2],1)
        U3 <- rgamma(n_true,rand_Theta[3],1)
        U4 <- rgamma(n_true,rand_Theta[4],1)
        U5 <- rgamma(n_true,rand_Theta[5],1)
        V1 <- (U1+U3)/(U5+U4)
        V2 <- (U2+U4)/(U5+U3)
        Z1 <- V1/(V1+1)
        Z2 <- V2/(V2+1)
        simu_data[k,,] <- cbind(Z1,Z2)
        # 
        # for (i in 1:nrow(true_data)) {
        #     simu_data[k,i,1] <- rbeta(1, rand_Theta[1]+rand_Theta[3], 
        #                               rand_Theta[4]+rand_Theta[5])
        #     simu_data[k,i,2] <- rbeta(1, rand_Theta[2]+rand_Theta[4], 
        #                               rand_Theta[3]+rand_Theta[5])            
        # }
        stored_Theta[k,] <- rand_Theta
    }
    
    list(Data = simu_data, Theta = stored_Theta)
}

# Generate a 5-dimensional g-and-k data set and parameter values from priors
gandk_quantile_samples <- function(z, theta) {
    
    A <- theta[1]
    B <- theta[2]
    c <- 0.8
    g <- theta[3]
    k <- theta[4]
    
    A + B * (1 + c * (1 - exp(- g * z)) / (1 + exp(- g * z))) * 
        (1 + z ** 2) ** k * z
}

GK_dataset_simu <- function(true_data, nb_iters) {
    # Initialize parameter values Theta = (A, B, g, k, rho)
    rand_Theta <- array(NA, 5)
    simu_data <- array(NA, c(nb_iters,nrow(true_data),5))
    stored_Theta <- array(NA, c(nb_iters,5))
    
    for (k in 1:nb_iters) {
        rand_Theta[1:4] <- runif(4, min = 0, max = 4)
        rand_Theta[5] <- runif(1, min = -0.5, max = 0.5)
        Sigma <- matrix(c(1,rand_Theta[5],0,0,0,
                          rand_Theta[5],1,rand_Theta[5],0,0,
                          0,rand_Theta[5],1,rand_Theta[5],0,
                          0,0,rand_Theta[5],1,rand_Theta[5],
                          0,0,0,rand_Theta[5],1),5,5)
        for (i in 1:nrow(true_data)) {
            rand_Z <- mvrnorm(1, mu = rep(0, 5), Sigma = Sigma)
            simu_data[k,i,] <- gandk_quantile_samples(rand_Z, rand_Theta)            
        }
        stored_Theta[k,] <- rand_Theta
    }
    
    list(Data = simu_data, Theta = stored_Theta)
}

# Compute weights via the V-statistic estimator
ES_distance <- function(true_data, simu_data) {
    # Initialize parameter values
    nb_iters <- nrow(simu_data)
    weights <- array(NA, nb_iters)
    distance <- array(NA, c(nb_iters,1))
    
    for (k in 1:nb_iters) {
        distance[k,] <- V_stat_fun(true_data, simu_data[k,,], 2)$V_stat
    }
    
    distance
}

# Compute weights via the KL divergence
KL_distance <- function(true_data, simu_data) {
    # Initialize parameter values
    nb_iters <- nrow(simu_data)
    weights <- array(NA, nb_iters)
    distance <- array(NA, c(nb_iters,1))
    
    for (k in 1:nb_iters) {
        distance[k,] <- KLx.divergence(true_data, simu_data[k,,], 1)
    }
    
    distance
}


# Compute weights via the Wasserstein distance (too slow)
WA_distance <- function(true_data, simu_data) {
    # Initialize parameter values
    nb_iters <- nrow(simu_data)
    weights <- array(NA, nb_iters)
    distance <- array(NA, c(nb_iters,1))
    
    for (k in 1:nb_iters) {
        distance[k,] <- swap_distance(t(pp(true_data)$coordinates), 
                                       t(pp(simu_data[k,,])$coordinates), 
                                       p=2, tolerance=1e-5)$distance
    }
    
    distance
}

# Compute weights via the MMD estimator
MMD_distance <- function(true_data, simu_data) {
    # Initialize parameter values
    nb_iters <- nrow(simu_data)
    weights <- array(NA, nb_iters)
    diff <- array(NA, c(nrow(true_data)*(nrow(true_data)-1)/2,1))
    distance <- array(NA, c(nb_iters,1))
    
    idx <- 0
    for (i in 1:nrow(true_data)) {
        for (j in i:nrow(true_data)) {
            idx <- idx + 1
            diff[idx] <- norm(as.matrix(true_data[i,]-true_data[j,]), 
                              type = "F")
        }
    }
    
    sigma_sq <- median(diff)**2
    
    for (k in 1:nb_iters) {
        distance[k,] <- MMD_dist_fun(true_data, simu_data[k,,], sigma_sq)$MMD
    }
    
    distance
}

# 2D Kernel Density Estimation (weighted) by Nicholas Hamilton
kde2d.weighted <- function (x, y, w, h, n = 100, lims = c(range(x), range(y))) {
    nx <- length(x)
    if (length(y) != nx) 
        stop("data vectors must be the same length")
    if (length(w) != nx & length(w) != 1)
        stop("weight vectors must be 1 or length of data")
    gx <- seq(lims[1], lims[2], length = n) # gridpoints x
    gy <- seq(lims[3], lims[4], length = n) # gridpoints y
    if (missing(h)) 
        h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
    if (missing(w)) 
        w <- numeric(nx)+1;
    h <- h/4
    ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
    ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
    # z is the density
    z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% 
        t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2])
    
    return(list(x = gx, y = gy, z = z))
}

### Sources
# Source for the V-statistic estimator of the Euclidean energy statistic
V_stat_source <- '
// Set namespaces
using namespace arma;
using namespace Rcpp;

// Load inpute variables
mat data_1 = as<mat>(data_1_r);
mat data_2 = as<mat>(data_2_r);
double order_num = as<double>(order_r);

// Get object dimensions
int obs_1 = data_1.n_rows;
int obs_2 = data_2.n_rows;

// Initialize three storage objects
double XX = 0;
double YY = 0;
double XY = 0;

// Compute required terms
for (int ii = 0; ii < obs_1; ii++) {
for (int jj = 0; jj < obs_1; jj++) {
XX += norm(data_1.row(ii)-data_1.row(jj),order_num);
}
}

for (int ii = 0; ii < obs_2; ii++) {
for (int jj = 0; jj < obs_2; jj++) {
YY += norm(data_2.row(ii)-data_2.row(jj),order_num);
}
}

for (int ii = 0; ii < obs_1; ii++) {
for (int jj = 0; jj < obs_2; jj++) {
XY += norm(data_1.row(ii)-data_2.row(jj),order_num);
}
}

// Create a result storage
double RESULT = 2*XY/(obs_1*obs_2) - XX/pow(obs_1,2) - YY/pow(obs_2,2);

return Rcpp::List::create(
Rcpp::Named("V_stat") = RESULT);
'

# Source for the U-statistic estimator of the radial kernel MMD
MMD_source <- '
// Set namespaces
using namespace arma;
using namespace Rcpp;

// Load inpute variables
mat data_1 = as<mat>(data_1_r);
mat data_2 = as<mat>(data_2_r);
double sigma_sq = as<double>(sigma_sq_r);

// Get object dimensions
int obs_1 = data_1.n_rows;
int obs_2 = data_2.n_rows;

// Initialize three storage objects
double XX = 0;
double YY = 0;
double XY = 0;

// Compute required terms
for (int ii = 0; ii < obs_1; ii++) {
for (int jj = 0; jj <= ii; jj++) {
XX += exp(-pow(norm(data_1.row(ii)-data_1.row(jj),2),2)/2/sigma_sq);
}
}

for (int ii = 0; ii < obs_2; ii++) {
for (int jj = 0; jj <= ii; jj++) {
YY += exp(-pow(norm(data_2.row(ii)-data_2.row(jj),2),2)/2/sigma_sq);
}
}

for (int ii = 0; ii < obs_1; ii++) {
for (int jj = 0; jj < obs_2; jj++) {
XY += exp(-pow(norm(data_1.row(ii)-data_2.row(jj),2),2)/2/sigma_sq);
}
}

// Create a result storage
double RESULT = 0;
RESULT = sqrt(2*XX/obs_1/(obs_1-1) + 2*YY/obs_2/(obs_2-1) - 2*XY/(obs_1*obs_2));

return Rcpp::List::create(Rcpp::Named("MMD") = RESULT);
'

# Computes the V-statistic estimator of the Euclidean energy statistic
# Input 2 data sets (matrices) and an order > 0
V_stat_fun <- cxxfunction(signature(data_1_r='numeric',
                                    data_2_r='numeric',
                                    order_r='numeric'),
                          V_stat_source, plugin = 'RcppArmadillo')

# Computes the U-statistic estimator of the Euclidean energy statistic
# Input 2 data sets (matrices) and a sigma_sq > 0
MMD_dist_fun <- cxxfunction(signature(data_1_r='numeric', 
                                 data_2_r='numeric', 
                                 sigma_sq_r='numeric'), 
                       MMD_source, plugin = 'RcppArmadillo')