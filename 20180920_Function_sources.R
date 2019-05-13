### Some C Codes
library(Rcpp)
library(RcppArmadillo)
library(inline)

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

# Source for the U-statistic estimator of the Euclidean energy statistic
U_stat_source <- '
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
for (int jj = 0; jj < ii; jj++) {
XX += norm(data_1.row(ii)-data_1.row(jj),order_num);
}
}

for (int ii = 0; ii < obs_2; ii++) {
for (int jj = 0; jj < ii; jj++) {
YY += norm(data_2.row(ii)-data_2.row(jj),order_num);
}
}

for (int ii = 0; ii < obs_1; ii++) {
for (int jj = 0; jj < obs_2; jj++) {
XY += norm(data_1.row(ii)-data_2.row(jj),order_num);
}
}

// Create a result storage
double RESULT = 2*XY/(obs_1*obs_2) - 2*XX/obs_1/(obs_1-1) - 2*YY/obs_2/(obs_2-1);

return Rcpp::List::create(
Rcpp::Named("U_stat") = RESULT);
'

# Source for the V-statistic estimator of the radial kernel MMD
V_radialMMD_source <- '
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
for (int jj = 0; jj < obs_1; jj++) {
XX += exp(-pow(norm(data_1.row(ii)-data_1.row(jj),2),2)/2/sigma_sq);
}
}

for (int ii = 0; ii < obs_2; ii++) {
for (int jj = 0; jj < obs_2; jj++) {
YY += exp(-pow(norm(data_2.row(ii)-data_2.row(jj),2),2)/2/sigma_sq);
}
}

for (int ii = 0; ii < obs_1; ii++) {
for (int jj = 0; jj < obs_2; jj++) {
XY += exp(-pow(norm(data_1.row(ii)-data_2.row(jj),2),2)/2/sigma_sq);
}
}

// Create a result storage
double RESULT = XX/pow(obs_1,2) + YY/pow(obs_2,2) - 2*XY/(obs_1*obs_2);

return Rcpp::List::create(
Rcpp::Named("V_radialMMD") = RESULT);
'

# Source for the U-statistic estimator of the radial kernel MMD
U_radialMMD_source <- '
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
double RESULT = 2*XX/obs_1/(obs_1-1) + 2*YY/obs_2/(obs_2-1) - 2*XY/(obs_1*obs_2);

return Rcpp::List::create(
Rcpp::Named("U_radialMMD") = RESULT);
'

### Load functions
# Computes the V-statistic estimator of the Euclidean energy statistic
# Input 2 data sets (matrices) and an order > 0
V_stat_fun <- cxxfunction(signature(data_1_r='numeric',
                                    data_2_r='numeric',
                                    order_r='numeric'),
                          V_stat_source, plugin = 'RcppArmadillo')

# Computes the U-statistic estimator of the Euclidean energy statistic
# Inpute 2 data sets (matrices) and an order > 0
U_stat_fun <- cxxfunction(signature(data_1_r='numeric',
                                    data_2_r='numeric',
                                    order_r='numeric'),
                          U_stat_source, plugin = 'RcppArmadillo')

# Computes the V-statistic estimator of the Euclidean energy statistic
# Input 2 data sets (matrices) and a sigma_sq > 0
V_radialMMD_fun <- cxxfunction(signature(data_1_r='numeric',
                                    data_2_r='numeric',
                                    sigma_sq_r='numeric'),
                               V_radialMMD_source, plugin = 'RcppArmadillo')

# Computes the U-statistic estimator of the Euclidean energy statistic
# Input 2 data sets (matrices) and a sigma_sq > 0
U_radialMMD_fun <- cxxfunction(signature(data_1_r='numeric',
                                         data_2_r='numeric',
                                         sigma_sq_r='numeric'),
                               U_radialMMD_source, plugin = 'RcppArmadillo')