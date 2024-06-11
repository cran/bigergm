#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 0
#endif
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

inline void set_seed(int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

// Function that simulates a between-block network.
// The first element of `coef_between` must be the edges parameter.
// [[Rcpp::export]]
arma::sp_mat simulate_between_network
(int numOfVertices,
 const Rcpp::List& list_feature_adjmat,
 const arma::vec& coef_between,
 const arma::vec& block_membership,
 bool directed, 
 int & seed
)
{
  // Number of covariates
  int numOfCovariates = list_feature_adjmat.length();
  // Initialize a sparse adjacency matrix for the between-block network
  arma::sp_mat between_adjmat = arma::sp_mat(numOfVertices, numOfVertices);
  // Prepare a sparse adjacency cube
  arma::field<arma::sp_mat> feature_cube(numOfCovariates);
  for (int p = 0; p < numOfCovariates; p++) {
    feature_cube(p) = Rcpp::as<arma::sp_mat>(list_feature_adjmat[p]);
  }
  
  // Necessary for R random number generator
  GetRNGstate();
  
#pragma omp parallel
{
  // Simulate between-block links
#pragma omp for
  for (int j = 0; j < numOfVertices; j++) {
    std::mt19937 generator(seed+j);
    std::uniform_real_distribution<double>  distr(0.0,1.0);
    
    for (int i = 0; i < numOfVertices; i++) {
      // Skip as many unnecessary calculations as possible in this nested loop, which makes the computation faster.
      if (block_membership[i] != block_membership[j] && ((directed && i != j) || (!directed && i < j))) {
        double x = distr(generator);
        double u = coef_between[0];
        for (int p = 0; p < numOfCovariates; p++) {
          double elem = feature_cube(p)(i, j);
          double elem_coef = coef_between[p+1];
          u += elem_coef * elem;
        }
        //std::printf("x: %f, Thread: %d, Loop: (%d, %d), u: %f\n", x, omp_get_thread_num(), i, j, u);
        if (u > log(x/(1-x))) {
          between_adjmat(i, j) = 1;
        }
      }
    }
  }
}
// This must be called after GetRNGstate before returning to R.
PutRNGstate();
return between_adjmat;
}

// [[Rcpp::export]]
arma::sp_mat simulate_between_network_covariates
(const int numOfVertices,
 const Rcpp::List & coef_between,
 const Rcpp::List & list_feature_adjmat,
 const arma::vec& block_membership,
 bool directed,
 int & seed
)
{ // Number of covariates
  int numOfCovariates = list_feature_adjmat.length();
  int numOfOptions = coef_between.length();
  // Initialize a sparse adjacency matrix for the between-block network
  arma::sp_mat between_adjmat = arma::sp_mat(numOfVertices, numOfVertices);
  // Prepare a sparse adjacency cube
  arma::field<arma::sp_mat> feature_cube(numOfCovariates);
  // Prepare a sparse probability cube
  arma::field<arma::sp_mat> probability_cube(numOfOptions);
  for (int p = 0; p < numOfOptions; p++) {
    probability_cube(p) = Rcpp::as<arma::sp_mat>(coef_between.at(p));
  }
  for (int p = 0; p < numOfCovariates; p++) {
    feature_cube(p) = Rcpp::as<arma::sp_mat>(list_feature_adjmat.at(p));
  }
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double>  distr(0.0,1.0);
  
  
  // Necessary for R random number generator
  GetRNGstate();
  for (int j = 0; j < numOfVertices; j++) {
    for (int i = 0; i < numOfVertices; i++) {
      // Skip as many unnecessary calculations as possible in this nested loop, which makes the computation faster.
      if (block_membership[i] != block_membership[j] && ((directed && i != j) || (!directed && i < j))) {
        int index = 0;
        double x = unif_rand();
        for (int p = 0; p < numOfCovariates; p++) {
          if(feature_cube(p)(i, j)){
            index += pow(2,p);
          }
        }
        double u = probability_cube(index).at(block_membership[i]-1,block_membership[j]-1);
        if (u > x) {
          // Rcpp::Rcout << "Between " + std::to_string(block_membership[i]) + " and " +std::to_string(block_membership[j]) + " with p = " + std::to_string(u)<< std::endl;
          // Rcpp::Rcout << "i = " + std::to_string(i+1)<< std::endl;
          // Rcpp::Rcout << "j = " + std::to_string(j+1)<< std::endl;
          // Rcpp::Rcout << "Cov 1 = " + std::to_string(feature_cube(0)(i, j))<< std::endl;
          // Rcpp::Rcout << "Cov 2 = " + std::to_string(feature_cube(1)(i, j))<< std::endl;
          between_adjmat(i, j) = 1;
        } 
      }
    }
  }
// This must be called after GetRNGstate before returning to R.
PutRNGstate();
return between_adjmat;
}

// [[Rcpp::export]]
arma::sp_mat simulate_between_network_no_covariates
(const int numOfVertices,
 const arma::sp_mat & coef_between,
 const arma::vec& block_membership,
 bool directed, 
 int & seed)
{ 
  // Initialize a sparse adjacency matrix for the between-block network
  arma::sp_mat between_adjmat = arma::sp_mat(numOfVertices, numOfVertices);
  
  // Necessary for R random number generator
  GetRNGstate();
#pragma omp parallel
{
  // Simulate between-block links
#pragma omp for 
  for (int j = 0; j < numOfVertices; j++) {
    std::uniform_real_distribution<double>  distr(0.0,1.0);
    std::mt19937 generator(seed+j);
    
    for (int i = 0; i < numOfVertices; i++) {
      // Skip as many unnecessary calculations as possible in this nested loop, which makes the computation faster.
      if (block_membership[i] != block_membership[j] && ((directed && i != j) || (!directed && i < j))) {
        double x = distr(generator);
        double u = coef_between.at(block_membership[i]-1,block_membership[j]-1);
        if (u > x) {
          // Rcpp::Rcout << "Between " + std::to_string(block_membership[i]) + " and " +std::to_string(block_membership[j]) + "with p = " + std::to_string(u)<< std::endl;
          between_adjmat(i, j) = 1;
        } 
      }
    }
  }
}
// This must be called after GetRNGstate before returning to R.
PutRNGstate();
return between_adjmat;
}

