// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]   

//'@importFrom Rcpp sourceCpp
//'@useDynLib MRT


#include "RcppArmadillo.h"

//BOOST
#include "boost/math/distributions/students_t.hpp"
#include "boost/program_options.hpp"
#include "boost/math/tools/roots.hpp"
#include "boost/math/special_functions/digamma.hpp"
#include <algorithm>

using namespace Rcpp;


Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


double phi(double d){
  double c = 1.345;
  double temp = 1.0;
  
  if(std::abs(d) >= c) {
    temp= c/std::abs(d);
  }
  
  return temp; 
}


double phi2(double d){
  double c = 4.685;
  double temp = 0.0;
  
  if(std::abs(d) < c) {
    temp = std::pow(1-std::pow(d/c,2.0),2.0);
  }
  
  return temp; 
}


double phi3(double d){
  double c = 3;
  double temp = 0.0;
  
  if(std::abs(d) < c) {
    temp = std::pow(1-std::pow(d/c,2.0),1.0);
  }
  
  return temp; 
}

// [[Rcpp::export]]
double mahalanobis_c(double y, double mu, double sigma)
{
  double delta = (y-mu)*(1./sigma)*(y-mu);
  return delta;
}

// [[Rcpp::export]]
double mahalanobis_HD(arma::vec y, arma::vec mu, arma::mat sigma)
{
  arma::vec d = y-mu;
  double delta = std::sqrt(as_scalar(d.t()*pinv(sigma)*d));
  return delta;
}


// [[Rcpp::export]]
double norm_c(double y, double mu, double sigma)
{
  double f = (1/(std::sqrt(arma::datum::pi*2.0*sigma))) * std::exp(-0.5*mahalanobis_c(y,mu,std::sqrt(sigma)));
 // return(f <= std::pow(1, -100)) ? 0.0 : f;
  return f;
}


arma::mat start_kmeans(arma::vec dat, int g){
  arma::mat params(g,3, arma::fill::zeros);
  arma::vec weights(1, arma::fill::ones);
  arma::uvec alloc(dat.n_elem, arma::fill::zeros);
  arma::uvec cluster_sizes(g, arma::fill::zeros);
  
  Rcpp::Environment base("package:stats"); 
  Rcpp::Function kmeans = base["kmeans"];  
  
  Rcpp::List res = kmeans(Rcpp::_["x"] = dat, Rcpp::_["centers"]  = g);
  
  arma::vec tmp_mu = res["centers"];
  params.col(1) = tmp_mu;
  
  
  for(int i=0; i<dat.n_elem;i++){
    arma::vec temp = abs(params.col(1)- dat(i)); 
    alloc.at(i) = arma::index_min(temp);
  }
  
  //cluster size and variance
  for(int i =0; i<g;i++){
    arma::uvec tmp1 = find(alloc == i);
    cluster_sizes.at(i) = tmp1.n_elem;
    params(i,2) = sum(pow(dat(tmp1)-params(i,1),2.0))/dat.n_elem;
   }
  
  params.col(0) = arma::conv_to< arma::vec >::from(cluster_sizes)/dat.n_elem;
  return(params);
}

//'@export
// [[Rcpp::export]]
Rcpp::List MMEst(const arma::vec & y, int g = 1, double tol = 0.00001, int max_iter = 100) {
  
  int n = y.n_elem;
  
  arma::mat d  = arma::zeros(g,n);
  arma::mat u  =  arma::zeros(g,n);
  arma::mat w  =  arma::ones(g,n);
  arma::mat tau  =  arma::zeros(g,n);

  //allocate double
  arma::vec mu = arma::zeros(g);
  arma::vec sigma = arma::zeros(g);
  arma::vec pi = arma::zeros(g);
  
  double diff = 1.0;
  
  arma::mat init_mat = arma::zeros(g,3);
  
  init_mat = start_kmeans(y,g);
  
  
  pi = init_mat.col(0);
  mu = init_mat.col(1);
  sigma = init_mat.col(2);
  
  //allocated int
  int k =0;
  
  //iterative fit for M-est
  while(diff>tol & k <max_iter ){
    
    for(int i=0;i<g;i++){
      
      double mu_i =mu(i);
      double sigma_i = sigma(i);
      
      arma::rowvec tmp = arma::zeros(1,n);
      std::transform(y.begin(), y.end(), tmp.begin(), [mu_i, sigma_i](double dat) { return mahalanobis_c(dat, mu_i, std::sqrt(sigma_i)); });
      d.row(i) = tmp;

      std::transform(y.begin(), y.end(), tmp.begin(), [mu_i, sigma_i](double dat) {return norm_c(dat, mu_i, sigma_i);} );
      tau.row(i) = pi(i)*tmp;
      
    }
    
    tau = tau.each_row() / sum(tau,0);
    
    std::transform(d.begin(), d.end(), w.begin(), [](double val) { return phi(val);} );
    u = w;

    for(int i=0;i<g;i++){
      mu(i) = sum(tau.row(i)%u.row(i)%y.t()) / sum(tau.row(i)%u.row(i)); 
      double mu_i = mu(i);
      arma::rowvec tmp = y.t() - mu_i;
      sigma(i) = sum(tau.row(i)%u.row(i)%(tmp)%(tmp))/  sum(tau.row(i)); //sum(tau.row(i));
      pi(i)= sum(tau.row(i))/accu(tau);
    }  
    ++k;
   }
  
  double out_thresh_const = 7;
  //double outly = -2*std::log(out_thresh_const*std::sqrt(2*3.141593*sigma)); 
  
  Rcpp::List temp_list = Rcpp::List::create(
    Rcpp::Named("mu")= export_vec(mu),
    Rcpp::Named("sigma") = export_vec(sigma),
    Rcpp::Named("pi") = export_vec(pi),
    Rcpp::Named("tau") = tau,
    Rcpp::Named("g") = g,
    Rcpp::Named("n_iter") =k
    
  );
  
  return temp_list;
  
  }


