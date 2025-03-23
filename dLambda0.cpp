//dLambda0.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List dLambda0 (const arma::mat& Z, const arma::mat& df_beta_no_eta, const arma::vec& beta, const arma::vec& failtime, 
               const arma::vec& eventime, const arma::vec& wgt, const arma::vec& nevent_wt,
               const arma::uvec& failtime_order, int& nsub, int& ntime){
  
  arma::mat S1(ntime, beta.n_elem);
  arma::vec S0(ntime);
  arma::uvec idxR;
  arma::vec wexpbZ = wgt % exp(Z*beta);
  
  for(int j=0; j<ntime; j++){
    idxR = arma::find(eventime>=failtime(j));
    S0(j) = sum(wexpbZ(idxR));
    S1.row(j) = wexpbZ(idxR).t() * Z.rows(idxR);
  }
  arma::vec dL0t = nevent_wt/S0;
  arma::mat eqdL0t(nsub, ntime), df_dL0t;
  arma::uvec id(nsub, arma::fill::ones);
  id = arma::cumsum(id);
  
  for(int j=0; j<ntime; j++){
    eqdL0t.col(j) = wgt % (id==failtime_order(j)) - dL0t(j) * (wexpbZ % (eventime>=failtime(j)));
  }
  
  df_dL0t = -(df_beta_no_eta * (S1.t()));
  df_dL0t = df_dL0t.each_row() % ((dL0t/S0).t());
  df_dL0t += eqdL0t.each_row() / (S0.t());
  
  return List::create(Named("dL0t")=dL0t, Named("df_dL0t")=df_dL0t);
}