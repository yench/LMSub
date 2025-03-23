//Rho.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
 arma::mat Rho(const arma::vec& H, int ctrl_num,  
               const arma::vec& case_times, const arma::vec& ctrl_times,
               const arma::vec& case_Z3,     const arma::vec& ctrl_Z3){

	 arma::mat Rho(ctrl_num, ctrl_num);

	 for(int i=0; i<(ctrl_num-1); i++){
		  for(int j=i+1; j<(ctrl_num); j++){
		    arma::uvec ind;
        	    ind = arma::find((case_times<ctrl_times(i)) && (case_times<ctrl_times(j)) &&
                                     (case_Z3==ctrl_Z3(i)) && (case_Z3==ctrl_Z3(j)));
  			Rho(i,j) = arma::prod(H(ind))-1;
		}
	}
	return Rho;
}