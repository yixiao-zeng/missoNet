#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List updateBeta(const arma::mat& Theta, arma::mat& B0, 
                const int& n, const arma::mat& xtx, const arma::mat& xty,
                const double& lamB, const double& eta, const double& tolin, const int& maxitrin) {
  arma::mat B = B0, V(size(B0), arma::fill::none), grad_f(size(B0), arma::fill::none), proxgrad_update_B(size(B0), arma::fill::none), Gt(size(B0), arma::fill::none);
  bool stationary = false;
  int itr = 1;
  double t = 1.0, c, lamBt, Q_update, F_update, Q1, Q2, Q3;
  
  while (itr <= maxitrin && !stationary) {
    t = t/eta;
    Q_update = 0.0;
    F_update = Q_update + 1.0;
    c = 1.0 * (itr - 2.0)/(itr + 1.0);
    V = B + c * (B - B0);
    grad_f = 2.0/n * (xtx * V - xty) * Theta;
    Q1 = 1.0/n * arma::trace((V.t() * xtx * V  - 2.0 * V.t() * xty) * Theta);
    
    while ( F_update > Q_update ) {
      t = t * eta;
      lamBt = lamB * t;
      proxgrad_update_B = V - t * grad_f;
      proxgrad_update_B.clean(lamBt);
      proxgrad_update_B.for_each([lamBt](arma::mat::elem_type& val){if (val > 0) {val -= lamBt;} else if (val < 0) {val += lamBt;}});
      
      F_update = 1.0/n * arma::trace((proxgrad_update_B.t() * xtx * proxgrad_update_B - 2.0 * proxgrad_update_B.t() * xty) * Theta);
      Gt = (V - proxgrad_update_B)/t;
      Q2 = t * arma::accu(grad_f % Gt);
      Q3 = t/2.0 * arma::accu(arma::pow(Gt, 2));
      Q_update = Q1 - Q2 + Q3;
    }
    
    if (arma::accu(arma::abs(proxgrad_update_B - B)) <= tolin) {
      stationary = true;
      B = proxgrad_update_B;
    } else {
      itr++;
      B0 = B;
      B = proxgrad_update_B;
    }
  }
  return List::create( _["it.final"] = itr, _["Bhat"] = B );
}


