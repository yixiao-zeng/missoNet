#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List updateBeta(arma::mat Theta, arma::mat B0,
                double lamB, double eta, double tolin, int maxitrin, const List& info)
{ int n = info["n"];
  arma::mat xtx = info["xtx"], xty = info["til.xty"];
  arma::mat xtyTheta = xty*Theta;
  arma::mat B = B0, V, grad_f, subgrad_update_B, B_sign, proxgrad_update_B, Gt;
  bool stationary = false;
  int itr = 1;
  double t = 1.0, c, Q_update, F_update, Q1, Q2, Q3;
  
  while (itr <= maxitrin && !stationary) {
    t = t/eta;
    Q_update = 0.0;
    F_update = Q_update + 1.0;
    c = 1.0*(itr-2.0)/(itr+1.0);
    V = B + c*(B - B0);
    grad_f = 2.0/n * (xtx*V*Theta-xtyTheta);
    Q1 = 1.0/n * arma::trace(V.t()*xtx*V*Theta - 2.0*V.t()*xtyTheta);
    
    while ( F_update > Q_update ) {
      t = t * eta;
      subgrad_update_B = V - t*grad_f;
      subgrad_update_B.clean(lamB*t);
      B_sign = arma::sign(subgrad_update_B);
      proxgrad_update_B = arma::abs(subgrad_update_B) - lamB*t;
      proxgrad_update_B = proxgrad_update_B % B_sign;
      
      Gt = (V - proxgrad_update_B)/t;
      
      F_update = 1.0/n * arma::trace(proxgrad_update_B.t()*xtx*proxgrad_update_B*Theta - 2.0*proxgrad_update_B.t()*xtyTheta);
      Q2 = t * arma::accu(grad_f % Gt);
      Q3 = t/2.0 * arma::accu(arma::pow(Gt,2));
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


