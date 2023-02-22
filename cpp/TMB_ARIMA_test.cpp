#include <TMB.hpp>

using namespace Eigen;
using namespace density;

template<class Type>
vector<Type> make_phi(vector<Type> pac) {
  // https://inla.r-inla-download.org/r-inla.org/doc/latent/ar.pdf
  if (pac.size() == 1) {
    return pac;
  }
  vector<Type> phi(pac);
  vector<Type> work(pac);
  Type a;
  for (int i = 0; i < pac.size()-1; i++) {
    a = phi[i + 1];
    vector<Type> phi_rev = phi.segment(0, i+1).reverse();
    work.segment(0, i+1) -= a * phi_rev;
    phi.segment(0, i+1) = work.segment(0, i+1);
  }
  return phi;
}

template<class Type>
Type invlink_phi(Type logit_phi) {
  return 1/(1 + exp(-logit_phi));
}
template<class Type>
vector<Type> invlink_phi(vector<Type> logit_phi) {
  return 1/(1 + exp(-logit_phi));
}

template<class Type>
Type AR(vector<Type> x, Type log_sigma, vector<Type> logit_phi,
        Type phi_sd = Type(1), Type sigma_sd = Type(1),
        Type constr_sd = 0.1, bool constr = true
) {
  Type sigma = exp(log_sigma);
  vector<Type> phi = make_phi(invlink_phi(logit_phi));
  Type res = -ARk(phi)(x/sigma);
  
  res += dnorm(sigma, Type(0), sigma_sd, true);
  res += sum(dnorm(logit_phi, Type(0), phi_sd, true));
  res += log_sigma;
  
  if (constr) {
    res += dnorm(sum(x), Type(0), constr_sd * x.size(), true);
  }
  
  return res;
}

template<class Type>
Type ARIMA_1d0(vector<Type> x, SparseMatrix<Type> D, Type log_sigma,
               vector<Type> logit_phi, int constr_type = 1, Type constr_sd = 0.001,
               Type phi_sd = Type(1),
               Type sigma_sd = Type(1)
) {
  vector<Type> x_diff = (D * x);
  Type res = AR(x_diff, log_sigma, logit_phi, phi_sd, sigma_sd, constr_sd, false);
  if (constr_type == 1) {
    res += dnorm(x[0]/exp(log_sigma), Type(0), constr_sd, true);
  } else if (constr_type == 2) {
    res += dnorm(sum(x)/exp(log_sigma), Type(0), x.size() * constr_sd, true);
  }
  
  return res;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_VECTOR(y);
  DATA_SPARSE_MATRIX(D);
  DATA_INTEGER(constr_type);
  DATA_SCALAR(constr_sd);
  DATA_MATRIX(IS_map);

  PARAMETER_VECTOR(mu);
  PARAMETER(log_sigma_mu);
  PARAMETER_VECTOR(logit_phi);
  PARAMETER(log_sigma_y);
  
  Type sigma_y = exp(log_sigma_y);
  
  Type nll = 0.0;
  nll -= ARIMA_1d0(mu, D, log_sigma_mu, logit_phi, constr_type, constr_sd);
  
  vector<Type> mu_IS = IS_map * mu;
  nll -= sum(dnorm(y, mu_IS, sigma_y, true));
  nll -= dnorm(sigma_y, Type(0), Type(2), true);
  nll -= log_sigma_y;
  
  vector<Type> phi = make_phi(invlink_phi(logit_phi));
  ADREPORT(phi);
  
  REPORT(mu);
  
  return nll; 
}
