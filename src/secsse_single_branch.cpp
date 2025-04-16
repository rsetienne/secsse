#include <cstdlib>    // std::getenv, std::atoi
#include <vector> 
#include <chrono>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]



template <typename NORMALIZER,
          typename ODE>
Rcpp::List calc_ll_single_branch(std::unique_ptr<ODE> od,
                                 const Rcpp::NumericVector& states,
                                 const Rcpp::NumericVector& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol,
                                 bool see_states) {

  auto t0 = std::min(forTime[0], forTime[1]);
  auto t1 = std::max(forTime[0], forTime[1]);
  
  auto T0 = std::chrono::high_resolution_clock::now();
  
  auto states_out = std::vector<double>(states.begin(), states.end());
  

  auto workhorse = secsse::Integrator<ODE, NORMALIZER>(std::move(od), method, atol, rtol);

  
  NORMALIZER norm;
  workhorse(states_out, t0, t1, norm);
    

  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);
  
  auto d = workhorse.size();

  const auto loglik = norm.loglik + secsse::normalize_loglik(std::begin(states_out) + d, 
                                                             std::end(states_out));
  
  const auto merge_branch = std::vector<double>(std::begin(states_out) + d, 
                                                std::end(states_out));
  
  return Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("merge_branch") = merge_branch,
                            Rcpp::Named("states") = states_out,
                            Rcpp::Named("duration") = DT.count());
}


// [[Rcpp::export]]
Rcpp::List calc_ll_single_branch_cpp(const std::string& rhs,
                                     const Rcpp::NumericVector& states,
                                     const Rcpp::NumericVector& forTime,
                                     const Rcpp::RObject& lambdas,
                                     const Rcpp::NumericVector& mus,
                                     const Rcpp::NumericMatrix& Q,
                                     const std::string& method,
                                     double atol,
                                     double rtol,
                                     bool see_states,
                                     bool use_normalization) {
  using namespace secsse;
  
  if (rhs == "ode_standard") {
    auto ll = Rcpp::as<Rcpp::NumericVector>(lambdas);
    
    return use_normalization ? 
    calc_ll_single_branch<odeintcpp::normalize>(std::make_unique<ode_standard<OdeVariant::normal_tree>>(ll, mus, Q), states, forTime, method, atol, rtol, see_states) : 
      calc_ll_single_branch<odeintcpp::no_normalization>(std::make_unique<ode_standard<OdeVariant::normal_tree>>(ll, mus, Q), states, forTime, method, atol, rtol, see_states);
    
  } else if (rhs == "ode_cla") {
    auto ll = Rcpp::as<Rcpp::List>(lambdas);
  

    return use_normalization ? 
      calc_ll_single_branch<odeintcpp::normalize>(std::make_unique<ode_cla<OdeVariant::normal_tree>>(ll, mus, Q), states, forTime, method, atol, rtol, see_states) :
      calc_ll_single_branch<odeintcpp::no_normalization>(std::make_unique<ode_cla<OdeVariant::normal_tree>>(ll, mus, Q), states, forTime, method, atol, rtol, see_states);
  
  } else {
    throw std::runtime_error("calc_ll_cpp: unknown rhs");
  }
}
