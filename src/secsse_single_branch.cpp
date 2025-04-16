#include <cstdlib>    // std::getenv, std::atoi
#include <vector> 
#include <chrono>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]

template <typename ODE>
Rcpp::List call_ll_single_branch_broken(std::unique_ptr<ODE> od,
                                        const Rcpp::NumericVector& states,
                                        const Rcpp::NumericVector& forTime,
                                        const std::string& method,
                                        double atol,
                                        double rtol,
                                        bool see_states) {
  
  auto t0_orig = std::min(forTime[0], forTime[1]);
  auto t1_orig = std::max(forTime[0], forTime[1]);
  
  auto T0 = std::chrono::high_resolution_clock::now();
  
  auto states_out = std::vector<double>(states.begin(), states.end());
  auto workhorse = secsse::Integrator<ODE>(std::move(od), method, atol, rtol, false);
  auto d = workhorse.size();
  
  double t0 = t0_orig;
  const double max_dt = 0.5;
  
  double t1 = max_dt;
  double loglik = 0.0;

 // std::cerr << "Hello, we are breaking up the branch!\n";
  while(true) {
   // std::cerr << t0 << " " << t1 << " " << loglik << "\n";
    workhorse(states_out, t0, t1);
    
    loglik += secsse::normalize_loglik(std::begin(states_out) + d, 
                                           std::end(states_out));
    
    if (t1 == t1_orig) break;
    
    auto dt = std::min(max_dt, t1_orig - t1);
    t0 = t1;
    t1 += dt;
  }

  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);

  const auto merge_branch = std::vector<double>(std::begin(states_out) + d, 
                                                std::end(states_out));
  
  return Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("merge_branch") = merge_branch,
                            Rcpp::Named("states") = states_out,
                            Rcpp::Named("duration") = DT.count());
}

template <typename ODE>
Rcpp::List call_ll_single_branch_obs(std::unique_ptr<ODE> od,
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
  

  
  auto workhorse = secsse::Integrator<ODE>(std::move(od), method, atol, rtol, true);
  
//  auto observe = odeintcpp::normalizing_observer< std::vector<double> >(100);
  
  double loglik_obs = 0.0;
  workhorse(states_out, t0, t1, loglik_obs);
  
  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);
  
  auto d = workhorse.size();
  
  const auto loglik = loglik_obs + secsse::normalize_loglik(std::begin(states_out) + d, 
                                                                  std::end(states_out));
  
  const auto merge_branch = std::vector<double>(std::begin(states_out) + d, 
                                                std::end(states_out));
  
  return Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("merge_branch") = merge_branch,
                            Rcpp::Named("states") = states_out,
                            Rcpp::Named("duration") = DT.count());
}




template <typename ODE>
Rcpp::List calc_ll_single_branch(std::unique_ptr<ODE> od,
                                 const Rcpp::NumericVector& states,
                                 const Rcpp::NumericVector& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol,
                                 bool see_states,
                                 bool use_normalization) {

  auto t0 = std::min(forTime[0], forTime[1]);
  auto t1 = std::max(forTime[0], forTime[1]);
  
/*  if (t1 > 0 && break_up == 1) {
    return call_ll_single_branch_broken(std::move(od), states, forTime, method, atol, rtol, see_states);
  } 
  
  if (t1 > 0 && break_up == 2) {
    return call_ll_single_branch_obs(std::move(od), states, forTime, method, atol, rtol, see_states);
  } */
  
  
  auto T0 = std::chrono::high_resolution_clock::now();
  
  auto states_out = std::vector<double>(states.begin(), states.end());
  
  auto workhorse = secsse::Integrator<ODE>(std::move(od), method, atol, rtol, use_normalization);
  
  workhorse(states_out, t0, t1);
  
  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);
  
  auto d = workhorse.size();

  const auto loglik = secsse::normalize_loglik(std::begin(states_out) + d, 
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
    return calc_ll_single_branch(std::make_unique<ode_standard<OdeVariant::normal_tree>>(ll, mus, Q), states, forTime, method, atol, rtol, see_states, use_normalization);
  } else if (rhs == "ode_cla") {
    auto ll = Rcpp::as<Rcpp::List>(lambdas);
  
    return calc_ll_single_branch(std::make_unique<ode_cla<OdeVariant::normal_tree>>(ll, mus, Q), states, forTime, method, atol, rtol, see_states, use_normalization);

  } else {
    throw std::runtime_error("calc_ll_cpp: unknown rhs");
  }
}
