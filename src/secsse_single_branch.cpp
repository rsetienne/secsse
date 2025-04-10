#include <cstdlib>    // std::getenv, std::atoi
#include <vector> 
#include <chrono>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]


template <typename ODE>
Rcpp::List calc_ll_single_branch(std::unique_ptr<ODE> od,
                                 const Rcpp::NumericVector& states,
                                 const Rcpp::NumericVector& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol,
                                 bool see_states,
                                 bool use_log_transform = false) {
  auto T0 = std::chrono::high_resolution_clock::now();
  
  auto states_out = std::vector<double>(states.begin(), states.end());
  auto states_start = states_out;
  auto workhorse = secsse::Integrator<ODE>(std::move(od), method, atol, rtol);

  auto t0 = std::min(forTime[0], forTime[1]);
  auto t1 = std::max(forTime[0], forTime[1]);

  // std::cerr << t0 << " " << t1 << "\n";
   
  workhorse(states_out, t0, t1);
  
  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);
  
  auto d = workhorse.size();
  
  auto temp = states_out;
  
  // to normalize the loglik, we have to exponentiate the states
  if (use_log_transform) {
    for (size_t i = 0; i < d; ++i) {
      states_out[i + d] = std::exp(states_out[i + d]);
    }
  }
  
 /* for (auto i : states_out) {
    std::cerr << i << " ";
 } std::cerr << "\n";*/
  
  
  const auto loglik = secsse::normalize_loglik(std::begin(states_out) + d, 
                                               std::end(states_out));
  
  const auto merge_branch = std::vector<double>(std::begin(states_out) + d, 
                                                std::end(states_out));
  
  return Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("merge_branch") = merge_branch,
                            Rcpp::Named("states") = states_out,
                            Rcpp::Named("duration") = DT.count());
}

Rcpp::List log_transform_analysis(const Rcpp::List& ll,
                                  const Rcpp::NumericVector& mus,
                                  const Rcpp::NumericMatrix& Q,
                                  const Rcpp::NumericVector& states,
                                  const Rcpp::NumericVector& forTime,
                                  const std::string& method,
                                  double atol,
                                  double rtol,
                                  bool see_states) {
  using namespace secsse;
  // first we integrate a very, very short timestep to get rid of zeros
  const static double mini_step = 1e-6;
  const static Rcpp::NumericVector shortForTime =  {0, mini_step};
  const Rcpp::NumericVector newForTime = {forTime[0], forTime[1] - mini_step};

  auto res = calc_ll_single_branch(std::make_unique<ode_cla<OdeVariant::normal_tree>>(ll, mus, Q),
                                   states, shortForTime, method, atol, rtol, see_states, false);
  
  Rcpp::NumericVector new_states = res["states"];

  std::vector<double> new_states_c(new_states.begin(), new_states.end());

  size_t d = new_states.size() / 2;
  for (size_t i = 0; i < d; ++i) {
    new_states[i + d] = std::log(new_states[i + d]);
  }

  std::vector<double> new_states2_c(new_states.begin(), new_states.end());
  
  // now we have an Rcpp::List with what we need
  return calc_ll_single_branch(std::make_unique<ode_cla_log<OdeVariant::normal_tree>>(ll, mus, Q),
                               new_states, newForTime, method, atol, rtol, see_states,
                               true);
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
                                     bool use_log_transform)
{
  using namespace secsse;
  
  if (rhs == "ode_standard") {
    auto ll = Rcpp::as<Rcpp::NumericVector>(lambdas);
    return calc_ll_single_branch(std::make_unique<ode_standard<OdeVariant::normal_tree>>(ll, mus, Q), states, forTime, method, atol, rtol, see_states);
  } else if (rhs == "ode_cla") {
    auto ll = Rcpp::as<Rcpp::List>(lambdas);
    
    if (!use_log_transform) {
      return calc_ll_single_branch(std::make_unique<ode_cla<OdeVariant::normal_tree>>(ll, mus, Q), states, forTime, method, atol, rtol, see_states);
    } else {
      return log_transform_analysis(ll, mus, Q, states, forTime, method, atol, rtol, see_states); 
    }
  } else {
    throw std::runtime_error("calc_ll_cpp: unknown rhs");
  }
}


