// Copyright 2023 Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#include <cstdlib>    // std::getenv, std::atoi
#include <vector> 
#include <chrono>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]


namespace secsse {

  // probably the cleanest way to retrieve RcppParallel's concurrency setting
  // set by RcppParallel::setThreadOptions(numThreads)
  size_t get_rcpp_num_threads() {
    auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
    return (nullptr == nt_env) 
      ? tbb::task_arena::automatic  // -1
      : static_cast<size_t>(std::atoi(nt_env));
  }

  template <typename ODE>
  Rcpp::List calc_ll(std::unique_ptr<ODE> od,
                     const Rcpp::IntegerVector& ances,
                     const Rcpp::NumericMatrix& states,
                     const Rcpp::NumericMatrix& forTime,
                     const std::string& method,
                     double atol,
                     double rtol,
                     bool see_states,
                     bool use_normalization)
  {
    auto num_threads = get_rcpp_num_threads();
    auto global_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, num_threads);

    auto T0 = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> tstates{};
    for (int i = 0; i < states.nrow(); ++i) {
      tstates.emplace_back(states.row(i).begin(), states.row(i).end());
    }
    const auto phy_edge = make_phy_edge_vector(rmatrix<const double>(forTime));
    auto inodes = find_inte_nodes(phy_edge, rvector<const int>(ances), tstates);

    calc_ll_res ll_res;
    if (use_normalization) {
      ll_res  = calc_ll(Integrator<ODE, odeintcpp::normalize>(      std::move(od), method, atol, rtol), inodes, tstates);
    } else {
      ll_res = calc_ll(Integrator<ODE, odeintcpp::no_normalization>(std::move(od), method, atol, rtol), inodes, tstates);
    }
         

    auto T1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> DT = (T1 - T0);
    Rcpp::NumericMatrix states_out;
    if (see_states) {
      // R side expect full states back.
      states_out = Rcpp::NumericMatrix(states.nrow(), states.ncol());
      for (int i = 0; i < states.nrow(); ++i) {
        std::copy(std::begin(tstates[i]), std::end(tstates[i]), 
                  states_out.row(i).begin());
      }
    }
    return Rcpp::List::create(Rcpp::Named("loglik") = ll_res.loglik,
                              Rcpp::Named("node_M") = ll_res.node_M,
                              Rcpp::Named("merge_branch") = ll_res.merge_branch,
                              Rcpp::Named("states") = states_out,
                              Rcpp::Named("duration") = DT.count());
  }


  template <typename ODE>
  Rcpp::NumericVector ct_condition(std::unique_ptr<ODE> od,
                                   const Rcpp::NumericVector& y,
                                   const double t,
                                   const std::string& method,
                                   double atol,
                                   double rtol,
                                   bool use_normalization)  {
    auto init_state = std::vector<double>(y.begin(), y.end());

    if (use_normalization) {
      odeintcpp::normalize norm;
      odeintcpp::integrate(method,
                           std::move(od),
                           &init_state,         // state vector
                           0.0,                 // t0
                           t,                   // t1
                           t * 0.01,
                           atol,
                           rtol,
                           norm);
    } else {
        odeintcpp::no_normalization norm;
        odeintcpp::integrate(method,
                             std::move(od),
                             &init_state,         // state vector
                             0.0,                 // t0
                             t,                   // t1
                             t * 0.01,
                             atol,
                             rtol,
                             norm);
    }
    return Rcpp::NumericVector(init_state.begin(), init_state.end());
  }
}  // namespace  secsse


// [[Rcpp::export]]
Rcpp::List calc_ll_cpp(const std::string& rhs,
                       const Rcpp::IntegerVector& ances,
                       const Rcpp::NumericMatrix& states,
                       const Rcpp::NumericMatrix& forTime,
                       const Rcpp::RObject& lambdas,
                       const Rcpp::NumericVector& mus,
                       const Rcpp::NumericMatrix& Q,
                       const std::string& method,
                       double atol,
                       double rtol,
                       bool is_complete_tree,
                       bool see_states,
                       bool use_normalization)
{
  using namespace secsse;  // remove 'secsse::' once deprecated code is removed
  if (rhs == "ode_standard") {
    auto ll = Rcpp::as<Rcpp::NumericVector>(lambdas);
    return is_complete_tree 
      ? calc_ll(std::make_unique<ode_standard<OdeVariant::complete_tree>>(ll, mus, Q), ances, states, forTime, method, atol, rtol, see_states, use_normalization)
      : calc_ll(std::make_unique<ode_standard<OdeVariant::normal_tree>>(ll, mus, Q), ances, states, forTime, method, atol, rtol, see_states, use_normalization);
  } 
  else if (rhs == "ode_cla") {
    auto ll = Rcpp::as<Rcpp::List>(lambdas);
    
    return is_complete_tree 
        ? calc_ll(std::make_unique<ode_cla<OdeVariant::complete_tree>>(ll, mus, Q), ances, states, forTime, method, atol, rtol, see_states, use_normalization)
        : calc_ll(std::make_unique<ode_cla<OdeVariant::normal_tree>>(ll, mus, Q), ances, states, forTime, method, atol, rtol, see_states, use_normalization);
  } 
  else {
    throw std::runtime_error("calc_ll_cpp: unknown rhs");
  }
}


// [[Rcpp::export]]
Rcpp::NumericVector ct_condition_cpp(const std::string rhs,
                                     const Rcpp::NumericVector& state,
                                     const double t,
                                     const Rcpp::RObject& lambdas,
                                     const Rcpp::NumericVector& mus,
                                     const Rcpp::NumericMatrix& Q,
                                     const std::string& method,
                                     double atol,
                                     double rtol,
                                     bool use_normalization) 
{
  using namespace secsse;  // remove '::secsse::' once deprecated code is removed
  if (rhs == "ode_standard") {
    auto ll = Rcpp::as<Rcpp::NumericVector>(lambdas);
    return ct_condition(std::make_unique<ode_standard<OdeVariant::ct_condition>>(ll, mus, Q), state, t, method, atol, rtol, use_normalization);
  } 
  else if (rhs == "ode_cla") {
    auto ll = Rcpp::as<Rcpp::List>(lambdas);

    return ct_condition(std::make_unique<ode_cla<OdeVariant::ct_condition>>(ll, mus, Q), state, t, method, atol, rtol,
                        use_normalization);
  } 
  else {
    throw std::runtime_error("ct_condition_cpp: unknown rhs");
  }
}


