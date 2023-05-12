// Copyright 2022 - 2023 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//
#include <vector>

#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include "odeint.h"    // NOLINT [build/include_subdir]
#include "util.h"      // NOLINT [build/include_subdir]

template <typename OD_TYPE>
double calc_ll(const Rcpp::NumericVector& ll,
               const Rcpp::NumericVector& mm,
               const Rcpp::NumericMatrix& Q,
               const std::vector<int>& ances,
               const std::vector< std::vector< double >>& for_time,
               std::vector<std::vector<double>>* states,
               Rcpp::NumericVector* merge_branch_out,
               Rcpp::NumericVector* nodeM_out,
               double absolute_tol,
               double relative_tol,
               std::string method) {
  OD_TYPE od(ll, mm, Q);
  size_t d = ll.size();

  long double loglik = 0.0;

  std::vector< double > mergeBranch(d);
  std::vector< double  > nodeN;
  std::vector< double  > nodeM;

  for (int a = 0; a < ances.size(); ++a) {
    int focal = ances[a];
    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, &desNodes, &timeInte);

    for (int i = 0; i < desNodes.size(); ++i) {
      int focal_node = desNodes[i];
      std::vector< double > y = (*states)[focal_node - 1];

      std::unique_ptr<OD_TYPE> od_ptr = std::make_unique<OD_TYPE>(od);

      odeintcpp::integrate(method,
                           std::move(od_ptr),       // ode class object
                           &y,                      // state vector
                           0.0,                     // t0
                           timeInte[i],             // t1
                           timeInte[i] * 0.01,
                           absolute_tol,
                           relative_tol);
      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
    }

    normalize_loglik_node(&nodeM, &loglik);
    normalize_loglik_node(&nodeN, &loglik);

    // code correct up till here.
    for (int i = 0; i < d; ++i) {
      mergeBranch[i] = nodeM[i + d] * nodeN[i + d] * ll[i];
    }
    normalize_loglik(&mergeBranch, &loglik);

    std::vector< double > newstate(d);
    for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());

    // -1 because of R conversion to C++ indexing
    (*states)[focal - 1] = newstate;
  }

  (*merge_branch_out) = Rcpp::NumericVector(mergeBranch.begin(),
                                         mergeBranch.end());
  (*nodeM_out) = Rcpp::NumericVector(nodeM.begin(), nodeM.end());

  return loglik;
}

// [[Rcpp::export]]
Rcpp::List calThruNodes_cpp(const Rcpp::NumericVector& ances,
                            const Rcpp::NumericMatrix& states_R,
                            const Rcpp::NumericMatrix& forTime_R,
                            const Rcpp::NumericVector& lambdas,
                            const Rcpp::NumericVector& mus,
                            const Rcpp::NumericMatrix& Q,
                            int num_threads,
                            double abstol,
                            double reltol,
                            std::string method,
                            bool is_complete_tree) {
  std::vector< std::vector< double >> states, forTime;

  numericmatrix_to_vector(states_R, &states);
  numericmatrix_to_vector(forTime_R, &forTime);

  Rcpp::NumericVector mergeBranch;
  Rcpp::NumericVector nodeM;

  double loglik;
  if (is_complete_tree) {
    loglik = calc_ll<ode_standard_ct>(lambdas,
                                   mus,
                                   Q,
                                   std::vector<int>(ances.begin(), ances.end()),
                                   forTime,
                                   &states,
                                   &mergeBranch,
                                   &nodeM,
                                   abstol,
                                   reltol,
                                   method);
  } else {
    loglik = calc_ll<ode_standard>(lambdas,
                                   mus,
                                   Q,
                                   std::vector<int>(ances.begin(), ances.end()),
                                   forTime,
                                   &states,
                                   &mergeBranch,
                                   &nodeM,
                                   abstol,
                                   reltol,
                                   method);
  }
  Rcpp::NumericMatrix states_out;
  vector_to_numericmatrix(states, &states_out);

  Rcpp::List output = Rcpp::List::create(Rcpp::Named("states") = states_out,
                                         Rcpp::Named("loglik") = loglik,
                                         Rcpp::Named("mergeBranch") =
                                                      mergeBranch,
                                         Rcpp::Named("nodeM") = nodeM);
  return output;
}

// [[Rcpp::export]]
Rcpp::NumericVector ct_condition(const Rcpp::NumericVector& y,
                                 const double t,
                                 const Rcpp::NumericVector& ll,
                                 const Rcpp::NumericVector& mm,
                                 const Rcpp::NumericMatrix& Q,
                                 const std::string& method,
                                 double atol,
                                 double rtol) {
  ode_standard_ct od(ll, mm, Q);

  std::vector<double> init_state(y.begin(), y.end());

  std::unique_ptr<ode_standard_ct> od_ptr =
     std::make_unique<ode_standard_ct>(od);

  odeintcpp::integrate(method,
                       std::move(od_ptr),   // ode class object
                       &init_state,          // state vector
                       0.0,                 // t0
                       t,                   // t1
                       t * 0.01,
                       atol,
                       rtol);

  Rcpp::NumericVector out;
  for (int i = 0; i < init_state.size(); ++i) {
    out.push_back(init_state[i]);
  }
  return out;
}
