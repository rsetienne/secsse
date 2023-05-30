//
//  Copyright (c) 2022 - 2023, Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <Rcpp.h>

#include "odeint.h"    // NOLINT [build/include_subdir]
#include "rhs.h"       // NOLINT [build/include_subdir]
#include "util.h"      // NOLINT [build/include_subdir]

//// continuous storage
storage calc_ll_full(const Rcpp::NumericVector& ll,
                     const Rcpp::NumericVector& mm,
                     const Rcpp::NumericMatrix& Q,
                     const std::vector<int>& ances,
                     const std::vector< std::vector< double >>& for_time,
                     const std::vector<std::vector<double>>& states,
                     double absolute_tol,
                     double relative_tol,
                     std::string method,
                     bool verbose) {
  size_t d = ll.size();

  std::vector< double > mergeBranch(d);
  std::vector< double  > nodeN;
  std::vector< double  > nodeM;

  storage master_storage;
  int update_freq = ances.size() / 20;
  if (update_freq < 1) update_freq = 1;
  if (verbose) Rcpp::Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcpp::Rcout << "*";

  for (int a = 0; a < ances.size(); ++a) {
    int focal = ances[a];

    if (a % update_freq == 0 && verbose) {
      Rcpp::Rcout << "**";
    }
    Rcpp::checkUserInterrupt();

    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, &desNodes, &timeInte);

    for (int i = 0; i < desNodes.size(); ++i) {
      int focal_node = desNodes[i];

      ode_standard_store od(ll, mm, Q);

      std::vector< double > y = states[focal_node - 1];

      std::vector< std::vector< double >> yvecs;
      std::vector<double> t_vals;

      std::unique_ptr<ode_standard_store> od_ptr =
                         std::make_unique<ode_standard_store>(od);
      odeintcpp::integrate_full(method,
                                std::move(od_ptr),    // ode class object
                                &y,                   // state vector
                                0.0,                  // t0
                                timeInte[i],          // t1
                                timeInte[i] * 0.01,
                                absolute_tol,
                                relative_tol,
                                &yvecs,
                                &t_vals);

      data_storage local_storage;
      for (size_t i = 0; i < yvecs.size(); ++i) {
        local_storage.add_entry(t_vals[i], yvecs[i]);
      }
      master_storage.add_entry(focal, focal_node, local_storage);
    }
  }

  return master_storage;
}

template <typename OD_TYPE>
storage calc_ll(const Rcpp::NumericVector& ll,
                const Rcpp::NumericVector& mm,
                const Rcpp::NumericMatrix& Q,
                const std::vector<int>& ances,
                const std::vector< std::vector< double >>& for_time,
                const std::vector<std::vector<double>>& states,
                double absolute_tol,
                double relative_tol,
                std::string method,
                int num_steps,
                bool verbose) {
  size_t d = ll.size();

  std::vector< double > mergeBranch(d);
  std::vector< double  > nodeN;
  std::vector< double  > nodeM;

  storage master_storage;
  int update_freq = ances.size() / 20;
  if (update_freq < 1) update_freq = 1;
  if (verbose) Rcpp::Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcpp::Rcout << "*";

  for (int a = 0; a < ances.size(); ++a) {
    int focal = ances[a];

    if (a % update_freq == 0 && verbose) {
      Rcpp::Rcout << "**";
    }
    Rcpp::checkUserInterrupt();

    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, &desNodes, &timeInte);

    for (int i = 0; i < desNodes.size(); ++i) {
      int focal_node = desNodes[i];

      data_storage local_storage;

      OD_TYPE od(ll, mm, Q);

      double t = 0.0;
      std::vector< double > y = states[focal_node - 1];
      local_storage.add_entry(t, y);
      double dt = timeInte[i] * 1.0 / num_steps;

      for (int j = 0; j < num_steps; ++j) {
        std::unique_ptr<OD_TYPE> od_ptr = std::make_unique<OD_TYPE>(od);
        odeintcpp::integrate(method,
                             std::move(od_ptr),  // ode class object
                             &y,                 // state vector
                             t,                  // t0
                             t + dt,             // t1
                             timeInte[i] * 0.01,
                             absolute_tol,
                             relative_tol);
        t += dt;
        local_storage.add_entry(t, y);
      }

      master_storage.add_entry(focal, focal_node, local_storage);
    }
  }

  return master_storage;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix calThruNodes_store_cpp(const Rcpp::NumericVector& ances,
                                           const Rcpp::NumericMatrix& states_R,
                                           const Rcpp::NumericMatrix& forTime_R,
                                           const Rcpp::NumericVector& lambdas,
                                           const Rcpp::NumericVector& mus,
                                           const Rcpp::NumericMatrix& Q,
                                           int num_threads,
                                           double abstol,
                                           double reltol,
                                           std::string method,
                                           bool is_complete_tree,
                                           int num_steps,
                                           bool verbose) {
  std::vector< std::vector< double >> states, forTime;

  numericmatrix_to_vector(states_R, &states);
  numericmatrix_to_vector(forTime_R, &forTime);

  storage found_results;

  if (num_steps > 0) {
    if (is_complete_tree) {
      found_results = calc_ll<ode_standard_ct>(lambdas,
                                               mus,
                                               Q,
                                               std::vector<int>(ances.begin(),
                                                                ances.end()),
                                               forTime,
                                               states,
                                               abstol,
                                               reltol,
                                               method,
                                               num_steps,
                                               verbose);
    } else {
      found_results = calc_ll<ode_standard>(lambdas,
                                            mus,
                                            Q,
                                            std::vector<int>(ances.begin(),
                                                             ances.end()),
                                            forTime,
                                            states,
                                            abstol,
                                            reltol,
                                            method,
                                            num_steps,
                                            verbose);
    }
  } else {
    found_results = calc_ll_full(lambdas,
                                 mus,
                                 Q,
                                 std::vector<int>(ances.begin(), ances.end()),
                                 forTime,
                                 states,
                                 abstol,
                                 reltol,
                                 method,
                                 verbose);
  }
  std::vector< std::vector< double >> prep_mat;
  for (auto i : found_results.data_) {
    std::vector< double > add;
    for (size_t j = 0; j < i.probabilities.t.size(); ++j) {
      add = {static_cast<double>(i.ances),
             static_cast<double>(i.focal_node),
             i.probabilities.t[j]};

      for (const auto& k : i.probabilities.probs[j]) {
        add.push_back(k);
      }

      prep_mat.push_back(add);
    }
  }

  Rcpp::NumericMatrix output;
  vector_to_numericmatrix(prep_mat, &output);

  return output;
}
