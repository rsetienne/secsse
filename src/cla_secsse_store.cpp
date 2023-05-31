// Copyright 2023 Thijs Janzen
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
#include "config.h"       // NOLINT [build/include_subdir]
#include "odeint.h"       // NOLINT [build/include_subdir]
#include "util.h"         // NOLINT [build/include_subdir]

#include <Rcpp.h>

storage calc_ll_cla_store_full(
    const Rcpp::List& ll,
    const Rcpp::NumericVector& mm,
    const Rcpp::NumericMatrix& Q,
    const std::vector<int>& ances,
    const std::vector< std::vector< double >>& for_time,
    const std::vector<std::vector<double>>& states,
    std::string method,
    double atol,
    double rtol,
    bool verbose)  {
  std::vector< std::vector< std::vector< double > >> ll_cpp;
  for (size_t i = 0; i < ll.size(); ++i) {
    Rcpp::NumericMatrix temp = ll[i];
    std::vector< std::vector< double >> temp2;
    for (size_t j = 0; j < temp.nrow(); ++j) {
      std::vector<double> row;
      for (size_t k = 0; k < temp.ncol(); ++k) {
        row.push_back(temp(j, k));
      }
      temp2.push_back(row);
    }
    ll_cpp.push_back(temp2);
  }

  std::vector<double> mm_cpp(mm.begin(), mm.end());

  std::vector< std::vector<double >> Q_cpp;
  numericmatrix_to_vector(Q, &Q_cpp);

  std::vector<double> y;

  std::vector<int> desNodes;
  std::vector<double> timeInte;

  storage master_storage;
  int update_freq = ances.size() / 20;
  if (update_freq < 1) update_freq = 1;
  if (verbose) Rcpp::Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcpp::Rcout << "*";

  for (int a = 0; a < ances.size(); ++a) {
    if (a % update_freq == 0) {
      if (verbose) Rcpp::Rcout << "**";
    }
    Rcpp::checkUserInterrupt();

    int focal = ances[a];

    find_desNodes(for_time, focal, &desNodes, &timeInte);

    int focal_node;
    for (int i = 0; i < desNodes.size(); ++i) {
      focal_node = desNodes[i];
      assert((focal_node) >= 0);
      assert((focal_node) < states.size());

      ode_cla_store local_od(ll_cpp, mm_cpp, Q_cpp);

      y = states[focal_node];
      std::vector< std::vector< double >> yvecs;
      std::vector<double> t_vals;

      std::unique_ptr<ode_cla_store> od_ptr =
           std::make_unique<ode_cla_store>(local_od);
      odeintcpp::integrate_full(method,
                                std::move(od_ptr),  // ode class object
                                &y,                 // state vector
                                0.0,                // t0
                                timeInte[i],        // t1
                                timeInte[i] * 0.01,
                                atol,
                                rtol,
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

storage calc_ll_cla_store(const Rcpp::List& ll,
                          const Rcpp::NumericVector& mm,
                          const Rcpp::NumericMatrix& Q,
                          const std::vector<int>& ances,
                          const std::vector< std::vector< double >>& for_time,
                          const std::vector<std::vector<double>>& states,
                          int num_steps,
                          std::string method,
                          double atol,
                          double rtol,
                          bool verbose = false)  {
  std::vector< std::vector< std::vector< double > >> ll_cpp;
  for (size_t i = 0; i < ll.size(); ++i) {
    Rcpp::NumericMatrix temp = ll[i];
    std::vector< std::vector< double >> temp2;
    for (size_t j = 0; j < temp.nrow(); ++j) {
      std::vector<double> row;
      for (size_t k = 0; k < temp.ncol(); ++k) {
        row.push_back(temp(j, k));
      }
      temp2.push_back(row);
    }
    ll_cpp.push_back(temp2);
  }

  std::vector<double> mm_cpp(mm.begin(), mm.end());

  std::vector< std::vector<double >> Q_cpp;
  numericmatrix_to_vector(Q, &Q_cpp);

  // temp, not used:
  ode_cla od(ll_cpp, mm_cpp, Q_cpp);

  std::vector<double> y;

  std::vector<int> desNodes;
  std::vector<double> timeInte;

  storage master_storage;
  int update_freq = ances.size() / 20;
  if (update_freq < 1) update_freq = 1;
  if (verbose) Rcpp::Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcpp::Rcout << "*";

  for (int a = 0; a < ances.size(); ++a) {
    if (a % update_freq == 0 && verbose) {
      Rcpp::Rcout << "**";
    }
    Rcpp::checkUserInterrupt();

    int focal = ances[a];

    find_desNodes(for_time, focal, &desNodes, &timeInte);

    int focal_node;
    for (int i = 0; i < desNodes.size(); ++i) {
      focal_node = desNodes[i];
      assert((focal_node) >= 0);
      assert((focal_node) < states.size());

      data_storage local_storage;

      ode_cla local_od(ll_cpp, mm_cpp, Q_cpp);

      double t = 0.0;
      y = states[focal_node];
      local_storage.add_entry(t, y);

      double dt = timeInte[i] * 1.0 / num_steps;
      for (int j = 0; j < num_steps; ++j) {
        std::unique_ptr<ode_cla> od_ptr = std::make_unique<ode_cla>(local_od);
        odeintcpp::integrate(method,
                             std::move(od_ptr),  // ode class object
                             &y,                 // state vector
                             bstime_t{t},                  // t0
                             bstime_t{t + dt},             // t1/
                             bstime_t{dt * 0.1},
                             atol,
                             rtol);
        t += dt;
        local_storage.add_entry(t, y);
      }

      master_storage.add_entry(focal, focal_node, local_storage);
    }
  }
  return master_storage;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cla_calThruNodes_store_cpp(
    const Rcpp::NumericVector& ances,
    const Rcpp::NumericMatrix& states_R,
    const Rcpp::NumericMatrix& forTime_R,
    const Rcpp::List& lambdas,
    const Rcpp::NumericVector& mus,
    const Rcpp::NumericMatrix& Q,
    std::string method,
    double atol,
    double rtol,
    bool is_complete_tree,
    int num_steps,
    bool verbose) {
  try {
    std::vector< std::vector< double >> states, forTime;
    numericmatrix_to_vector(states_R, &states);
    numericmatrix_to_vector(forTime_R, &forTime);

    storage found_results;

    if (num_steps > 0) {
      found_results = calc_ll_cla_store(lambdas,
                                        mus,
                                        Q,
                                        std::vector<int>(ances.begin(),
                                                         ances.end()),
                                        forTime,
                                        states,
                                        num_steps,
                                        method,
                                        atol,
                                        rtol,
                                        verbose);
    } else {
      found_results = calc_ll_cla_store_full(lambdas,
                                             mus,
                                             Q,
                                             std::vector<int>(ances.begin(),
                                                              ances.end()),
                                             forTime,
                                             states,
                                             method,
                                             atol,
                                             rtol,
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
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
