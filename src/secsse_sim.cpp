//
//  Copyright (c) 2022 - 2023, Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#include <Rcpp.h>
#include "secsse_sim2.h"   // NOLINT [build/include_subdir]
#include <string>
#include <algorithm>
#include <numeric>
#include <vector>


namespace util {  // collection of left-overs

using vec = std::vector< double >;
using nummat = Rcpp::NumericMatrix;
// Transpose Rcpp::NumericMatrix into
// std::vector<std::vector<double>>
void numericmatrix_to_vector(const nummat& m, std::vector<vec>* v) {
  (*v) = std::vector< vec >(m.nrow(), vec(m.ncol(), 0.0));
  for (int i = 0; i < m.nrow(); ++i) {
    vec row(m.ncol(), 0.0);
    for (int j = 0; j < m.ncol(); ++j) {
      row[j] = m(i, j);
    }
    (*v)[i] = row;
  }
  return;
}

void vector_to_numericmatrix(const std::vector< vec>& v, nummat* m) {
  size_t n_rows = v.size();
  size_t n_cols = v[0].size();
  (*m) = nummat(n_rows, n_cols);
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      (*m)(i, j) = v[i][j];
    }
  }
  return;
}

num_mat_mat list_to_nummatmat(const Rcpp::List& lambdas_R) {
  num_mat_mat out(lambdas_R.size());
  for (int m = 0; m < lambdas_R.size(); ++m) {
    nummat entry_R = lambdas_R[m];
    num_mat entry_cpp(entry_R.nrow(), vec(entry_R.ncol(), 0.0));
    for (int i = 0; i < entry_R.nrow(); ++i) {
      for (int j = 0; j < entry_R.ncol(); ++j) {
        entry_cpp[i][j] = entry_R(i, j);
      }
    }
    out[m] = entry_cpp;
  }
  return out;
}

}  // namespace util

// [[Rcpp::export]]
Rcpp::List secsse_sim_cpp(const std::vector<double>& m_R,
                          const Rcpp::List& lambdas_R,
                          const Rcpp::NumericMatrix& q_R,
                          double max_time,
                          double max_species,
                          bool max_species_extant,
                          double min_species,
                          const std::vector<double>& init_states,
                          std::string condition,
                          int num_concealed_states,
                          bool non_extinction,
                          bool verbose,
                          int max_tries,
                          int seed,
                          const std::vector<double>& conditioning_vec,
                          bool return_tree_size_hist,
                          bool start_at_crown) {
  try {
  num_mat q;
  util::numericmatrix_to_vector(q_R, &q);

  num_mat_mat lambdas = util::list_to_nummatmat(lambdas_R);

  secsse_sim sim(m_R,
                 lambdas,
                 q,
                 max_time,
                 max_species,
                 max_species_extant,
                 init_states,
                 non_extinction,
                 seed,
                 start_at_crown);
  std::array<double, 6> tracker = {0, 0, 0, 0, 0, 0};
  std::vector<int> tree_size_hist;
  if (return_tree_size_hist) tree_size_hist = std::vector<int>(max_tries, -1);
  int cnt = 0;
  while (true) {
    sim.run();
    cnt++;
    if (return_tree_size_hist) sim.update_tree_size_hist(&tree_size_hist[cnt]);

    if (sim.num_species() >= min_species) {
      sim.check_conditioning(condition,
                             num_concealed_states,
                             m_R.size(),
                             conditioning_vec);

      if (sim.run_info == done) {
        break;
      } else {
        tracker[ sim.run_info ]++;
      }
    } else {    // if not reached minimum size
      if (sim.run_info == extinct) {
        tracker[extinct]++;
      } else {
        tracker[ 5 ]++;
      }
    }

    if (verbose) {
      if (cnt % 1000 == 0) {
        Rcpp::Rcout << "extinct: " << tracker[extinct] << " "
                    << "large: "   << tracker[overshoot] << " "
                    << "cond: "    << tracker[conditioning] << " "
                    << "small: "   << tracker[5] << "\n";
      }
    }
    if (cnt > max_tries) {
      break;
    }
    Rcpp::checkUserInterrupt();
    if (!non_extinction && sim.run_info == extinct) break;
  }

  Rcpp::NumericMatrix ltable_for_r;      // extract and return
  util::vector_to_numericmatrix(sim.extract_ltable(), &ltable_for_r);

  auto traits = sim.get_traits();
  auto init = sim.get_initial_state();

  Rcpp::List output =
    Rcpp::List::create(Rcpp::Named("ltable") = ltable_for_r,
                       Rcpp::Named("traits") = traits,
                       Rcpp::Named("initial_state") = init,
                       Rcpp::Named("tracker") = tracker,
                       Rcpp::Named("hist_tree_size") = tree_size_hist);
  return output;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
