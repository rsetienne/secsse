//
//  Copyright (c) 2022 - 2023, Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include "config.h"       // NOLINT [build/include_subdir]
#include "odeint.h"       // NOLINT [build/include_subdir]
#include "rhs.h"          // NOLINT [build/include_subdir]
#include "util.h"         // NOLINT [build/include_subdir]

#include <vector>
#include <Rcpp.h>



template <typename ODE_TYPE>
double calc_ll_cla(const Rcpp::List& ll,
                   const Rcpp::NumericVector& mm,
                   const Rcpp::NumericMatrix& Q,
                   const std::vector<int>& ances,
                   const std::vector< std::vector< double >>& for_time,
                   std::vector<std::vector<double>>* states,
                   Rcpp::NumericVector* merge_branch_out,
                   Rcpp::NumericVector* nodeM_out,
                   const std::string& method,
                   double absolute_tol,
                   double relative_tol) {
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

  ODE_TYPE od(ll_cpp, mm_cpp, Q_cpp);

  size_t d = od.get_d();

  std::vector<double> mergeBranch(d);
  std::vector<double> nodeN;
  std::vector<double> nodeM;

  int max_ances = *std::max_element(ances.begin(), ances.end());
  std::vector< double > add(states[0].size(), 0.0);
  while (max_ances > (*states).size()) {
    (*states).push_back(add);
  }
  (*states).push_back(add);

  std::vector< double > logliks(ances.size());
  std::vector<double> y;

  std::vector<int> desNodes;
  std::vector<double> timeInte;
  long double loglik = 0;
  for (int a = 0; a < ances.size(); ++a) {
   
    int focal = ances[a];
    find_desNodes(for_time, focal, &desNodes, &timeInte);

    int focal_node;
    for (size_t i = 0; i < desNodes.size(); ++i) {
      focal_node = desNodes[i];
      if (focal_node < 0) throw "focal_node < 0";
      if (focal_node >= states->size()) throw "focal_node > states.size";

      y = (*states)[focal_node];
      
  //    std::cerr << focal << " " << focal_node << " " << timeInte[i] << "\n";

      std::unique_ptr<ODE_TYPE> od_ptr = std::make_unique<ODE_TYPE>(od);
      odeintcpp::integrate(method,
                           std::move(od_ptr),     // ode class object
                           &y,                    // state vector
                           0.0,                   // t0
                           timeInte[i],           // t1
                           timeInte[i] * 0.1,
                           absolute_tol,
                           relative_tol);

      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
    }

    normalize_loglik_node(&nodeM, &loglik);
    normalize_loglik_node(&nodeN, &loglik);

  //  std::cerr << loglik << " ";
    
    mergeBranch = std::vector<double>(d, 0.0);

    for (size_t i = 0; i < d; ++i) {
      for (size_t j = 0; j < d; ++j) {
        for (size_t k = 0; k < d; ++k) {
          if (ll_cpp[i][j][k] != 0.0) {
           mergeBranch[i] += ll_cpp[i][j][k] * (nodeN[j + d] * nodeM[k + d] +
                                                nodeM[j + d] * nodeN[k + d]);
          }
        }
      }
      mergeBranch[i] *= 0.5;
    }

    normalize_loglik(&mergeBranch, &loglik);

  //  std::cerr << loglik << "\n";
    
    std::vector<double> newstate(d);
    for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());

    if (focal_node < 0) throw "focal_node < 0";
    if (focal_node >= states->size()) throw "focal_node > states.size";

    (*states)[focal] = newstate;
  }

  for (int i = 0; i < mergeBranch.size(); ++i) {
    (*merge_branch_out).push_back(mergeBranch[i]);
  }
  for (int i = 0; i < nodeM.size(); ++i) {
    (*nodeM_out).push_back(nodeM[i]);
  }
   
  return loglik;
}

// [[Rcpp::export]]
Rcpp::NumericVector ct_condition_cla(const Rcpp::NumericVector& y,
                                 double t,
                                 const Rcpp::List& ll,
                                 const Rcpp::NumericVector& mm,
                                 const Rcpp::NumericMatrix& Q,
                                 std::string method,
                                 double atol,
                                 double rtol) {
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

  ode_cla_e od(ll_cpp, mm_cpp, Q_cpp);

  std::vector<double> init_state(y.begin(), y.end());

  std::unique_ptr<ode_cla_e> od_ptr = std::make_unique<ode_cla_e>(od);
  odeintcpp::integrate(method,
                       std::move(od_ptr),    // ode class object
                       &init_state,          // state vector
                       0.0,                  // t0
                       t,                    // t1
                       t * 0.01,
                       atol,
                       rtol);

  Rcpp::NumericVector out;
  for (int i = 0; i < init_state.size(); ++i) {
    out.push_back(init_state[i]);
  }
  return out;
}


//' function to do cpp stuff
//' @param ances ances
//' @param states_R states_R
//' @param forTime_R fr
//' @param lambdas l
//' @param mus mus
//' @param Q Q
//' @param method method
//' @param atol atol
//' @param rtol rtol
//' @param is_complete_tree ss
//' @export
// [[Rcpp::export]]
Rcpp::List cla_calThruNodes_cpp(const Rcpp::NumericVector& ances,
                                const Rcpp::NumericMatrix& states_R,
                                const Rcpp::NumericMatrix& forTime_R,
                                const Rcpp::List& lambdas,
                                const Rcpp::NumericVector& mus,
                                const Rcpp::NumericMatrix& Q,
                                std::string method,
                                double atol,
                                double rtol,
                                bool is_complete_tree) {
try {
  std::vector< std::vector< double >> states, forTime;
  numericmatrix_to_vector(states_R, &states);
  numericmatrix_to_vector(forTime_R, &forTime);

  Rcpp::NumericVector mergeBranch;
  Rcpp::NumericVector nodeM;

  double loglik = 0.0;
  if (is_complete_tree) {
    loglik = calc_ll_cla< ode_cla_d >(lambdas,
                                      mus,
                                      Q,
                                      std::vector<int>(ances.begin(),
                                                       ances.end()),
                                      forTime,
                                      &states,
                                      &mergeBranch,
                                      &nodeM,
                                      method, atol, rtol);
  } else {
    loglik = calc_ll_cla< ode_cla >(lambdas,
                                    mus,
                                    Q,
                                    std::vector<int>(ances.begin(),
                                                     ances.end()),
                                    forTime,
                                    &states,
                                    &mergeBranch,
                                    &nodeM,
                                    method, atol, rtol);
  }
  Rcpp::NumericMatrix states_out;
  vector_to_numericmatrix(states, &states_out);
  Rcpp::List output = Rcpp::List::create(Rcpp::Named("states") = states_out,
                                         Rcpp::Named("loglik") = loglik,
                                         Rcpp::Named("mergeBranch") =
                                           mergeBranch,
                                         Rcpp::Named("nodeM") = nodeM);
  return output;
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}
