#include <vector>
#include "odeint.h"
#include "util.h"

#include <Rcpp.h>
using namespace Rcpp;

using vec_t = std::vector< std::vector< std::vector< double > >>;

vec_t list_to_cpp(const Rcpp::List& ll) {
  vec_t ll_cpp;
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
  return ll_cpp;
}


double calc_ll_cla_timezones(const Rcpp::List& params,
                             const Rcpp::NumericVector& crit_t,
                             const std::vector<int>& ances,
                             const std::vector< std::vector< double >>& for_time,
                             std::vector<std::vector<double>>& states,
                             Rcpp::NumericVector& merge_branch_out,
                             Rcpp::NumericVector& nodeM_out,
                             const std::string& method,
                             double absolute_tol,
                             double relative_tol) {
  
  // we first have to load all the odes
  ode_transition_compound<ode_cla> master_od;
  
  size_t d;
  vec_t ll_cpp;
  for (int i = 0; i < crit_t.size(); ++i) {
    Rcpp::List loc = params[i];
    Rcpp::List ll = loc[0];
    NumericVector mm = loc[1];
    NumericMatrix qq = loc[2];
    
    ll_cpp = list_to_cpp(ll);
    
    std::vector<double> mm_cpp(mm.begin(), mm.end());
    d = mm_cpp.size();
    
    std::vector< std::vector<double >> q_cpp;
    numericmatrix_to_vector(qq, q_cpp);
    
    ode_cla local(ll_cpp, mm_cpp, q_cpp);
    
    master_od.add_entry(local, crit_t[i]);
  }
  
  std::vector<double> mergeBranch(d);
  std::vector<double> nodeN;
  std::vector<double> nodeM;
  
  int max_ances = *std::max_element(ances.begin(), ances.end());
  std::vector< double > add(states[0].size(), 0.0);
  while (max_ances > states.size()) {
    states.push_back(add);
  }
  states.push_back(add);
  
  std::vector< double > logliks(ances.size());
  std::vector<double> y;
  
  std::vector<int> desNodes;
  std::vector<double> timeInte;
  long double loglik = 0;
  
  for (int a = 0; a < ances.size(); ++a) {
    
    int focal = ances[a];
    find_desNodes(for_time, focal, desNodes, timeInte);
    
    int focal_node;
    //  Rcpp::Rcout << a << " ";
    for (int i = 0; i < desNodes.size(); ++i) {
      focal_node = desNodes[i];
      assert((focal_node) >= 0);
      assert((focal_node) < states.size());
      
      y = states[focal_node];
      
      
      std::unique_ptr<ode_transition_compound<ode_cla>> od_ptr = 
        std::make_unique<ode_transition_compound<ode_cla>>(master_od);
      odeintcpp::integrate(method,
                           std::move(od_ptr), // ode class object
                           y, // state vector
                           0.0, // t0
                           timeInte[i], //t1
                                   timeInte[i] * 0.1,
                                   absolute_tol,
                                   relative_tol); // t1
      
      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
    }
    
    normalize_loglik_node(nodeM, loglik); //Rcout << "nodeM: " << loglik<< "\n";
    normalize_loglik_node(nodeN, loglik); //Rcout << "nodeN: " << loglik<< "\n";
    
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
    
    normalize_loglik(mergeBranch, loglik);
    
    std::vector<double> newstate(d);
    for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    
    assert((focal) >= 0);
    assert((focal) < states.size());
    states[focal] = newstate;
  }
  
  for (int i = 0; i < mergeBranch.size(); ++i) {
    merge_branch_out.push_back(mergeBranch[i]);
  }
  for (int i = 0; i < nodeM.size(); ++i) {
    nodeM_out.push_back(nodeM[i]);
  }
  
  return loglik;
}

ode_cla get_ode_cla(const Rcpp::List& params,
                    int index) {
  Rcpp::List loc = params[index];
  Rcpp::List ll = loc[0];
  NumericVector mm = loc[1];
  NumericMatrix qq = loc[2];
  
  vec_t ll_cpp = list_to_cpp(ll);
  
  std::vector<double> mm_cpp(mm.begin(), mm.end());
  std::vector< std::vector<double >> q_cpp;
  numericmatrix_to_vector(qq, q_cpp);
  
  ode_cla local(ll_cpp, mm_cpp, q_cpp);
  return local;
}

double calc_ll_cla_timezones_break(const Rcpp::List& params,
                                   const std::vector<double>& crit_t,
                                   const std::vector<int>& ances,
                                   const std::vector< std::vector< double >>& for_time,
                                   std::vector<std::vector<double>>& states,
                                   const std::vector<double>& node_heights,
                                   Rcpp::NumericVector& merge_branch_out,
                                   Rcpp::NumericVector& nodeM_out,
                                   const std::string& method,
                                   double absolute_tol,
                                   double relative_tol) {
  
  Rcpp::List loc = params[0];
  Rcpp::List ll = loc[0];
  NumericVector mm = loc[1];
  NumericMatrix qq = loc[2];
  
  vec_t ll_cpp = list_to_cpp(ll);
  
  std::vector<double> mm_cpp(mm.begin(), mm.end());
  size_t d = mm_cpp.size();

  std::vector<double> mergeBranch(d);
  std::vector<double> nodeN;
  std::vector<double> nodeM;
  
  int max_ances = *std::max_element(ances.begin(), ances.end());
  std::vector< double > add(states[0].size(), 0.0);
  while (max_ances > states.size()) {
    states.push_back(add);
  }
  states.push_back(add);
  
  std::vector< double > logliks(ances.size());
  std::vector<double> y;
  
  std::vector<int> desNodes;
  std::vector<double> timeInte;
  long double loglik = 0;
  
  for (int a = 0; a < ances.size(); ++a) {
    
    int focal = ances[a];
    find_desNodes(for_time, focal, desNodes, timeInte);
    
    int focal_node;
    double end_t = get_node_height(focal, node_heights);
    
    //  Rcpp::Rcout << a << " ";
    for (int i = 0; i < desNodes.size(); ++i) {
      focal_node = desNodes[i];
      assert((focal_node) >= 0);
      assert((focal_node) < states.size());
      
      y = states[focal_node];
      double start_t = get_node_height(focal_node, node_heights);
      
      // now we need to integrate the branch.
      // First, we have to find the integrator at the start
      int start_index = get_time_index(start_t, crit_t);
      int end_index = get_time_index(end_t, crit_t);
      
      if (start_index == end_index) {
        // normal integration along branch:
        ode_cla local = get_ode_cla(params, start_index);
        std::unique_ptr<ode_cla> od_ptr = std::make_unique<ode_cla>(local);
        
        odeintcpp::integrate(method,
                             std::move(od_ptr), // ode class object
                             y, // state vector
                             0.0, // t0
                             timeInte[i], //t1
                                     timeInte[i] * 0.1,
                                     absolute_tol,
                                     relative_tol); // t1
      } else if (end_index - start_index == 5) {
        // we traverse a transition line, this means we have to break up
        // the integration into two parts:
        
        double t0 = start_t;
        double t1 = crit_t[end_index - 1];
        double t2 = end_t; 
        // rescale these so both start at 0:
        // 
        double t1_1 = t1 - t0;
        double t2_1 = t2 - t1;
        
        ode_cla od_start = get_ode_cla(params, start_index);
        std::unique_ptr<ode_cla> od_ptr1 = std::make_unique<ode_cla>(od_start);
        odeintcpp::integrate(method,
                             std::move(od_ptr1), // ode class object
                             y, // state vector
                             0.0, // t0
                             t1_1, //t1
                             t1_1 * 0.1,
                             absolute_tol,
                             relative_tol); // t1
        // and then the second part:
        ode_cla od_end = get_ode_cla(params, end_index);
        std::unique_ptr<ode_cla> od_ptr2 = std::make_unique<ode_cla>(od_end);
        odeintcpp::integrate(method,
                             std::move(od_ptr2), // ode class object
                             y, // state vector
                             0.0, // t0
                             t2_1, //t1
                             t2_1 * 0.1,
                             absolute_tol,
                             relative_tol); // t1
      } else {
        // we traverse one or more transition lines
        std::vector< double > time_points;
        std::vector< int > param_indicators;
        for (int i = 0; i < crit_t.size(); ++i) {
          if (crit_t[i] > start_t && crit_t[i] < end_t) {
            time_points.push_back(crit_t[i]);
            param_indicators.push_back(i);
          }
        }
        time_points.push_back(end_t);
        param_indicators.push_back(end_index);
        
        double t0 = start_t;
        for (int i = 0; i < time_points.size(); ++i) {
          
          double t1 = time_points[i];
          
          ode_cla local_od = get_ode_cla(params, param_indicators[i]);

          std::unique_ptr<ode_cla> od_ptr = std::make_unique<ode_cla>(local_od);
          odeintcpp::integrate(method,
                               std::move(od_ptr), // ode class object
                               y, // state vector
                               0.0, // t0
                               t1 - t0, //t1
                               (t1 - t0) * 0.1,
                               absolute_tol,
                               relative_tol);
          t0 = t1;
        }
      }
      
      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
    }
    
    normalize_loglik_node(nodeM, loglik); //Rcout << "nodeM: " << loglik<< "\n";
    normalize_loglik_node(nodeN, loglik); //Rcout << "nodeN: " << loglik<< "\n";
    
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
    
    normalize_loglik(mergeBranch, loglik);
    
    std::vector<double> newstate(d);
    for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    
    assert((focal) >= 0);
    assert((focal) < states.size());
    states[focal] = newstate;
  }
  
  for (int i = 0; i < mergeBranch.size(); ++i) {
    merge_branch_out.push_back(mergeBranch[i]);
  }
  for (int i = 0; i < nodeM.size(); ++i) {
    nodeM_out.push_back(nodeM[i]);
  }
  
  return loglik;
}


// [[Rcpp::export]]
Rcpp::List cla_calThruNodes_timezones_cpp(const Rcpp::NumericVector& ances,
                                          const Rcpp::NumericMatrix& states_R,
                                          const Rcpp::NumericMatrix& forTime_R,
                                          const Rcpp::List& params,
                                          const Rcpp::NumericVector& crit_t,
                                          const Rcpp::NumericMatrix& node_heights_R,
                                          std::string method,
                                          double atol,
                                          double rtol,
                                          bool is_complete_tree) {
  
  try {
    std::vector< std::vector< double >> states, forTime;
    numericmatrix_to_vector(states_R, states);
    numericmatrix_to_vector(forTime_R, forTime);
    
    NumericVector mergeBranch;
    NumericVector nodeM;
    
    std::vector<double> node_heights = make_node_heights(node_heights_R);
  
    double loglik = 0.0;
    if (is_complete_tree) {
      Rcpp::stop("complete tree conditioning is not available with timezones");
    } else {
      loglik = calc_ll_cla_timezones_break(params,
                                           std::vector<double>(crit_t.begin(), crit_t.end()),
                                           std::vector<int>(ances.begin(), ances.end()),
                                           forTime,
                                           states,
                                           node_heights,
                                           mergeBranch,
                                           nodeM,
                                           method, atol, rtol);
    }
    
    NumericMatrix states_out;
    vector_to_numericmatrix(states, states_out);
    
    Rcpp::List output = Rcpp::List::create( Named("states") = states_out,
                                            Named("loglik") = loglik,
                                            Named("mergeBranch") = mergeBranch,
                                            Named("nodeM") = nodeM);
    return output;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}