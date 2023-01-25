#include <vector>
#include <cassert>
#include <Rcpp.h>
using namespace Rcpp;

#include "odeint.h"
#include "util.h"


double calc_ll_timezone(const Rcpp::List& params,
                        const NumericVector& crit_t,
                        const std::vector<int>& ances,
                        const std::vector< std::vector< double >>& for_time,
                        std::vector<std::vector<double>>& states,
                        Rcpp::NumericVector& merge_branch_out,
                        Rcpp::NumericVector& nodeM_out,
                        double absolute_tol,
                        double relative_tol,
                        std::string method) {
  
  // we first have to load all the odes
  ode_transition_compound<ode_standard> master_od;
  
   for (int i = 0; i < crit_t.size(); ++i) {
    Rcpp::List loc = params[i];
    NumericVector ll = loc[0];
    NumericVector mm = loc[1];
    NumericMatrix qq = loc[2];
    assert(ll.size() > 0);
    assert(mm.size() > 0);
    assert(qq.nrow() > 0);
    
    ode_standard local(ll, mm, qq);
    
    master_od.add_entry(local, crit_t[i]);
  }
  
  // std::cerr << "size = " << master_od.get_size() << "\n";
  
  Rcpp::List loc = params[0];
  NumericVector ll1_r = loc[0];
  
  std::vector<double> ll1(ll1_r.begin(), ll1_r.end());
  
  size_t d = ll1.size();
  
  long double loglik = 0.0;
  
  std::vector< double > mergeBranch(d);
  std::vector< double  > nodeN;
  std::vector< double  > nodeM;
  
  for (int a = 0; a < ances.size(); ++a) {
    int focal = ances[a];
    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, desNodes, timeInte);
    
    for (int i = 0; i < desNodes.size(); ++i) {
      int focal_node = desNodes[i];
      std::vector< double > y = states[focal_node - 1];
   //   std::cerr << timeInte[i] << " ";
      
      std::unique_ptr<ode_transition_compound<ode_standard>> od_ptr = 
        std::make_unique<ode_transition_compound<ode_standard>>(master_od);
      
      odeintcpp::integrate(method, 
                           std::move(od_ptr), // ode class object
                           y,// state vector
                           0.0,// t0
                           timeInte[i], //t1
                           timeInte[i] * 0.01,
                           absolute_tol,
                           relative_tol); // t1
      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
      
      for (auto i : y) {
        if (std::isnan(i)) {
          Rcpp::stop("found nan");
        }
      }
    }
    
    normalize_loglik_node(nodeM, loglik);
    normalize_loglik_node(nodeN, loglik);
    
    // code correct up till here.
    for (int i = 0; i < d; ++i) {
      mergeBranch[i] = nodeM[i + d] * nodeN[i + d] * ll1[i];
    }
    normalize_loglik(mergeBranch, loglik);
    
    std::vector< double > newstate(d);
    for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    
    states[focal - 1] = newstate; // -1 because of R conversion to C++ indexing
  }
  
  merge_branch_out = NumericVector(mergeBranch.begin(), mergeBranch.end());
  nodeM_out = NumericVector(nodeM.begin(), nodeM.end());
  
  return loglik;
}

ode_standard make_ode_standard(const Rcpp::List& params,
                               int index) {
  Rcpp::List loc = params[index];
  NumericVector ll = loc[0];
  NumericVector mm = loc[1];
  NumericMatrix qq = loc[2];
  ode_standard local(ll, mm, qq);
  return local;
}


double calc_ll_timezone_break(const Rcpp::List& params,
                              const std::vector<double>& crit_t,
                              const std::vector<int>& ances,
                              const std::vector< std::vector< double >>& for_time,
                              std::vector<std::vector<double>>& states,
                              const std::vector<double>& node_heights,
                              Rcpp::NumericVector& merge_branch_out,
                              Rcpp::NumericVector& nodeM_out,
                              double absolute_tol,
                              double relative_tol,
                              std::string method) {
  
  long double loglik = 0.0;
  
  Rcpp::List loc = params[0];
  NumericVector ll1_r = loc[0];
  std::vector<double> ll1(ll1_r.begin(), ll1_r.end());
  int d = ll1.size();
  
  std::vector< double > mergeBranch(d);
  std::vector< double  > nodeN;
  std::vector< double  > nodeM;
  
  for (int a = 0; a < ances.size(); ++a) {
    int focal = ances[a];
    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, desNodes, timeInte);
    double end_t = get_node_height(focal, node_heights);
    
    for (int i = 0; i < desNodes.size(); ++i) {
      int focal_node = desNodes[i];
     
      std::vector< double > y = states[focal_node - 1];
    
      double start_t = get_node_height(focal_node, node_heights);
      // timeInte[i] = end_t - start_t
      
      // now we need to integrate the branch.
      // First, we have to find the integrator at the start
      int start_index = get_time_index(start_t, crit_t);
      int end_index = get_time_index(end_t, crit_t);
      
      if (start_index == end_index) {
        // simple integration, easy peasy:
        ode_standard od = make_ode_standard(params, start_index);
        std::unique_ptr<ode_standard> od_ptr = std::make_unique<ode_standard>(od);
        odeintcpp::integrate(method, 
                             std::move(od_ptr), // ode class object
                             y,// state vector
                             0.0,// t0
                             timeInte[i], //t1
                             timeInte[i] * 0.01,
                             absolute_tol,
                             relative_tol); // t1
        
      }  else {
        // we traverse multiple timezones
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
          
          ode_standard local_od = make_ode_standard(params, param_indicators[i]);
          
          std::unique_ptr<ode_standard> od_ptr = std::make_unique<ode_standard>(local_od);
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
      
      for (auto i : y) {
        if (std::isnan(i)) {
          Rcpp::stop("found nan");
        }
      }
    }
    
    normalize_loglik_node(nodeM, loglik);
    normalize_loglik_node(nodeN, loglik);
    
    // code correct up till here.
    for (int i = 0; i < d; ++i) {
      mergeBranch[i] = nodeM[i + d] * nodeN[i + d] * ll1[i];
    }
    normalize_loglik(mergeBranch, loglik);
    
    std::vector< double > newstate(d);
    for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    
    states[focal - 1] = newstate; // -1 because of R conversion to C++ indexing
  }
  
  merge_branch_out = NumericVector(mergeBranch.begin(), mergeBranch.end());
  nodeM_out = NumericVector(nodeM.begin(), nodeM.end());
  
  return loglik;
}

// [[Rcpp::export]]
Rcpp::List calThruNodes_timezones_cpp(const NumericVector& ances,
                                      const NumericMatrix& states_R,
                                      const NumericMatrix& forTime_R,
                                      const List& params,
                                      const NumericVector& crit_t,
                                      const NumericMatrix& node_heights_R,
                                      int num_threads,
                                      double abstol,
                                      double reltol,
                                      std::string method,
                                      bool is_complete_tree) {
  
  
  std::vector< std::vector< double >> states, forTime;
  
  numericmatrix_to_vector(states_R, states);
  numericmatrix_to_vector(forTime_R, forTime);
  
  NumericVector mergeBranch;
  NumericVector nodeM;
  
  std::vector<double> node_heights = make_node_heights(node_heights_R);
  
//  std::cerr << "welcome! we are incpp now\n";
  
  double loglik;
  if (is_complete_tree) {
    Rcpp::stop("Complete tree conditioning is not yet implemented for timezones");
  } else {
    loglik = calc_ll_timezone_break(params,
                                    std::vector<double>(crit_t.begin(), crit_t.end()),
                                   std::vector<int>(ances.begin(), ances.end()),
                                   forTime,
                                   states,
                                   node_heights,
                                   mergeBranch,
                                   nodeM,
                                   abstol,
                                   reltol,
                                   method);
  }
  
  NumericMatrix states_out;
  vector_to_numericmatrix(states, states_out);
  
  Rcpp::List output = Rcpp::List::create( Named("states") = states_out,
                                          Named("loglik") = loglik,
                                          Named("mergeBranch") = mergeBranch,
                                          Named("nodeM") = nodeM);
  return output;
}