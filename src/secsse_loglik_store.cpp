#include <vector>

#include <Rcpp.h>
using namespace Rcpp;

#include "odeint.h"
#include "util.h"

template <typename OD_TYPE>
storage calc_ll(const Rcpp::NumericVector& ll,
                const Rcpp::NumericVector& mm,
                const Rcpp::NumericMatrix& Q,
                const std::vector<int>& ances,
                const std::vector< std::vector< double >>& for_time,
                std::vector<std::vector<double>>& states,
                Rcpp::NumericVector& merge_branch_out,
                Rcpp::NumericVector& nodeM_out,
                double absolute_tol,
                double relative_tol,
                std::string method,
                int num_steps) {
  
  
  size_t d = ll.size();
  
  std::vector< double > mergeBranch(d);
  std::vector< double  > nodeN;
  std::vector< double  > nodeM;
  
  
  storage master_storage;
  int update_freq = ances.size() / 20;
  if(update_freq < 1) update_freq = 1;
  Rcout << "0--------25--------50--------75--------100\n";
  Rcout << "*";
  
  
  for (int a = 0; a < ances.size(); ++a) {
    int focal = ances[a];
    
    if (a % update_freq == 0) {
      Rcout << "**";
    }
    Rcpp::checkUserInterrupt();
    
    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, desNodes, timeInte);
    
    
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
                             std::move(od_ptr), // ode class object
                             y,// state vector
                             t,// t0
                             t + dt, //t1
                             timeInte[i] * 0.01,
                             absolute_tol,
                             relative_tol); // t1
        t += dt;
        local_storage.add_entry(t, y);
      }
      master_storage.add_entry(focal, focal_node, local_storage);
    }
  }
  
  return master_storage;
}


//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calThruNodes_store_cpp(const NumericVector& ances,
                                           const NumericMatrix& states_R,
                                           const NumericMatrix& forTime_R,
                                           const NumericVector& lambdas,
                                           const NumericVector& mus,
                                           const NumericMatrix& Q,
                                           int num_threads,
                                           double abstol,
                                           double reltol,
                                           std::string method,
                                           bool is_complete_tree,
                                           int num_steps) {
  
  std::vector< std::vector< double >> states, forTime;
  
  numericmatrix_to_vector(states_R, states);
  numericmatrix_to_vector(forTime_R, forTime);
  
  NumericVector mergeBranch;
  NumericVector nodeM;
  
  storage found_results;
  
  if (is_complete_tree) {
    found_results = calc_ll<ode_standard_ct>(lambdas,
                                             mus,
                                             Q,
                                             std::vector<int>(ances.begin(), ances.end()),
                                             forTime,
                                             states,
                                             mergeBranch,
                                             nodeM,
                                             abstol,
                                             reltol,
                                             method,
                                             num_steps);
  } else {
    found_results = calc_ll<ode_standard>(lambdas,
                                          mus,
                                          Q,
                                          std::vector<int>(ances.begin(), ances.end()),
                                          forTime,
                                          states,
                                          mergeBranch,
                                          nodeM,
                                          abstol,
                                          reltol,
                                          method,
                                          num_steps);
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
  vector_to_numericmatrix(prep_mat, output);
  
  return output;
}