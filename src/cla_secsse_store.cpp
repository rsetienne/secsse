#include <vector>
#include "odeint.h"
#include "util.h"

#include <Rcpp.h>
using namespace Rcpp;

using matrix = std::vector< std::vector< double >>;


struct data_storage {
  std::vector<double> t;
  std::vector<std::vector<double>> probs;
  
  void add_entry(double time, std::vector<double> prob) {
    t.push_back(time);
    probs.push_back(prob);
  }
};

struct entry {
  int ances;
  int focal_node;
  data_storage probabilities;
  
  entry(int a, int fn, const data_storage& probs) : 
    ances(a), focal_node(fn), probabilities(probs)
  {};
};

struct storage {
  std::vector< entry > data_;
  
  void add_entry(int a, int fn, const data_storage& p) {
    data_.push_back(entry(a, fn, p));
  }
};

double calc_prob_a(const std::vector<double>& y) {
  double sum = y[6] + y[7] + y[8] + y[9] + y[10] + y[11];
  auto rel_prob =  (y[6] + y[7] + y[8]) / sum;
  if (rel_prob < 0.0) {
    std::cerr << rel_prob << " " << "full: ";
    for (auto i : y) {
      std::cerr << i << " ";
    } std::cerr << "\n";
  }
  return rel_prob;
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
                         double rtol)  {
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
  numericmatrix_to_vector(Q, Q_cpp);
  
  // temp, not used:
  ode_cla_d od(ll_cpp, mm_cpp, Q_cpp);
  size_t d = od.get_d();
  
  std::vector<double> mergeBranch(d);
  std::vector<double> nodeN;
  std::vector<double> nodeM;
  
  std::vector< double > logliks(ances.size());
  std::vector<double> y;
  
  std::vector<int> desNodes;
  std::vector<double> timeInte;

  storage master_storage;
  int update_freq = ances.size() / 20;
  if(update_freq < 1) update_freq = 1;
  Rcout << "0--------25--------50--------75--------100\n";
  Rcout << "*";
  
 // std::cerr << "starting ances loop\n"; force_output();
  
  for (int a = 0; a < ances.size(); ++a) {
    if (a % update_freq == 0) {
      Rcout << "**";
    }
    Rcpp::checkUserInterrupt();
    
    //Rcpp::Rcout << a << " " << ances.size() << "\n"; force_output();
    
    int focal = ances[a];
    
    find_desNodes(for_time, focal, desNodes, timeInte);
    
    int focal_node;
    //  Rcpp::Rcout << a << " ";
    for (int i = 0; i < desNodes.size(); ++i) {
      focal_node = desNodes[i];
      assert((focal_node) >= 0);
      assert((focal_node) < states.size());
      
      // Rcpp::Rcout << focal_node << " " << states.size() << "\n"; force_output();
      
      y = states[focal_node];
    
      std::vector< double > dxdt(y.size(), 0.0);
      
      double t = 0.0;
      data_storage local_storage;

      ode_cla_d local_od(ll_cpp, mm_cpp, Q_cpp);
      
      double dt = timeInte[i] * 1.0 / num_steps;
      
      local_storage.add_entry(t, y);
      
      for (int j = 0; j < num_steps; ++j) {
      
        std::unique_ptr<ode_cla_d> od_ptr = std::make_unique<ode_cla_d>(local_od);
        odeintcpp::integrate(method,
                             std::move(od_ptr), // ode class object
                             y, // state vector
                             t, // t0
                             t + dt, //t1
                             dt * 0.1,
                             atol,
                             rtol); // t1

        t += dt;
       // Rcpp::Rcout << j << " adding\n"; force_output();
        local_storage.add_entry(t, y);
      }
      
      master_storage.add_entry(focal, focal_node, local_storage);
    }
  }
  return master_storage;
}

Rcpp::NumericMatrix cla_calThruNodes_store_cpp(const Rcpp::NumericVector& ances,
                                      const Rcpp::NumericMatrix& states_R,
                                      const Rcpp::NumericMatrix& forTime_R,
                                      const Rcpp::List& lambdas,
                                      const Rcpp::NumericVector& mus,
                                      const Rcpp::NumericMatrix& Q,
                                      std::string method,
                                      double atol,
                                      double rtol,
                                      bool is_complete_tree,
                                      int num_steps) {
  
  try {
    std::vector< std::vector< double >> states, forTime;
    numericmatrix_to_vector(states_R, states);
    numericmatrix_to_vector(forTime_R, forTime);
    
    NumericVector mergeBranch;
    NumericVector nodeM;
    storage found_results = calc_ll_cla_store(lambdas,
                                 mus,
                                 Q,
                                 std::vector<int>(ances.begin(), ances.end()),
                                 forTime,
                                 states,
                                 num_steps,
                                 method,
                                 atol,
                                 rtol);
     
     
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
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}