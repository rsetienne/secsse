#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <tuple>

#include <thread>
#include <chrono>

#include "odeint.h"
#include "util.h"

#include <cmath>

#include <RcppParallel.h>


#include "threaded_ll.h"

struct combine_states_cla {
  
  combine_states_cla(int d, const ode_cla& od) : d_(d), od_(od) {}
  
  state_vec operator()(const std::tuple< state_vec, state_vec >& input_states) {
    
    state_vec nodeN =  std::get<0>(input_states);
    state_vec nodeM =  std::get<1>(input_states);
    
    double ll1 = nodeN.back(); nodeN.pop_back();
    double ll2 = nodeM.back(); nodeM.pop_back();
    
    state_vec mergeBranch = std::vector<double>(d_, 0.0);
    
    for (size_t i = 0; i < d_; ++i) {
      for (size_t j = 0; j < d_; ++j) {
        for (size_t k = 0; k < d_; ++k) {
          
          double a = od_.get_l(i, j, k);
          
          if (a != 0.0) {
            double mult = (nodeN[j + d_] * nodeM[k + d_] +
                           nodeM[j + d_] * nodeN[k + d_]);
            
            mergeBranch[i] += a * mult;
            
            if (j+d_ > nodeN.size()) mergeBranch[i] = -1;
            if (j+d_ > nodeM.size()) mergeBranch[i] = -2;
            if (k+d_ > nodeN.size()) mergeBranch[i] = -3;
            if (k+d_ > nodeM.size()) mergeBranch[i] = -4;
            
          }
        }
      }
    }
    
    double loglik = ll1 + ll2;
    
    // next line does get called:
    normalize_loglik(mergeBranch, loglik);
    // output_vec(nodeN); // somehow doesn't get called
    // output_vec(nodeM); // somehow doesn't get called
    
    state_vec newstate(d_);
    for (int i = 0; i < d_; ++i) {
      newstate[i] = nodeM[i];
    }
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    newstate.push_back(loglik);
    
    return newstate;
  }
  
  size_t d_;
  ode_cla od_;
};

//' cla log likelihood using tbb flow
//' @description cla loglik
//' @param ances ances
//' @param states_R states
//' @param forTime_r ff
//' @param lambdas list of lambdas
//' @param mus vector of mus
//' @param Q q
//' @param num_threads num threads
//' @export
// [[Rcpp::export]]
Rcpp::List calc_cla_ll_threaded(const Rcpp::NumericVector& ances,
                                const Rcpp::NumericMatrix& states_R,
                                const Rcpp::NumericMatrix& forTime_R,
                                const Rcpp::List& lambdas_R,
                                const Rcpp::NumericVector& mus_R,
                                const Rcpp::NumericMatrix& Q,
                                int num_threads) {
  try {
    
    // let's see if the list stuff works!
    /* for (int i = 0; i < lambdas.size(); ++i) {
     NumericMatrix temp = lambdas[i];
     for (int j = 0; j < temp.nrow(); ++j) {
     for (int k = 0; k < temp.ncol(); ++k) {
     std::cerr << temp(j, k) << " ";  
     }  std::cerr << "\n";
     }
     std::cerr << "\n\n";
     std::cerr << i << "\n";
     
     }
     */
    
    std::vector< std::vector< double >> states_cpp, for_time_cpp, Q_cpp;
    numericmatrix_to_vector(states_R, states_cpp);
    numericmatrix_to_vector(forTime_R, for_time_cpp);
    numericmatrix_to_vector(Q, Q_cpp);
    
    std::vector< int > ances_cpp(ances.begin(), ances.end());
    
    std::vector<double> mus_cpp(mus_R.begin(), mus_R.end());
    
    std::vector<std::vector<std::vector<double>>> lambdas_cpp;
    
    list_to_vector(lambdas_R, lambdas_cpp);
    
    ode_cla od_(lambdas_cpp, mus_cpp, Q_cpp);
    
    threaded_ll<ode_cla, combine_states_cla> ll_calc(od_, ances_cpp, for_time_cpp, states_cpp, num_threads);
    
   // Rcout << "set up threaded_ll object\n"; force_output();
    
    return ll_calc.calc_ll();
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}