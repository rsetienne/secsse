#include <Rcpp.h>


#include "secsse_sim.h"
#include "util.h"


num_mat_mat list_to_nummatmat(const Rcpp::List& lambdas_R) {
  
  num_mat_mat out(lambdas_R.size());
  for (size_t m = 0; m < lambdas_R.size(); ++m) {
    Rcpp::NumericMatrix entry_R = lambdas_R[m];
    num_mat entry_cpp(entry_R.nrow(), std::vector<double>(entry_R.ncol(), 0.0));
    for(size_t i = 0; i < entry_R.nrow(); ++i) {
      for (size_t j = 0; j < entry_R.ncol(); ++j) {
        entry_cpp[i][j] = entry_R(i, j);
      }
    }
    out[m] = entry_cpp;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List secsse_sim_cpp(const std::vector<double>& m_R,
                          const Rcpp::List& lambdas_R,
                          const Rcpp::NumericMatrix& q_R,
                          double max_time,
                          double max_species,
                          const std::vector<double>& init_states,
                          std::vector<double> conditioning) {
  
  //std::cerr << "loading data from R\n"; force_output();
  num_mat q; 
  numericmatrix_to_vector(q_R, q);
  
  num_mat_mat lambdas = list_to_nummatmat(lambdas_R);
  
  if (conditioning[0] == -1) conditioning.clear(); // "none"
  
  
  secsse_sim sim(m_R, 
                 lambdas,
                 q,
                 max_time,
                 max_species,
                 init_states);

  while (true) {
      sim.run(); 
      if (sim.check_num_traits(conditioning)) {
        break;
      }
  }
  //extract and return
  
  Rcpp::NumericMatrix ltable_for_r;
  vector_to_numericmatrix(sim.extract_ltable(), ltable_for_r);
  
  auto traits = sim.get_traits();
  auto init = sim.get_initial_state();
  
  Rcpp::List output = Rcpp::List::create( Rcpp::Named("ltable") = ltable_for_r,
                                          Rcpp::Named("traits") = traits,
                                          Rcpp::Named("initial_state") = init);
  return output;
}