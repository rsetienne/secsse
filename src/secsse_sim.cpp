#include <Rcpp.h>


#include "secsse_sim.h"
#include "util.h"

//' @export
// [[Rcpp::export]]
Rcpp::List secsse_sim_cpp(const std::vector<double>& m_R,
                          const Rcpp::List& lambdas_R,
                          const Rcpp::NumericMatrix& q_R,
                          double max_time,
                          double max_species) {
  
  num_mat q; 
  numericmatrix_to_vector(q_R, q);
  
  num_mat_mat lambdas;
  list_to_vector(lambdas_R, lambdas);
  
  secsse_sim sim(m_R, 
                 lambdas,
                 q,
                 max_time,
                 max_species);
  
  while (true) {
      sim.run(-1); // use random trait.
      // check num traits
      int obs_num_traits = sim.get_num_traits();
      
      if (m_R.size() == obs_num_traits) {
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