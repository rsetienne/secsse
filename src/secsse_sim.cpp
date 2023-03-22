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
                          std::vector<double> conditioning_vec,
                          bool non_extinction,
                          bool verbose,
                          int max_tries) {
  
  //std::cerr << "loading data from R\n"; force_output();
  num_mat q; 
  numericmatrix_to_vector(q_R, q);
  
  num_mat_mat lambdas = list_to_nummatmat(lambdas_R);
  
  if (conditioning_vec[0] == -1) conditioning_vec.clear(); // "none"
  
  
  secsse_sim sim(m_R, 
                 lambdas,
                 q,
                 max_time,
                 max_species,
                 init_states,
                 non_extinction);
  std::array<int, 5> tracker = {0, 0, 0, 0, 0};
  int cnt = 0;
  while (true) {
        sim.run(); 
        sim.check_num_traits(conditioning_vec);
          
        if (sim.run_info != done) {
          cnt++;
          tracker[ sim.run_info ]++;
           if (verbose) {
             if (cnt % 1000 == 0) {
              Rcpp::Rcout << "extinct: " << tracker[extinct] << " "
                          << "large: "   << tracker[overshoot] << " " 
                          << "cond: "    << tracker[conditioning] << "\n"; 
            }
          }
        } else {
          break;
        }

        if (cnt > max_tries) {
          break;
        }
        Rcpp::checkUserInterrupt();
        if (!non_extinction && sim.run_info == extinct) break;
  }
  //extract and return
 // if (verbose) Rcpp::Rcout << "preparing to convert to R output\n";
  Rcpp::NumericMatrix ltable_for_r;
  vector_to_numericmatrix(sim.extract_ltable(), ltable_for_r);
  
  auto traits = sim.get_traits();
  auto init = sim.get_initial_state();
  
  
  
  Rcpp::List output = Rcpp::List::create( Rcpp::Named("ltable") = ltable_for_r,
                                          Rcpp::Named("traits") = traits,
                                          Rcpp::Named("initial_state") = init,
                                          Rcpp::Named("tracker") = tracker);
 // if (verbose) Rcpp::Rcout << "handing over to R\n";
  return output;
}