#include <vector>

#include <Rcpp.h>
using namespace Rcpp;

#include "odeint.h"
#include "util.h"

double calc_ll(const Rcpp::NumericVector& ll,
               const Rcpp::NumericVector& mm,
               const Rcpp::NumericMatrix& Q,
               const std::vector<int>& ances,
               const std::vector< std::vector< double >>& for_time,
               std::vector<std::vector<double>>& states,
               Rcpp::NumericVector& merge_branch_out,
               Rcpp::NumericVector& nodeM_out,
               double absolute_tol,
               double relative_tol,
               std::string method) {

  ode_standard od(ll, mm, Q);
  size_t d = ll.size();

  long double loglik = 0.0;

  std::vector< double > mergeBranch(d);
  std::vector< double  > nodeN;
  std::vector< double  > nodeM;

  for (int i = 0; i < ances.size(); ++i) {
    int focal = ances[i];
    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, desNodes, timeInte);


    for (int i = 0; i < desNodes.size(); ++i) {
      int focal_node = desNodes[i];
      std::vector< double > y = states[focal_node - 1];
      
      std::unique_ptr<ode_standard> od_ptr = std::make_unique<ode_standard>(od);
     
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
    }

    normalize_loglik_node(nodeM, loglik);
    normalize_loglik_node(nodeN, loglik);

    // code correct up till here.
    for (int i = 0; i < d; ++i) {
      mergeBranch[i] = nodeM[i + d] * nodeN[i + d] * ll[i];
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
Rcpp::List calThruNodes_cpp(const NumericVector& ances,
                            const NumericMatrix& states_R,
                            const NumericMatrix& forTime_R,
                            const NumericVector& lambdas,
                            const NumericVector& mus,
                            const NumericMatrix& Q,
                            int num_threads,
                            double abstol,
                            double reltol,
                            std::string method) {

  std::vector< std::vector< double >> states, forTime;
  numericmatrix_to_vector(states_R, states);
  numericmatrix_to_vector(forTime_R, forTime);

  NumericVector mergeBranch;
  NumericVector nodeM;

  double loglik = calc_ll(lambdas,
                          mus,
                          Q,
                          std::vector<int>(ances.begin(), ances.end()),
                          forTime,
                          states,
                          mergeBranch,
                          nodeM,
                          abstol,
                          reltol,
                          method);

  NumericMatrix states_out;
  vector_to_numericmatrix(states, states_out);

  Rcpp::List output = Rcpp::List::create( Named("states") = states_out,
                                          Named("loglik") = loglik,
                                          Named("mergeBranch") = mergeBranch,
                                          Named("nodeM") = nodeM);
  return output;
}
