#include <vector>
#include "odeint.h"
#include "util.h"

#include <Rcpp.h>
using namespace Rcpp;

double calc_ll_cla(const Rcpp::List& ll,
                   const Rcpp::NumericVector& mm,
                   const Rcpp::NumericMatrix& Q,
                   const std::vector<int>& ances,
                   const std::vector< std::vector< double >>& for_time,
                   std::vector<std::vector<double>>& states,
                   Rcpp::NumericVector& merge_branch_out,
                   Rcpp::NumericVector& nodeM_out) {

  std::vector< std::vector< std::vector< double > >> ll_;
  for (int i = 0; i < ll.size(); ++i) {
    Rcpp::NumericMatrix tt = ll[i];
    std::vector< std::vector< double >> entry;
    numericmatrix_to_vector(tt, entry);
    ll_.push_back(entry);
  }


  MyOde_cla od(ll, mm, Q);

  size_t d = ll.size();


  std::vector<double> mergeBranch(d);
  std::vector<double> nodeN;
  std::vector<double> nodeM;

  int max_ances = *std::max_element(ances.begin(), ances.end());
  while (max_ances > states.size()) {
    std::vector< double > add(states[0].size(), 0.0);
    states.push_back(add);
  }

  std::vector< double > logliks(ances.size());
  for (int a = 0; a < ances.size(); ++a) {

    double loglik = 0.0;

    int focal = ances[a];
    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, desNodes, timeInte);

    std::vector<double> y;
    int focal_node;
    for (int i = 0; i < desNodes.size(); ++i) {
      focal_node = desNodes[i];

      if (focal_node - 1 < 0) {
        Rcout << focal_node - 1 << "\n"; force_output();
        Rcpp::stop("focal_node - 1 < 0");
      }
      if (focal_node - 1 > states.size()) {
        Rcout << focal_node - 1 << "\n"; force_output();
        Rcpp::stop("focal_node - 1 > states.size()");
      }
      y = states[focal_node - 1];

      bno::integrate(od, y, 0.0, timeInte[i], 0.1 * timeInte[i]);

      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
    }

    normalize_loglik_node(nodeM, loglik, d);
    normalize_loglik_node(nodeN, loglik, d);

    mergeBranch = std::vector<double>(d, 0.0);
    for (int i = 0; i < ll_.size(); ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          mergeBranch[i] += ll_[i][j][k] * (nodeN[j + d] * nodeM[k + d] +
                                            nodeM[j + d] * nodeN[k + d]); // cross_M_N[j][k];
        }
      }
      mergeBranch[i] *= 0.5;
    }

    normalize_loglik(mergeBranch, loglik);


    std::vector<double> newstate(d);
    for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());

    if (focal - 1 < 0) {
      Rcout << "focal - 1 < 0"; force_output();
      Rcpp::stop("focal - 1 < 0");
    }
    if (focal - 1 > states.size()) {
      Rcout << "focal - 1 > states.size()"; force_output();
      Rcpp::stop("focal - 1 > states.size()");
    }

    states[focal - 1] = newstate; // -1 because of R conversion to C++ indexing
    logliks[a] = loglik;
  }

  merge_branch_out = NumericVector(mergeBranch.begin(), mergeBranch.end());
  nodeM_out = NumericVector(nodeM.begin(), nodeM.end());

  auto sum_loglik = std::accumulate(logliks.begin(), logliks.end(), 0.0);

  return sum_loglik;
}

// [[Rcpp::export]]
Rcpp::List cla_calThruNodes_cpp(const Rcpp::NumericVector& ances,
                                const Rcpp::NumericMatrix& states_R,
                                const Rcpp::NumericMatrix& forTime_R,
                                const Rcpp::List& lambdas,
                                const Rcpp::NumericVector& mus,
                                const Rcpp::NumericMatrix& Q) {

try {
  std::vector< std::vector< double >> states, forTime;
  numericmatrix_to_vector(states_R, states);
  numericmatrix_to_vector(forTime_R, forTime);

  NumericVector mergeBranch;
  NumericVector nodeM;

 // Rcout << "welcome into cla_calThruNodes_cpp\n"; force_output();

  double loglik = calc_ll_cla(lambdas,
                     mus,
                     Q,
                     std::vector<int>(ances.begin(), ances.end()),
                     forTime,
                     states,
                     mergeBranch,
                     nodeM);

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