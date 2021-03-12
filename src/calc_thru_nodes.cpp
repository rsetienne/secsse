#include <vector>

#include <Rcpp.h>
using namespace Rcpp;

#include "odeint.h"
#include "util.h"

/*
#include <RcppParallel.h>

double calc_ll_threaded(const Rcpp::NumericVector& ll,
                        const Rcpp::NumericVector& mm,
                        const Rcpp::NumericMatrix& Q,
                        const std::vector<int>& ances,
                        const std::vector< std::vector< double >>& for_time,
                        std::vector<std::vector<double>>& states,
                        Rcpp::NumericVector& merge_branch_out,
                        Rcpp::NumericVector& nodeM_out,
                        int num_threads);*/

double calc_ll(const Rcpp::NumericVector& ll,
               const Rcpp::NumericVector& mm,
               const Rcpp::NumericMatrix& Q,
               const std::vector<int>& ances,
               const std::vector< std::vector< double >>& for_time,
               std::vector<std::vector<double>>& states,
               Rcpp::NumericVector& merge_branch_out,
               Rcpp::NumericVector& nodeM_out) {

  MyOde od(ll, mm, Q);
  size_t d = ll.size();

  double loglik = 0.0;

  std::vector<double> mergeBranch(d);
  std::vector<double> nodeN;
  std::vector<double> nodeM;

  for (int i = 0; i < ances.size(); ++i) {
    int focal = ances[i];
    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, desNodes, timeInte);


    for (int i = 0; i < desNodes.size(); ++i) {
      int focal_node = desNodes[i];
      std::vector<double> y = states[focal_node - 1];
      bno::integrate(od, y, 0.0, timeInte[i], 0.1 * timeInte[i]);
     // odeintcpp::integrate(*od, y, 0.0, timeInte[i]);
      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
    }

    normalize_loglik_node(nodeM, loglik, d);
    normalize_loglik_node(nodeN, loglik, d);

    // code correct up till here.
    for (int i = 0; i < d; ++i) {
      mergeBranch[i] = nodeM[i + d] * nodeN[i + d] * ll[i];
    }
    normalize_loglik(mergeBranch, loglik);

    std::vector<double> newstate(d);
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
                            int num_threads) {

  std::vector< std::vector< double >> states, forTime;
  numericmatrix_to_vector(states_R, states);
  numericmatrix_to_vector(forTime_R, forTime);

  NumericVector mergeBranch;
  NumericVector nodeM;

  double loglik;
 // if (num_threads == 1) {
    loglik = calc_ll(lambdas,
                          mus,
                          Q,
                          std::vector<int>(ances.begin(), ances.end()),
                          forTime,
                          states,
                          mergeBranch,
                          nodeM);
  /*} else {
    loglik = calc_ll_threaded(lambdas,
                     mus,
                     Q,
                     std::vector<int>(ances.begin(), ances.end()),
                     forTime,
                     states,
                     mergeBranch,
                     nodeM,
                     num_threads);
  }*/

  NumericMatrix states_out;
  vector_to_numericmatrix(states, states_out);

  Rcpp::List output = Rcpp::List::create( Named("states") = states_out,
                                          Named("loglik") = loglik,
                                          Named("mergeBranch") = mergeBranch,
                                          Named("nodeM") = nodeM);
  return output;
}

/*
class phy_node {
public:
  std::vector<double> operator()(const std::vector<double>& v) {
    double loglik = v.back();
    std::vector<double> y = v;
    v.pop_back(); // remove log likelihood entry
    std::vector<double> output = calc_ll(y);
  }

private:
  MyOde od;
  std::vector<double> dt;
  std::vector<double> calc_ll(std::vector<double> states);
};
*/

/*
double calc_ll_threaded(const Rcpp::NumericVector& ll,
                        const Rcpp::NumericVector& mm,
                        const Rcpp::NumericMatrix& Q,
                        const std::vector<int>& ances,
                        const std::vector< std::vector< double >>& for_time,
                        std::vector<std::vector<double>>& states,
                        Rcpp::NumericVector& merge_branch_out,
                        Rcpp::NumericVector& nodeM_out,
                        int num_threads) {

  MyOde od(ll, mm, Q);
  size_t d = ll.size();

  double loglik = 0.0;

  std::vector<double> mergeBranch(d);
  std::vector<double> nodeN;
  std::vector<double> nodeM;

  // first, we make a tbb graph

  tbb::flow::graph g;

  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);

  for (int i = 0; i < ances.size(); ++i) {
    int focal = ances[i];
    std::vector<int> desNodes;
    std::vector<double> timeInte;
    find_desNodes(for_time, focal, desNodes, timeInte);


    for (int i = 0; i < desNodes.size(); ++i) {
      int focal_node = desNodes[i];
      std::vector<double> y = states[focal_node - 1];
      bno::integrate(od, y, 0.0, timeInte[i], 0.1 * timeInte[i]);
      // odeintcpp::integrate(*od, y, 0.0, timeInte[i]);
      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
    }

    normalize_loglik_node(nodeM, loglik, d);
    normalize_loglik_node(nodeN, loglik, d);

    for (int i = 0; i < d; ++i) {
      mergeBranch[i] = nodeM[i + d] * nodeN[i + d] * ll[i];
    }
    normalize_loglik(mergeBranch, loglik);

    std::vector<double> newstate(d);
    for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());

    states[focal - 1] = newstate; // -1 because of R conversion to C++ indexing
  }

  merge_branch_out = NumericVector(mergeBranch.begin(), mergeBranch.end());
  nodeM_out = NumericVector(nodeM.begin(), nodeM.end());

  return loglik;
}



*/
