#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <tuple>

#include <thread>
#include <chrono>

#include "odeint.h"
#include "util.h"


#include <RcppParallel.h>

using state_vec = std::vector<double>; 

// std::vector< std::vector< double > > ref_states;

struct update_state {
  
  update_state(double dt, int id,
               const MyOde_cla& od, 
               int d) : dt_(dt), id_(id), od_(od), d_(d) {}
  
  
  state_vec operator()(const state_vec& input) {
    state_vec current_state = input;
    double loglik = current_state.back();
    current_state.pop_back();
    bno::integrate(od_, current_state, 0.0, dt_, 0.1 * dt_);
    normalize_loglik_node(current_state, loglik, d_);
    current_state.push_back(loglik);
  
    return current_state;
  }
  
  double dt_;
  int id_;
  MyOde_cla od_;
  int d_;
};

struct combine_states {
  
  combine_states(int d, const MyOde_cla& od) : d_(d), od_(od) {}
  
  state_vec operator()(const std::tuple< state_vec, state_vec >& input_states) {
    state_vec vec1 =  std::get<0>(input_states);
    state_vec vec2 =  std::get<1>(input_states);
    
    double ll1 = vec1.back(); vec1.pop_back();
    double ll2 = vec2.back(); vec2.pop_back();
    
    state_vec mergeBranch(d_, 0.0);
    
    int ll_size = od_.get_ll_size();
    
    for (int i = 0; i < ll_size; ++i) {
      for (int j = 0; j < d_; ++j) {
        for (int k = 0; k < d_; ++k) {
          mergeBranch[i] += od_.get_l(i, j, k) * (vec1[j + d_] * vec2[k + d_] +
            vec2[j + d_] * vec1[k + d_]); // cross_M_N[j][k];
        }
      }
      mergeBranch[i] *= 0.5;
    }
    
    double loglik = ll1 + ll2;
    normalize_loglik(mergeBranch, loglik);
    
    state_vec newstate(d_);
    for (int i = 0; i < d_; ++i) {
      newstate[i] = vec2[i];
    }
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    newstate.push_back(loglik);
    
    return newstate;
  }
  
  size_t d_;
  MyOde_cla od_;
};

class collect_ll {
  state_vec &my_ll;
public:
  collect_ll( state_vec &ll ) : my_ll(ll) {}
  state_vec operator()( const state_vec& v ) {
    my_ll = v;
    return my_ll;
  }
};

Rcpp::List calc_ll_cla_threaded_cpp( const Rcpp::List& ll,
                                     const Rcpp::NumericVector& mm,
                                     const Rcpp::NumericMatrix& Q,
                                     const std::vector<int>& ances,
                                     const std::vector< std::vector< double >>& for_time,
                                     std::vector<std::vector<double>>& states,
                                     int num_threads) {
  
  // https://xinhuang.github.io/posts/2015-07-27-use-tbb-to-generate-dynamic-dependency-graph-for-computation.html
  // nodes have to be created as pointers...
  std::vector< std::vector< std::vector< double > >> ll_;
  for (int i = 0; i < ll.size(); ++i) {
    Rcpp::NumericMatrix tt = ll[i];
    std::vector< std::vector< double >> entry;
    numericmatrix_to_vector(tt, entry);
    ll_.push_back(entry);
  }
  
  MyOde_cla od(ll, mm, Q);
  size_t d = Q.size();
  
  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
  
  // first, we make a tbb graph
  
  int num_tips = ances.size() + 1;
  
  using state_node = tbb::flow::function_node< state_vec, state_vec>;
  using merge_node = tbb::flow::function_node< std::tuple<state_vec, state_vec>, state_vec>;
  using join_node  = tbb::flow::join_node< std::tuple<state_vec, state_vec>, tbb::flow::queueing >;
  
  tbb::flow::graph g;
  // connect flow graph
  tbb::flow::broadcast_node<double> start(g);
  
  std::vector< state_node* > state_nodes;
  std::vector< merge_node* > merge_nodes;
  std::vector< join_node*  > join_nodes;
  
 
  for (int i = 0; i < states.size() + 1; ++i) {
    double dt = get_dt(for_time, i); // TODO: find dt!
    auto new_node = new state_node(g, tbb::flow::unlimited, update_state(dt, i, od, d));
    state_nodes.push_back(new_node);
  }
  
  for (int i = 0; i < ances.size(); ++i) {
    std::vector<int> connections = find_connections(for_time, ances[i]);
    
    auto new_join = new join_node(g);
    join_nodes.push_back(new_join);
    tbb::flow::make_edge(*state_nodes[connections[0]], std::get<0>(join_nodes.back()->input_ports()));
    tbb::flow::make_edge(*state_nodes[connections[1]], std::get<1>(join_nodes.back()->input_ports()));
    
    auto new_merge_node = new merge_node(g, tbb::flow::unlimited, combine_states(d, od));
    merge_nodes.push_back(new_merge_node);
    
    tbb::flow::make_edge(*join_nodes.back(), *merge_nodes.back());
    
    tbb::flow::make_edge(*merge_nodes.back(), *state_nodes[ ances[i] ]);
  }
  
  state_vec output;
  tbb::flow::function_node< state_vec, state_vec> collect( g, tbb::flow::serial, collect_ll(output) );
  tbb::flow::make_edge(*merge_nodes.back(), collect);
  
  
  state_vec nodeM;
  std::vector<int> connections = find_connections(for_time, ances.back());
  tbb::flow::function_node< state_vec, state_vec> collect_nodeM( g, tbb::flow::serial, collect_ll(nodeM) );
  tbb::flow::make_edge(*state_nodes[connections[1]], collect_nodeM);
  
  nodeM.pop_back();
  
  for (int i = 0; i < num_tips; ++i) {
    tbb::flow::broadcast_node< state_vec > input(g);
    
    tbb::flow::make_edge(input, *state_nodes[i]);
    
    std::vector<double> startvec = states[i];
    startvec.push_back(0.0);
    
    input.try_put(startvec);
  }  
  
  g.wait_for_all();

  double loglikelihood = output.back();
  
  NumericVector mergeBranch;
  for (int i = d; i < (d + d); ++i) {
    mergeBranch.push_back(output[d]);
  }
  
  return Rcpp::List::create(Rcpp::Named("mergeBranch") = mergeBranch,
                            Rcpp::Named("nodeM") = nodeM,
                            Rcpp::Named("loglik") = loglikelihood);
}



// [[Rcpp::export]]
Rcpp::List calc_cla_ll_threaded(const Rcpp::NumericVector& ances,
                                const Rcpp::NumericMatrix& states_R,
                                const Rcpp::NumericMatrix& forTime_R,
                                const Rcpp::List& lambdas,
                                const Rcpp::NumericVector& mus,
                                const Rcpp::NumericMatrix& Q,
                                int num_threads) {
  try {
    
    std::vector< std::vector< double >> states_cpp, for_time_cpp;
    numericmatrix_to_vector(states_R, states_cpp);
    numericmatrix_to_vector(forTime_R, for_time_cpp);
    
    std::vector< int > ances_cpp(ances.begin(), ances.end());
    
    return calc_ll_cla_threaded_cpp(lambdas, mus, Q, ances_cpp, 
                                    for_time_cpp, states_cpp, num_threads);
    
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}