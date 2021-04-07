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
               const MyOde& od, 
               int d) : dt_(dt), id_(id), od_(od), d_(d) {}
  
  
  state_vec operator()(const state_vec& input) {
    state_vec current_state = input;
    double loglik = current_state.back();
    current_state.pop_back();
    bno::integrate(od_, current_state, 0.0, dt_, 0.1 * dt_);
    normalize_loglik_node(current_state, loglik, d_);
    current_state.push_back(loglik);

    //std::cerr << "state_node: " << id_ << "\n";
   
    return current_state;
  }
  
  double dt_;
  int id_;
  MyOde od_;
  int d_;
};

struct combine_states {
  
  combine_states(int d, const MyOde& od) : d_(d), od_(od) {}
  
  state_vec operator()(const std::tuple< state_vec, state_vec >& input_states) {
    state_vec vec1 =  std::get<0>(input_states);
    state_vec vec2 =  std::get<1>(input_states);
    
    double ll1 = vec1.back(); vec1.pop_back();
    double ll2 = vec2.back(); vec2.pop_back();
    
    state_vec mergeBranch(d_);
    for (int i = 0; i < d_; ++i) {
      mergeBranch[i] = vec2[i + d_] * vec1[i + d_] * od_.get_l(i);
    }
    
    double loglik = ll1 + ll2;
    normalize_loglik(mergeBranch, loglik);
    
    state_vec newstate(d_);
    for (int i = 0; i < d_; ++i) {
      newstate[i] = vec2[i];
    }
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    newstate.push_back(loglik);
    
   /* for (auto i : newstate) {
      std::cerr << i << " ";
   } std::cerr << "\n";*/
    
    return newstate;
  }
  
  size_t d_;
  MyOde od_;
};

class collect_ll {
  double &my_ll;
public:
  collect_ll( double &ll ) : my_ll(ll) {}
  double operator()( const state_vec& v ) {
   // my_sum += get<0>(v) + get<1>(v);
    my_ll = v.back();
    return my_ll;
  }
};

double calc_ll_threaded_cpp2(const Rcpp::NumericVector& ll,
                            const Rcpp::NumericVector& mm,
                            const Rcpp::NumericMatrix& Q,
                            const std::vector<int>& ances,
                            const std::vector< std::vector< double >>& for_time,
                            std::vector<std::vector<double>>& states,
                            int num_threads) {
  
  // https://xinhuang.github.io/posts/2015-07-27-use-tbb-to-generate-dynamic-dependency-graph-for-computation.html
  // nodes have to be created as pointers...
  
  MyOde od(ll, mm, Q);
  size_t d = ll.size();
  
  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
  
  // first, we make a tbb graph
  
  int num_tips = ances.size() + 1;
  
  using start_node = tbb::flow::broadcast_node<state_vec>;
  using state_node = tbb::flow::function_node< state_vec, state_vec>;
  using merge_node = tbb::flow::function_node< std::tuple<state_vec, state_vec>, state_vec>;
  using join_node  = tbb::flow::join_node< std::tuple<state_vec, state_vec>, tbb::flow::queueing >;
  
  tbb::flow::graph g;
  // connect flow graph
  tbb::flow::broadcast_node<double> start(g);
  
  std::vector< state_node* > state_nodes;
  std::vector< merge_node* > merge_nodes;
  std::vector< join_node*  > join_nodes;
  
//  Rcout << "starting loading state_nodes\n"; force_output();
  
  for (int i = 0; i < states.size() + 1; ++i) {
    double dt = get_dt(for_time, i); // TODO: find dt!
    auto new_node = new state_node(g, tbb::flow::unlimited, update_state(dt, i, od, d));
    state_nodes.push_back(new_node);
  //  Rcout << i << " " << dt << "\n"; force_output();
  }
  
//p  Rcout << "starting going through ances\n"; force_output();
  
  for (int i = 0; i < ances.size(); ++i) {
    std::vector<int> connections = find_connections(for_time, ances[i]);
    
 /*   Rcout << ances[i] << " " << connections[0] << " " << connections[1] << "\n"; force_output();
    if (connections[0] > state_nodes.size()) {
      Rcout << "connections[0] > state_nodes.size()\n"; force_output();
    }
    
    if (connections[1] > state_nodes.size()) {
      Rcout << "connections[1] > state_nodes.size()\n"; force_output();
    }
    */
    auto new_join = new join_node(g);
    join_nodes.push_back(new_join);
    tbb::flow::make_edge(*state_nodes[connections[0]], std::get<0>(join_nodes.back()->input_ports()));
    tbb::flow::make_edge(*state_nodes[connections[1]], std::get<1>(join_nodes.back()->input_ports()));
    
    auto new_merge_node = new merge_node(g, tbb::flow::unlimited, combine_states(d, od));
    merge_nodes.push_back(new_merge_node);
    
    tbb::flow::make_edge(*join_nodes.back(), *merge_nodes.back());
    
  /*  if (ances[i] > state_nodes.size()) {
      Rcout << "a > state_nodes.size()\n"; force_output();
  }*/
    tbb::flow::make_edge(*merge_nodes.back(), *state_nodes[ ances[i] ]);
  }
  
  double loglik = 0.0;
  tbb::flow::function_node< state_vec, double> collect( g, tbb::flow::serial, collect_ll(loglik) );
  tbb::flow::make_edge(*merge_nodes.back(), collect);
  
  Rcout << "done creating flowgraph\n"; force_output();
  
 for (int i = 0; i < num_tips; ++i) {
   tbb::flow::broadcast_node< state_vec > input(g);
   
   tbb::flow::make_edge(input, *state_nodes[i]);
   
   std::vector<double> startvec = states[i];
   startvec.push_back(0.0);
  /* for (auto j : startvec) {
    std::cerr << j << " ";
  } std::cerr << "\n";*/
   
   input.try_put(startvec);
 }  
 
 g.wait_for_all();
//   std::cerr << loglik << "\n";
  
  return loglik;
}


//' ll threaded
//' @param ll lambdas
//' @param mm mus
//' @param Q qs
//' @param ances vector of ances
//' @param for_time fortime
//' @param states states matrix
//' @param merge_branch_out
//' @param node_M out
//' @param num_threads
//' @return log likelihood
//' @export
// [[Rcpp::export]]
double calc_ll_threaded(const Rcpp::NumericVector& ll,
                        const Rcpp::NumericVector& mm,
                        const Rcpp::NumericMatrix& Q,
                        const Rcpp::NumericVector& ances,
                        const Rcpp::NumericMatrix& for_time,
                        const Rcpp::NumericMatrix& states,
                        int num_threads) {
  try {
  std::vector< int > ances_cpp(ances.begin(), ances.end());
  
  std::vector< std::vector< double >> for_time_cpp;
  numericmatrix_to_vector(for_time, for_time_cpp);
  
  std::vector< std::vector< double >> states_cpp;
  numericmatrix_to_vector(states, states_cpp);
  
  
 // Rcout << "just before calc_ll_threaded_cpp\n"; force_output();
  
  double result = calc_ll_threaded_cpp2(ll, mm, Q, 
                       ances_cpp, for_time_cpp, states_cpp, num_threads);
  return result;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
  
  
  


/*
struct my_node {
  double t;
  my_node(double time) : t(time) {}
  
  double operator()(int input_val) {
    std::this_thread::sleep_for(std::chrono::milliseconds(30));
   // std::cerr << t << " " << input_val << "\n";
    return t * input_val;
  }
};

struct my_sum {
  double operator()( std::tuple< double, double > v ) {
      return std::get<0>(v) + std::get<1>(v);
  }
};



class my_sum2 {
  double &my_sum;
public:
  my_sum2( double &s ) : my_sum(s) {}
  double operator()( std::tuple< double, double > v ) {
    my_sum += std::get<0>(v) + std::get<1>(v);
    return my_sum;
  }
};



//' test flow graph for testing
//' @param num_threads number of threads
//' @param init_vals initial values
//' @return sum of product
//' @export
// [[Rcpp::export]]
double test_flow_graph(int num_threads,
                       const NumericVector& init_vals) {

  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
  
  double result = 0.0;
  tbb::flow::graph g; 
  // connect flow graph
  tbb::flow::broadcast_node<double> input1(g);
  tbb::flow::broadcast_node<double> input2(g);
  
  std::vector<double> times = {0.0 , 1.0};
  
  tbb::flow::function_node<double, double> node1(g, tbb::flow::unlimited, my_node(times[0]));
  tbb::flow::function_node<double, double> node2(g, tbb::flow::unlimited, my_node(times[1]));
  
  tbb::flow::join_node< std::tuple<double, double>, tbb::flow::queueing > join(g);
  tbb::flow::function_node< std::tuple<double, double>, double > summer (g, tbb::flow::serial, my_sum());
  
  tbb::flow::function_node<double, double> result_node(g, tbb::flow::unlimited, collect(result));
  
  tbb::flow::make_edge(input1, node1);
  tbb::flow::make_edge(input2, node2);
  tbb::flow::make_edge(node1, std::get<0>(join.input_ports()));
  tbb::flow::make_edge(node2, std::get<1>(join.input_ports()));
  tbb::flow::make_edge(join, summer);
  tbb::flow::make_edge(summer, result_node);
  
  input1.try_put(init_vals[0]);
  input2.try_put(init_vals[1]);
  
  g.wait_for_all();
  
  return result;
}

//' test flow graph for testing
//' @param num_threads number of threads
//' @param init_vals initial values
//' @return sum of product
//' @export
// [[Rcpp::export]]

double test_flow_graph2(int num_threads,
                       const NumericVector& init_vals) {
  
  
  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
  
  double result = 0.0;
  tbb::flow::graph g; 
  
  std::vector<double> times = {0.0 , 1.0};
  
  using internal_node = tbb::flow::function_node<double, double>;
  
  std::vector< internal_node > nodes;
  for (int i = 0; i < times.size(); ++i) {
    nodes.push_back(internal_node(g, tbb::flow::unlimited, my_node(times[i])));
  }
  tbb::flow::join_node< std::tuple<double, double>, tbb::flow::queueing > join(g);
  tbb::flow::function_node< std::tuple<double, double>, double > summer (g, tbb::flow::serial, my_sum2(result));
  
  tbb::flow::make_edge(nodes[0], std::get<0>(join.input_ports()));
  tbb::flow::make_edge(nodes[1], std::get<1>(join.input_ports()));
  tbb::flow::make_edge(join, summer);
  
  for (int i = 0; i < init_vals.size(); ++i) {
    tbb::flow::broadcast_node<double> input(g);
    tbb::flow::make_edge(input, nodes[i]);
    input.try_put(init_vals[i]);
  }
  
  g.wait_for_all();
  
  return result;
}
  
  */
  