//
//  Copyright (c) 2022 - 2023, Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <thread>
#include <chrono>
#include <type_traits>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <set>

#include "config.h"       // NOLINT [build/include_subdir]
#include "odeint.h"       // NOLINT [build/include_subdir]
#include "rhs.h"          // NOLINT [build/include_subdir]
#include "util.h"         // NOLINT [build/include_subdir]

#include <vector>
#include <Rcpp.h>

namespace orig_cla {

template <typename ODE_TYPE>
double calc_ll_cla(const Rcpp::List& ll,
                   const Rcpp::NumericVector& mm,
                   const Rcpp::NumericMatrix& Q,
                   const std::vector<int>& ances,
                   const std::vector< std::vector< double >>& for_time,
                   std::vector<std::vector<double>>* states,
                   Rcpp::NumericVector* merge_branch_out,
                   Rcpp::NumericVector* nodeM_out,
                   const std::string& method,
                   double absolute_tol,
                   double relative_tol) {
  std::vector< std::vector< std::vector< double > >> ll_cpp;
  for (int i = 0; i < ll.size(); ++i) {
    Rcpp::NumericMatrix temp = ll[i];
    std::vector< std::vector< double >> temp2;
    for (int j = 0; j < temp.nrow(); ++j) {
      std::vector<double> row;
      for (int k = 0; k < temp.ncol(); ++k) {
        row.push_back(temp(j, k));
      }
      temp2.push_back(row);
    }
    ll_cpp.push_back(temp2);
  }

  std::vector<double> mm_cpp(mm.begin(), mm.end());

  std::vector< std::vector<double >> Q_cpp;
  numericmatrix_to_vector(Q, &Q_cpp);

  ODE_TYPE od(ll_cpp, mm_cpp, Q_cpp);

  size_t d = od.get_d();

  std::vector<double> mergeBranch(d);
  std::vector<double> nodeN;
  std::vector<double> nodeM;

  int max_ances = *std::max_element(ances.begin(), ances.end());
  std::vector< double > add(states[0].size(), 0.0);
  while (max_ances > static_cast<int>((*states).size())) {
    (*states).push_back(add);
  }
  (*states).push_back(add);

  std::vector< double > logliks(ances.size());
  std::vector<double> y;

  std::vector<int> desNodes;
  std::vector<double> timeInte;
  long double loglik = 0;
  for (size_t a = 0; a < ances.size(); ++a) {
   
    int focal = ances[a];
    find_desNodes(for_time, focal, &desNodes, &timeInte);

    int focal_node = 0;
    for (size_t i = 0; i < desNodes.size(); ++i) {
      focal_node = desNodes[i];
      if (focal_node < 0) throw "focal_node < 0";
      if (focal_node >= static_cast<int>(states->size())) throw "focal_node > states.size";

      y = (*states)[focal_node];
      
  //    std::cerr << focal << " " << focal_node << " " << timeInte[i] << "\n";

      std::unique_ptr<ODE_TYPE> od_ptr = std::make_unique<ODE_TYPE>(od);
      odeintcpp::integrate(method,
                           std::move(od_ptr),     // ode class object
                           &y,                    // state vector
                           0.0,                   // t0
                           timeInte[i],           // t1
                           timeInte[i] * 0.1,
                           absolute_tol,
                           relative_tol);

      if (i == 0) nodeN = y;
      if (i == 1) nodeM = y;
    }

    normalize_loglik_node(&nodeM, &loglik);
    normalize_loglik_node(&nodeN, &loglik);

  //  std::cerr << loglik << " ";
    
    mergeBranch = std::vector<double>(d, 0.0);

    for (size_t i = 0; i < d; ++i) {
      for (size_t j = 0; j < d; ++j) {
        for (size_t k = 0; k < d; ++k) {
          if (ll_cpp[i][j][k] != 0.0) {
           mergeBranch[i] += ll_cpp[i][j][k] * (nodeN[j + d] * nodeM[k + d] +
                                                nodeM[j + d] * nodeN[k + d]);
          }
        }
      }
      mergeBranch[i] *= 0.5;
    }

    normalize_loglik(&mergeBranch, &loglik);

  //  std::cerr << loglik << "\n";
    
    std::vector<double> newstate(d);
    for (size_t i = 0; i < d; ++i) newstate[i] = nodeM[i];
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());

    if (focal_node < 0) throw "focal_node < 0";
    if (focal_node >= static_cast<int>(states->size())) throw "focal_node > states.size";

    (*states)[focal] = newstate;
  }

  for (size_t i = 0; i < mergeBranch.size(); ++i) {
    (*merge_branch_out).push_back(mergeBranch[i]);
  }
  for (size_t i = 0; i < nodeM.size(); ++i) {
    (*nodeM_out).push_back(nodeM[i]);
  }
   
  return loglik;
}

}


namespace fiddled_cla {


// probably the cleanest way to retrieve RcppParallel's concurrency setting
// set by RcppParallel::setThreadOptions(numThreads)
inline size_t get_rcpp_num_threads() {
  auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  return (nullptr == nt_env) 
    ? tbb::task_arena::automatic  // -1
  : static_cast<size_t>(std::atoi(nt_env));
}


using state_ptr = std::vector<double>*;

struct des_node_t {
  state_ptr state = nullptr;
  double time = 0;   // branch length to ancestor
};

struct inte_node_t {
  state_ptr ances_state = nullptr;
  des_node_t desc[2];
};
using inte_nodes_t = std::vector<inte_node_t>;


inte_nodes_t find_inte_nodes(std::vector<std::vector<double>>& phy_edge, const std::vector<int>& ances, std::vector<std::vector<double>>* states) {
  std::sort(std::begin(phy_edge), std::end(phy_edge), [](auto& a, auto& b) {
    return a[0] < b[0];
  });
  auto comp = [](auto& edge, int val) { return edge[0] < val; };
  auto res = inte_nodes_t{ances.size()};
  for (size_t i = 0; i < ances.size(); ++i) { //tbb::parallel_for<size_t>(0, ances.size(), 1, [&](size_t i) {
    const auto focal = ances[i];
    auto& inode = res[i];
    inode.ances_state = &(*states)[focal - 1];
    // ances node shall be set to 'all NA' on the R side, 'all nan' on the C/C++ side.
    assert(std::all_of(std::begin(*inode.ances_state), std::end(*inode.ances_state), [](const auto& val) { return std::isnan(val); }));
    inode.ances_state->clear();   // NA is not nan
    
    auto it0 = std::lower_bound(std::begin(phy_edge), std::end(phy_edge), focal, comp);
    auto it1 = std::lower_bound(it0 + 1, std::end(phy_edge), focal, comp);
    assert((it0 != phy_edge.end()) && (it1 != phy_edge.end()));
    
    // easy to overlook: the sequence matters for creating the 'merged' branch.
    // imposes some pre-condition that is nowere to find :(
    if ((*it0)[1] > (*it1)[1]) {
      std::swap(*it0, *it1);
    }
    inode.desc[0] = { &(*states)[(*it0)[1] - 1], (*it0)[2] };
    inode.desc[1] = { &(*states)[(*it1)[1] - 1], (*it1)[2] };
  };
  return res;
}


template <typename RaIt>
double normalize_loglik(RaIt first, RaIt last) {
  const auto sabs = std::accumulate(first, last, 0.0, [](const auto& s, const auto& x) {
    return s + std::abs(x); 
  });
  if (sabs <= 0.0) [[unlikely]] return 0.0;
  const auto fact = 1.0 / sabs;
  for (; first != last; ++first) *first *= fact;
  return std::log(sabs);
}

// some SFINAE magic
// Primary template handles all types not supporting the operation.
template <typename, template <typename> class, typename = std::void_t<>>
struct detect : std::false_type {};
              
              // Specialization recognizes/validates only types supporting the archetype.
              template <typename T, template <typename> class Op>
              struct detect<T, Op, std::void_t<Op<T>>> : std::true_type {};
              
              template <typename OD_TYPE>
              using const_ode_callop = decltype(static_cast<void(OD_TYPE::*)(const std::vector<double>&, std::vector<double>&, const double) const>(&OD_TYPE::operator()));
              
              
              template <typename OD_TYPE>
              class Integrator {
              public:
                Integrator(std::unique_ptr<OD_TYPE>&& od, const std::string& method, double atol, double rtol) : 
                od_(std::move(od)),
                method_(method),
                atol_(atol),
                rtol_(rtol)
                {}
                
                auto operator()(std::vector<double>& state, double time) const {
                  if constexpr (detect<OD_TYPE, const_ode_callop>::value) {
                    // ode rhs is const - we can reuse
                    odeintcpp::integrate(method_,
                                         od_.get(),                  // ode class object
                                         &state,
                                         0.0,                         // t0
                                         time,                        // t1
                                         time * 0.01,                 // initial dt
                                         atol_,
                                         rtol_);
                  }
                  else {
                    // ode rhs is mutable - we must create a fresh copy
                    odeintcpp::integrate(method_,
                                         std::make_unique<OD_TYPE>(*od_.get()),  // copy
                                         &state,
                                         0.0,                         // t0
                                         time,                        // t1
                                         time * 0.01,                 // initial dt
                                         atol_,
                                         rtol_);
                  }
                }
                
              private:
                std::unique_ptr<OD_TYPE> od_;
                const std::string method_;
                const double atol_;
                const double rtol_;
              };
              
              
              template <typename OD_TYPE>
              double calc_ll_cla(const Rcpp::List& ll_R,
                             const Rcpp::NumericVector& mm_R,
                             const Rcpp::NumericMatrix& Q_R,
                             const std::vector<int>& ances,
                             std::vector< std::vector<double>>& phy_edge,  // mutable
                             std::vector<std::vector<double>>* states,
                             Rcpp::NumericVector* merge_branch_out,
                             Rcpp::NumericVector* nodeM_out,
                             double absolute_tol,
                             double relative_tol,
                             std::string method) {
                auto num_threads = get_rcpp_num_threads();
                auto global_control = tbb::global_control{tbb::global_control::max_allowed_parallelism, num_threads};
                
                auto ll_cpp = list_to_vector(ll_R);
                auto mm_cpp = std::vector<double>(mm_R.begin(), mm_R.end());
                auto q_cpp = num_mat_to_vec(Q_R);
                
                auto integrator = Integrator{std::make_unique<OD_TYPE>(ll_cpp, mm_cpp, q_cpp), method, absolute_tol, relative_tol};
                const size_t d = mm_cpp.size();
                
#ifdef __cpp_lib_atomic_float
                std::atomic<double> global_loglik{0.0};
#else
                std::mutex mutex;   // no RMW for std::atomic<double>
                double global_loglik = 0.0;
#endif
                
                auto inodes = find_inte_nodes(phy_edge, ances, states);
                auto is_dirty = [](const auto& inode) {
                  return inode.ances_state->empty() && (inode.desc[0].state->empty() || inode.desc[1].state->empty());
                };
                
                for (auto first = std::begin(inodes); first != std::end(inodes) ;) {
                  auto last = std::partition(first, std::end(inodes), std::not_fn(is_dirty));
                  tbb::parallel_for_each(first, last, [&](auto& inode) {
                    std::vector<double> y[2];
                    double loglik[2];
                    tbb::parallel_for<size_t>(0, 2, 1, [&](size_t i) {
                      auto& dnode = inode.desc[i];
                      y[i] = *dnode.state;    // copy of state vector
                      integrator(y[i], dnode.time);
                      loglik[i] = normalize_loglik(std::begin(y[i]) + d, std::end(y[i]));
                    });
                    auto& mergebranch = *inode.ances_state;
                    mergebranch.resize(2 * d);
                    for (size_t i = 0; i < d; ++i) {
                      for (size_t j = 0; j < d; ++j) {
                        for (size_t k = 0; k < d; ++k) {
                          if (ll_cpp[i][j][k] != 0.0) {
                            mergebranch[i] += ll_cpp[i][j][k] * (y[0][j + d] * y[1][k + d] +
                              y[1][j + d] * y[0][k + d]);
                          }
                        }
                      }
                      mergebranch[i] *= 0.5;
                    }
                    loglik[0] += normalize_loglik(std::begin(mergebranch) + d, std::end(mergebranch));
#ifdef __cpp_lib_atomic_float
                    global_loglik.fetch_add(inode.desc[0].time_ll + inode.desc[1].time_ll);
#else               
{
  std::lock_guard<std::mutex> _{mutex};
  global_loglik += loglik[0] + loglik[1];
}
#endif        
                  });
                  first = last;
                }
                
                const auto& root_node = inodes.back();    // the last calculted
                const auto& last_merge = *root_node.ances_state;
                (*merge_branch_out) = Rcpp::NumericVector(std::begin(last_merge) + d, std::end(last_merge));
                std::vector<double> last_M{ *root_node.desc[1].state };
                integrator(last_M, root_node.desc[1].time);
                normalize_loglik(std::begin(last_M) + d, std::end(last_M));
                (*nodeM_out) = Rcpp::NumericVector(std::begin(last_M), std::end(last_M));
                return global_loglik;
              }

}

using namespace fiddled_cla;


// [[Rcpp::export]]
Rcpp::List cla_calThruNodes_cpp(const Rcpp::NumericVector& ances,
                                const Rcpp::NumericMatrix& states_R,
                                const Rcpp::NumericMatrix& forTime_R,
                                const Rcpp::List& lambdas,
                                const Rcpp::NumericVector& mus,
                                const Rcpp::NumericMatrix& Q,
                                std::string method,
                                double atol,
                                double rtol,
                                bool is_complete_tree) {
try {
  std::vector< std::vector< double >> states, forTime;
  numericmatrix_to_vector(states_R, &states);
  numericmatrix_to_vector(forTime_R, &forTime);

  Rcpp::NumericVector mergeBranch;
  Rcpp::NumericVector nodeM;

  double loglik = 0.0;
  if (is_complete_tree) {
    loglik = calc_ll_cla< ode_cla_d >(lambdas,
                                      mus,
                                      Q,
                                      std::vector<int>(ances.begin(),
                                                       ances.end()),
                                      forTime,
                                      &states,
                                      &mergeBranch,
                                      &nodeM,
                                      atol, rtol, method);
  } else {
    loglik = calc_ll_cla< ode_cla >(lambdas,
                                    mus,
                                    Q,
                                    std::vector<int>(ances.begin(),
                                                     ances.end()),
                                    forTime,
                                    &states,
                                    &mergeBranch,
                                    &nodeM,
                                    atol, rtol, method);
  }
  Rcpp::NumericMatrix states_out;
  vector_to_numericmatrix(states, &states_out);
  Rcpp::List output = Rcpp::List::create(Rcpp::Named("states") = states_out,
                                         Rcpp::Named("loglik") = loglik,
                                         Rcpp::Named("mergeBranch") =
                                           mergeBranch,
                                         Rcpp::Named("nodeM") = nodeM);
  return output;
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}


// [[Rcpp::export]]
Rcpp::NumericVector ct_condition_cla(const Rcpp::NumericVector& y,
                                     double t,
                                     const Rcpp::List& ll,
                                     const Rcpp::NumericVector& mm,
                                     const Rcpp::NumericMatrix& Q,
                                     std::string method,
                                     double atol,
                                     double rtol) {
  std::vector< std::vector< std::vector< double > >> ll_cpp;
  for (int i = 0; i < ll.size(); ++i) {
    Rcpp::NumericMatrix temp = ll[i];
    std::vector< std::vector< double >> temp2;
    for (int j = 0; j < temp.nrow(); ++j) {
      std::vector<double> row;
      for (int k = 0; k < temp.ncol(); ++k) {
        row.push_back(temp(j, k));
      }
      temp2.push_back(row);
    }
    ll_cpp.push_back(temp2);
  }
  
  std::vector<double> mm_cpp(mm.begin(), mm.end());
  
  std::vector< std::vector<double >> Q_cpp;
  numericmatrix_to_vector(Q, &Q_cpp);
  
  ode_cla_e od(ll_cpp, mm_cpp, Q_cpp);
  
  std::vector<double> init_state(y.begin(), y.end());
  
  std::unique_ptr<ode_cla_e> od_ptr = std::make_unique<ode_cla_e>(od);
  odeintcpp::integrate(method,
                       std::move(od_ptr),    // ode class object
                       &init_state,          // state vector
                       0.0,                  // t0
                       t,                    // t1
                       t * 0.01,
                       atol,
                       rtol);
  
  Rcpp::NumericVector out;
  for (size_t i = 0; i < init_state.size(); ++i) {
    out.push_back(init_state[i]);
  }
  return out;
}
