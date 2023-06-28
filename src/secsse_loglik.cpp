// Copyright 2022 - 2023 Thijs Janzen and Hanno Hildenbrandt
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//
#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <thread>
#include <chrono>
#include <type_traits>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <set>
#include "config.h"    // NOLINT [build/include_subdir]
#include "odeint.h"    // NOLINT [build/include_subdir]
#include "rhs.h"       // NOLINT [build/include_subdir]
#include "util.h"      // NOLINT [build/include_subdir]


namespace orig {

  template <typename OD_TYPE>
  double calc_ll(const Rcpp::NumericVector& ll,
                const Rcpp::NumericVector& mm,
                const Rcpp::NumericMatrix& Q,
                const std::vector<int>& ances,
                const std::vector< std::vector< double >>& for_time,
                std::vector<std::vector<double>>* states,
                Rcpp::NumericVector* merge_branch_out,
                Rcpp::NumericVector* nodeM_out,
                double absolute_tol,
                double relative_tol,
                std::string method) {
    OD_TYPE od(ll, mm, Q);
    size_t d = ll.size();
    
    long double loglik = 0.0;
    
    std::vector< double > mergeBranch(d);
    std::vector< double  > nodeN;
    std::vector< double  > nodeM;
    
    for (int a = 0; a < ances.size(); ++a) {
      int focal = ances[a];
      std::vector<int> desNodes;
      std::vector<double> timeInte;
      find_desNodes(for_time, focal, &desNodes, &timeInte);
      
      for (int i = 0; i < desNodes.size(); ++i) {
        int focal_node = desNodes[i];
        std::vector< double > y = (*states)[focal_node - 1];
        
        std::unique_ptr<OD_TYPE> od_ptr = std::make_unique<OD_TYPE>(od);
        
        odeintcpp::integrate(method,
                             std::move(od_ptr),       // ode class object
                             &y,                      // state vector
                             0.0,                     // t0
                             timeInte[i],             // t1
                             timeInte[i] * 0.01,
                             absolute_tol,
                             relative_tol);
        if (i == 0) nodeN = y;
        if (i == 1) nodeM = y;
      }
      normalize_loglik_node(&nodeM, &loglik);
      normalize_loglik_node(&nodeN, &loglik);
      
      // code correct up till here.
      for (int i = 0; i < d; ++i) {
        mergeBranch[i] = nodeM[i + d] * nodeN[i + d] * ll[i];
      }
      normalize_loglik(&mergeBranch, &loglik);
      
      std::vector< double > newstate(d);
      for (int i = 0; i < d; ++i) newstate[i] = nodeM[i];
      newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
      
      // -1 because of R conversion to C++ indexing
      (*states)[focal - 1] = newstate;
    }
    
    (*merge_branch_out) = Rcpp::NumericVector(mergeBranch.begin(),
    mergeBranch.end());
    (*nodeM_out) = Rcpp::NumericVector(nodeM.begin(), nodeM.end());
    
    return loglik;
  }

}


namespace fiddled {

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
  double calc_ll(const Rcpp::NumericVector& ll,
                 const Rcpp::NumericVector& mm,
                 const Rcpp::NumericMatrix& Q,
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
    auto integrator = Integrator{std::make_unique<OD_TYPE>(ll, mm, Q), method, absolute_tol, relative_tol};
    const size_t d = ll.size();
    
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
          mergebranch[i] =y[1][i];
          mergebranch[i + d] = y[1][i + d] * y[0][i + d] * ll[i];
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

using namespace fiddled;


// [[Rcpp::export]]
Rcpp::List calThruNodes_cpp(const Rcpp::NumericVector& ances,
                            const Rcpp::NumericMatrix& states_R,
                            const Rcpp::NumericMatrix& forTime_R,
                            const Rcpp::NumericVector& lambdas,
                            const Rcpp::NumericVector& mus,
                            const Rcpp::NumericMatrix& Q,
                            int num_threads,  // unused
                            double abstol,
                            double reltol,
                            std::string method,
                            bool is_complete_tree) {
  std::vector< std::vector< double >> states, forTime;
  
  numericmatrix_to_vector(states_R, &states);
  numericmatrix_to_vector(forTime_R, &forTime);
  
  Rcpp::NumericVector mergeBranch;
  Rcpp::NumericVector nodeM;
  
  auto T0 = std::chrono::high_resolution_clock::now();
  double loglik;
  if (is_complete_tree) {
    loglik = calc_ll<ode_standard_ct>(lambdas,
                                      mus,
                                      Q,
                                      std::vector<int>(ances.begin(), ances.end()),
                                      forTime,
                                      &states,
                                      &mergeBranch,
                                      &nodeM,
                                      abstol,
                                      reltol,
                                      method);
  } else {
    loglik = calc_ll<ode_standard>(lambdas,
                                   mus,
                                   Q,
                                   std::vector<int>(ances.begin(), ances.end()),
                                   forTime,
                                   &states,
                                   &mergeBranch,
                                   &nodeM,
                                   abstol,
                                   reltol,
                                   method);
  }
  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);
  Rcpp::NumericMatrix states_out;
  vector_to_numericmatrix(states, &states_out);
  
  Rcpp::List output = Rcpp::List::create(Rcpp::Named("states") = states_out,
                                         Rcpp::Named("loglik") = loglik,
                                         Rcpp::Named("mergeBranch") = mergeBranch,
                                         Rcpp::Named("duration") = DT.count(),
                                         Rcpp::Named("nodeM") = nodeM);
  return output;
}

// [[Rcpp::export]]
Rcpp::NumericVector ct_condition(const Rcpp::NumericVector& y,
                                 const double t,
                                 const Rcpp::NumericVector& ll,
                                 const Rcpp::NumericVector& mm,
                                 const Rcpp::NumericMatrix& Q,
                                 const std::string& method,
                                 double atol,
                                 double rtol) {
  ode_standard_ct od(ll, mm, Q);
  
  std::vector<double> init_state(y.begin(), y.end());
  
  std::unique_ptr<ode_standard_ct> od_ptr =
    std::make_unique<ode_standard_ct>(od);
  
  odeintcpp::integrate(method,
                       std::move(od_ptr),   // ode class object
                       &init_state,          // state vector
                       0.0,                 // t0
                       t,                   // t1
                       t * 0.01,
                       atol,
                       rtol);
  
  Rcpp::NumericVector out;
  for (size_t i = 0; i < init_state.size(); ++i) {
    out.push_back(init_state[i]);
  }
  return out;
}