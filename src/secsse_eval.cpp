// Copyright 2023 Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#include <cstdlib>    // std::getenv, std::atoi 
#include <vector>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]


namespace secsse {

  template <typename ODE>
  Rcpp::List eval(std::unique_ptr<ODE> od,
                  const Rcpp::IntegerVector& ances,
                  const Rcpp::NumericMatrix& states,
                  const Rcpp::NumericMatrix& forTime,
                  const std::string& method,
                  double atol,
                  double rtol,
                  size_t num_steps)
  {
    auto num_threads = get_rcpp_num_threads();
    auto global_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, num_threads);

    auto T0 = std::chrono::high_resolution_clock::now();

    // calculate valid (ancestral) states by means of calc_ll
    std::vector<std::vector<double>> tstates{};
    for (int i = 0; i < states.nrow(); ++i) {
      tstates.emplace_back(states.row(i).begin(), states.row(i).end());
    }
    const auto phy_edge = make_phy_edge_vector(rmatrix<const double>(forTime));
    auto inodes = find_inte_nodes(phy_edge, rvector<const int>(ances), tstates);
    auto integrator = Integrator<ODE>(std::move(od), method, atol, rtol);
    calc_ll(integrator, inodes, tstates);

    // integrate over each edge
    auto snodes = inodes_t<storing::inode_t>(std::begin(inodes),
                                             std::end(inodes));
    tbb::parallel_for_each(std::begin(snodes), std::end(snodes), 
                           [&](auto& snode) {
      tbb::parallel_for(0, 2, [&](size_t i) {
        integrator(snode.desc[i], num_steps);
      });
    });
    // convert to Thijs's data layout:
    // rows of [ances, focal, t, [probs]]
    const size_t nrow = 2 * snodes.size() * (num_steps + 1);
    const size_t ncol = 3 + 2 * integrator.size();
    Rcpp::NumericMatrix out(nrow, ncol);
    size_t row_index = 0;
    auto sptr_to_ridx = [&](state_ptr sptr) { 
      return static_cast<double>(std::distance(tstates.data(), sptr) + 1); };
    for (size_t i = 0; i < snodes.size(); ++i) {
      for (auto d : {0, 1}) {
        for (size_t j = 0; j < (num_steps + 1); ++j, ++row_index) {
          auto& p = snodes[i].desc[d].storage[j];
          auto row = out.row(row_index);
          row[0] = sptr_to_ridx(snodes[i].state); 
          row[1] = sptr_to_ridx(snodes[i].desc[d].state);
          row[2] = p.t;
          for (size_t k = 0; k < 2 * integrator.size(); ++k) {
            row[3 + k] = p.state[k];
          }
        }
      }
    }
    Rcpp::NumericMatrix states_out;
    states_out = Rcpp::NumericMatrix(states.nrow(), states.ncol());
    for (int i = 0; i < states.nrow(); ++i) {
      std::copy(std::begin(tstates[i]), std::end(tstates[i]), 
                states_out.row(i).begin());
    }
    auto T1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> DT = (T1 - T0);
    return Rcpp::List::create(Rcpp::Named("output") = out, 
                              Rcpp::Named("states") = states_out,
                              Rcpp::Named("duration") = DT.count());
  }

}


// [[Rcpp::export]]
Rcpp::List eval_cpp(const std::string& rhs,
                    const Rcpp::IntegerVector& ances,
                    const Rcpp::NumericMatrix& states,
                    const Rcpp::NumericMatrix& forTime,
                    const Rcpp::RObject& lambdas,
                    const Rcpp::NumericVector& mus,
                    const Rcpp::NumericMatrix& Q,
                    const std::string& method,
                    double atol,
                    double rtol,
                    bool is_complete_tree,
                    size_t num_steps)
{
  using namespace secsse;
  if (rhs == "ode_standard") {
    auto ll = Rcpp::as<Rcpp::NumericVector>(lambdas);
    return is_complete_tree 
      ? eval(std::make_unique<ode_standard<OdeVariant::complete_tree>>(ll, mus, Q), ances, states, forTime, method, atol, rtol, num_steps)
      : eval(std::make_unique<ode_standard<OdeVariant::normal_tree>>(ll, mus, Q), ances, states, forTime, method, atol, rtol, num_steps);
  } 
  else if (rhs == "ode_cla") {
    auto ll = Rcpp::as<Rcpp::List>(lambdas);
    return is_complete_tree 
      ? eval(std::make_unique<ode_cla<OdeVariant::complete_tree>>(ll, mus, Q), ances, states, forTime, method, atol, rtol, num_steps)
      : eval(std::make_unique<ode_cla<OdeVariant::normal_tree>>(ll, mus, Q), ances, states, forTime, method, atol, rtol, num_steps);
  }
  else {
    throw std::runtime_error("eval_cpp: unknown rhs");
  }
}
