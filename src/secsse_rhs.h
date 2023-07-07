//  Copyright (c) 2021 - 2023, Thijs Janzen
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#include <Rcpp.h>
#include <RcppParallel.h>
#include <type_traits>
#include <vector>


namespace secsse {

  template <typename T> using const_rvector = RcppParallel::RVector<const T>;
  template <typename T> using const_rmatrix = RcppParallel::RMatrix<const T>;
  template <typename T> using const_rmatrix_row = typename const_rmatrix<T>::Row;
  template <typename T> using const_rmatrix_col = typename const_rmatrix<T>::Column;

  template <typename T> using mutable_rvector = RcppParallel::RVector<T>;
  template <typename T> using mutable_rmatrix = RcppParallel::RMatrix<T>;
  template <typename T> using mutable_rmatrix_row = typename mutable_rmatrix<T>::Row;
  template <typename T> using mutable_rmatrix_col = typename mutable_rmatrix<T>::Column;


  // some SFINAE magic
  template <typename, template <typename> class, typename = std::void_t<>>
  struct detect : std::false_type {};
                
  template <typename T, template <typename> class Op>
  struct detect<T, Op, std::void_t<Op<T>>> : std::true_type {};
  
  template <typename ODE>
  using const_rhs_callop = decltype(static_cast<void(ODE::*)(const std::vector<double>&, std::vector<double>&, const double) const>(&ODE::operator()));
  

  enum class OdeVariant {
    normal_tree,
    complete_tree,
    ct_condition
  };


  template <OdeVariant variant>
  class ode_standard {
    const_rvector<double> l_;
    const_rvector<double> m_;
    const_rmatrix<double> q_;

  public:
    ode_standard(const Rcpp::NumericVector& l,
                 const Rcpp::NumericVector& m,
                 const Rcpp::NumericMatrix& q)
    : l_(l), m_(m), q_(q) {
    }

    size_t size() const noexcept { return l_.size(); }

    void mergebranch(const std::vector<double>& N, const std::vector<double>& M, std::vector<double>& out) const {
      const auto d = size();
      assert(2 * d == out.size());
      for (size_t i = 0; i < d; ++i) {
        out[i] = M[i];
        out[i + d] = M[i + d] * N[i + d] * l_[i];
      }
    }

    void operator()(const std::vector<double> &x,
                   std::vector<double> &dxdt,   // NOLINT [runtime/references]
                   const double /* t */) const
    {
      const auto d = size();
      if  constexpr (variant == OdeVariant::normal_tree) {
        // normal tree
        for (size_t i = 0; i < d; ++i) {
          const double t0 = l_[i] + m_[i];
          const double t1 = l_[i] * x[i];
          double dx0 = m_[i] + (t1 - t0) * x[i];
          double dxd = (2 * t1 - t0) * x[i + d];
          auto q = q_.row(i);
          for (size_t j = 0; j < d; ++j) {
            dx0 += (x[j] - x[i]) * q[j];
            dxd += (x[j + d] - x[i + d]) * q[j];
          }
          dxdt[i] = dx0;
          dxdt[i + d] = dxd;
        }
      }
      else if constexpr (variant == OdeVariant::complete_tree || variant == OdeVariant::ct_condition) {
        // complete tree including extinct branches or conditioning
        for (size_t i = 0; i < d; ++i) {
          double dx0 = (m_[i] - (l_[i] * x[i])) * (1 - x[i]);
          double dxd = -(l_[i] + m_[i]) * x[i + d];
          auto q = q_.row(i);
          for (size_t j = 0; j < d; ++j) {
            dx0 += (x[j] - x[i]) * q[j];
            dxd += (x[j + d] - x[i + d]) * q[j];
          }
          dxdt[i] = dx0;
          dxdt[i + d] = dxd;
        }
      }
    }
  };


  namespace {

    struct cla_precomp_t {
      std::vector<std::vector<std::vector<double>>> ll;
      std::vector<std::vector<std::pair<size_t, size_t>>> kb;
      std::vector<double> lambda_sum;
    };

    auto ode_cla_precomp(const Rcpp::List& Rll) {
      auto res = cla_precomp_t{};
      for (int i = 0; i < Rll.size(); ++i) {
        // we all love deeply nested loops...
        const_rmatrix<double> mr(Rcpp::as<Rcpp::NumericMatrix>(Rll[i]));
        auto& mc = res.ll.emplace_back();
        auto& kbm = res.kb.emplace_back();
        auto& ls = res.lambda_sum.emplace_back(0.0);
        for (size_t j = 0; j < mr.nrow(); ++j) {
          mc.emplace_back(mr.row(j).begin(), mr.row(j).end());
          auto& b = kbm.emplace_back(0, mc[j].size());
          for (; (mc[j][b.first] == 0.0) && (b.first <= b.second); ++b.first);        // first non-zero
          for (; (mc[j][b.second - 1] == 0.0) && (b.second > b.first); --b.second);   // last non-zero
          for (size_t k = 0; k < mc[j].size(); ++k) {
            ls += mc[j][k];
          }
        }
      }
      return res;
    }

  }


  template <OdeVariant variant>
  class ode_cla {
    // used for normal tree
    const const_rvector<double> m_;
    const const_rmatrix<double> q_;
    const cla_precomp_t prec_;

  public:
    ode_cla(const Rcpp::List ll,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericMatrix& q)
    : m_(m), q_(q), prec_(ode_cla_precomp(ll)) {
    }

    size_t size() const noexcept { return m_.size(); }

    void mergebranch(const std::vector<double>& N, const std::vector<double>& M, std::vector<double>& out) const {
      const auto d = size();
      assert(2 * d == out.size());
      for (size_t i = 0; i < d; ++i) {
        out[i] = M[i];
        out[i + d] = 0.0;
        for (size_t j = 0; j < d; ++j) {
          for (size_t k = 0; k < d; ++k) {
            out[i + d] += prec_.ll[i][j][k] * (N[j + d] * M[k + d] + M[j + d] * N[k + d]);
          }
        }
        out[i + d] *= 0.5;
      }
    }

    void operator()(const std::vector<double> &x,
                  std::vector<double> &dxdt,
                  const double /* t */) const
    {
      const auto d = size();
      if constexpr (variant == OdeVariant::normal_tree) {
        for (size_t i = 0; i < d; ++i) {
          double dx0 = 0.0;
          double dxd = 0.0;
          auto q = q_.row(i);
          const auto& kb = prec_.kb[i];
          for (size_t j = 0; j < d; ++j) {
            for (size_t k = kb[j].first; k < kb[j].second; ++k) {
              const double ll = prec_.ll[i][j][k];
              dx0 += ll * (x[j] * x[k]);
              dxd += ll * (x[j] * x[k + d] + x[j + d] * x[k]);
            }
            dx0 += (x[j] - x[i]) * q[j];
            dxd += (x[j + d] - x[i + d]) * q[j];
          }
          dxdt[i] = dx0 + m_[i] - (prec_.lambda_sum[i] + m_[i]) * x[i];
          dxdt[i + d] = dxd - (prec_.lambda_sum[i] + m_[i]) * x[i + d];
        }
      }
      else if constexpr (variant == OdeVariant::complete_tree) {
        // complete tree including extinct branches
        for (size_t i = 0; i < d; ++i) {
          double dxd = -(prec_.lambda_sum[i] + m_[i]) * x[i + d];
          auto q = q_.row(i);
          for (size_t j = 0; j < d; ++j) {
            dxd += (x[j + d] - x[i + d]) * q[j];
          }
          dxdt[i + d] = dxd;
        }
      }
      else if constexpr (variant == OdeVariant::ct_condition) {
        for (size_t i = 0; i < d; ++i) {
          double dx0 = m_[i] * (1 - x[i]);
          auto q = q_.row(i);
          const auto& kb = prec_.kb[i];
          for (size_t j = 0; j < d; ++j) {
            dx0 += (x[j] - x[i]) * q[j];
            for (size_t k = kb[j].first; k < kb[j].second; ++k) {
              dx0 += prec_.ll[i][j][k] * (x[j] * x[k] - x[i]);
            }
          }
          dxdt[i] = dx0;
        }
      }
    }
  };

} // namespace secsse
