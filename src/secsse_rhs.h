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

  template <typename T>
  using rvector = RcppParallel::RVector<T>;

  template <typename T>
  using rmatrix = RcppParallel::RMatrix<T>;

  
  template <typename T>
  class vector_view_t {
  public:
    vector_view_t(T* data, size_t n) : first_(data), n_(n) {};

    size_t size() const noexcept { return n_; }
    T* begin() noexcept { return first_; }
    T* end() noexcept { return first_ + n_; }
    T& operator[](size_t i) { return *(first_ + i); }
    void advance(size_t s) noexcept { first_ += s; }

  private:
    T* first_ = nullptr;
    size_t n_ = 0;
  };


  enum class OdeVariant {
    normal_tree,
    complete_tree,
    ct_condition
  };


  inline auto flat_q_matrix(const Rcpp::NumericMatrix& rq) {
    assert(rq.nrow() == rq.ncol());
    const auto d = static_cast<size_t>(rq.nrow());
    auto q = std::vector<double>(d * d);
    auto qv = vector_view_t<double>{q.data(), d};
    for (size_t i = 0; i < d; ++i, qv.advance(d)) {
      auto qrow = rq.row(i);
      for (size_t j = 0; j < d; ++j) {
        qv[j] = qrow[j];
      }
    }
    return q;
  }


  template <OdeVariant variant>
  class ode_standard {
    rvector<const double> l_;
    rvector<const double> m_;
    const std::vector<double> q_;

  public:
    ode_standard(const Rcpp::NumericVector& l,
                 const Rcpp::NumericVector& m,
                 const Rcpp::NumericMatrix& q)
    : l_(l), m_(m), q_(flat_q_matrix(q)) {
    }

    size_t size() const noexcept { return l_.size(); }

    void mergebranch(const std::vector<double>& N,
                     const std::vector<double>& M,
                     std::vector<double>& out) const {
      const auto d = size();
      assert(2 * d == out.size());
      for (size_t i = 0; i < d; ++i) {
        out[i] = M[i];
        out[i + d] = M[i + d] * N[i + d] * l_[i];
      }
    }

    void operator()(const std::vector<double>& x,
                    std::vector<double>& dxdt,   // NOLINT [runtime/references]
                    const double /* t */) const
    {
      const auto d = size();
      if  constexpr (variant == OdeVariant::normal_tree) {
        // normal tree
        auto qv = vector_view_t<const double>{q_.data(), d};
        for (size_t i = 0; i < d; ++i, qv.advance(d)) {
          const double t0 = l_[i] + m_[i];
          const double t1 = l_[i] * x[i];
          double dx0 = m_[i] + (t1 - t0) * x[i];
          double dxd = (2 * t1 - t0) * x[i + d];
          for (size_t j = 0; j < d; ++j) {
            dx0 += (x[j] - x[i]) * qv[j];
            dxd += (x[j + d] - x[i + d]) * qv[j];
          }
          dxdt[i] = dx0;
          dxdt[i + d] = dxd;
        }
      }
      else if constexpr ((variant == OdeVariant::complete_tree) || 
                         (variant == OdeVariant::ct_condition)) {
        // complete tree including extinct branches or conditioning
        auto qv = vector_view_t<const double>{q_.data(), d};
        for (size_t i = 0; i < d; ++i, qv.advance(d)) {
          double dx0 = (m_[i] - (l_[i] * x[i])) * (1 - x[i]);
          double dxd = -(l_[i] + m_[i]) * x[i + d];
          for (size_t j = 0; j < d; ++j) {
            dx0 += (x[j] - x[i]) * qv[j];
            dxd += (x[j + d] - x[i + d]) * qv[j];
          }
          dxdt[i] = dx0;
          dxdt[i + d] = dxd;
        }
      }
    }
  };

  struct ode_cla_precomp_t {
    std::vector<double> ll;               // flat, transposed ll matrices
    std::vector<std::vector<size_t>> nz;  // indices of non-zero values
    std::vector<double> lambda_sum;

    explicit ode_cla_precomp_t(const Rcpp::List& Rll) {
      const auto n = Rll.size();
      auto probe = Rcpp::as<Rcpp::NumericMatrix>(Rll[0]);
      assert(probe.nrow() == probe.ncol());
      const auto d = static_cast<size_t>(probe.nrow());
      ll.resize(n * d * d, 0.0);
      nz.resize(n * d, {});
      auto llv = vector_view_t<double>{ll.data(), d};
      auto nzv = nz.begin();
      for (int i = 0; i < Rll.size(); ++i) {
        // we all love deeply nested loops...
        rmatrix<const double> mr(Rcpp::as<Rcpp::NumericMatrix>(Rll[i]));
        auto& ls = lambda_sum.emplace_back(0.0);
        for (size_t j = 0; j < mr.nrow(); ++j, llv.advance(d), ++nzv) {
          for (size_t k = 0; k < d; ++k) {
            if (0.0 != (llv[k] = mr.row(j)[k])) {
              nzv->push_back(k);
              ls += llv[k];
            }
          }
        }
      }
    }
  };


  template <OdeVariant variant>
  class ode_cla {
    const rvector<const double> m_;
    const std::vector<double> q_;   // flat, transposed q matrix
    const ode_cla_precomp_t prec_;

  public:

    ode_cla(const Rcpp::List ll,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericMatrix& q)
    : m_(m), q_(flat_q_matrix(q)), prec_(ll) {
    }

    size_t size() const noexcept { return m_.size(); }

    void mergebranch(const std::vector<double>& N,
                     const std::vector<double>& M,
                     std::vector<double>& out) const {
      const auto d = size();
      assert(2 * d == out.size());
      auto llv = vector_view_t<const double>(prec_.ll.data(), d);
      for (size_t i = 0; i < d; ++i) {
        out[i] = M[i];
        out[i + d] = 0.0;
        for (size_t j = 0; j < d; ++j, llv.advance(d)) {
          for (size_t k = 0; k < d; ++k) {
            out[i + d] += llv[k] * (N[j + d] * M[k + d] + M[j + d] * N[k + d]);
          }
        }
        out[i + d] *= 0.5;
      }
    }

    void operator()(const std::vector<double>& x,
                    std::vector<double>& dxdt,
                    const double /* t */ ) const
    {
      
      const auto d = size();
      if constexpr (variant == OdeVariant::normal_tree) {
        auto llv = vector_view_t<const double>(prec_.ll.data(), d);
        auto nzv = prec_.nz.begin();
        auto qv = vector_view_t<const double>{q_.data(), d};
        for (size_t i = 0; i < d; ++i, qv.advance(d)) {
          double dx0 = 0.0;
          double dxd = 0.0;
          for (size_t j = 0; j < d; ++j, llv.advance(d), ++nzv) {
            for (auto k : *nzv) {
              dx0 += llv[k] * (x[j] * x[k]);
              dxd += llv[k] * (x[j] * x[k + d] + x[j + d] * x[k]);
            }
            dx0 += (x[j] - x[i]) * qv[j];
            dxd += (x[j + d] - x[i + d]) * qv[j];
          }
          dxdt[i] = dx0 + m_[i] - (prec_.lambda_sum[i] + m_[i]) * x[i];
          dxdt[i + d] = dxd - (prec_.lambda_sum[i] + m_[i]) * x[i + d];
        }
      }
      else if constexpr (variant == OdeVariant::complete_tree) {
        // complete tree including extinct branches
        auto qv = vector_view_t<const double>{q_.data(), d};
        for (size_t i = 0; i < d; ++i, qv.advance(d)) {
          double dxd = -(prec_.lambda_sum[i] + m_[i]) * x[i + d];
          for (size_t j = 0; j < d; ++j) {
            dxd += (x[j + d] - x[i + d]) * qv[j];
          }
          dxdt[i + d] = dxd;
        }
      }
      else if constexpr (variant == OdeVariant::ct_condition) {
        auto llv = vector_view_t<const double>(prec_.ll.data(), d);
        auto nzv = prec_.nz.begin();
        auto qv = vector_view_t<const double>{q_.data(), d};
        for (size_t i = 0; i < d; ++i, qv.advance(d)) {
          double dx0 = m_[i] * (1 - x[i]);
          for (size_t j = 0; j < d; ++j, llv.advance(d), ++nzv) {
            dx0 += (x[j] - x[i]) * qv[j];
            for (auto k : *nzv) {
              dx0 += llv[k] * (x[j] * x[k] - x[i]);
            }
          }
          dxdt[i] = dx0;
        }
      }
    }
  };
} // namespace secsse
