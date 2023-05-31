//
//  Copyright (c) 2021 - 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#include "Rcpp.h"                     // NOLINT [build/include_subdir]
#include "util.h"                     // NOLINT [build/include_subdir]


class ode_standard {
public:
  ode_standard(const std::vector<double>& l,
               const std::vector<double>& m,
               const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q) {
    d = l.size();
  }
  
  ode_standard(const Rcpp::NumericVector& l,
               const Rcpp::NumericVector& m,
               const Rcpp::NumericMatrix& q) {
    l_ = std::vector<double>(l.begin(), l.end());
    m_ = std::vector<double>(m.begin(), m.end());
    numericmatrix_to_vector(q, &q_);
    d = l_.size();
  }
  
  void operator()(const std::vector< double > &x,
                std::vector<  double > &dxdt,   // NOLINT [runtime/references]
                const double /* t */) {
    for (size_t i = 0; i < d; ++i) {
      if (l_[i] != 0.0) {
        dxdt[i] = m_[i] - (l_[i] + m_[i]) * x[i] +
          l_[i] * x[i] * x[i];
        long double FF3 = -1.0 * l_[i] - m_[i] + 2 * l_[i] * x[i];
        dxdt[i + d] = FF3 * x[i + d];
      } else {
        dxdt[i] = - 1.0 * m_[i] * x[i] + m_[i];
        dxdt[i + d] = -1.0 * m_[i] * x[i + d];
      }
      
      for (size_t j = 0; j < d; ++j) {
        long double diff_e = x[j] - x[i];
        dxdt[i] += diff_e * q_[i][j];
        
        long double diff_d = x[j + d] - x[i + d];
        dxdt[i + d] += diff_d * q_[i][j];
      }
    }
    return;
  }
  
  double get_l(int index) const {
    return l_[index];
  }
  
  size_t get_d() const {
    return d;
  }
  
private:
  std::vector< double > l_;
  std::vector< double > m_;
  std::vector< std::vector< double >> q_;
  size_t d;
};

class ode_standard_ct {
public:
  ode_standard_ct(const std::vector<double>& l,
                  const std::vector<double>& m,
                  const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q) {
    d = l.size();
  }
  
  ode_standard_ct(const Rcpp::NumericVector& l,
                  const Rcpp::NumericVector& m,
                  const Rcpp::NumericMatrix& q) {
    l_ = std::vector<double>(l.begin(), l.end());
    m_ = std::vector<double>(m.begin(), m.end());
    numericmatrix_to_vector(q, &q_);
    d = l_.size();
  }
  
  void operator()(const std::vector< double > &x ,
                std::vector<  double > &dxdt,   // NOLINT [runtime/references]
                const double /* t */) {
    for (int i = 0; i < d; ++i) {
      long double diff_1 = (m_[i] - (l_[i] * x[i]));
      dxdt[i] =  diff_1 * (1 - x[i]);
      dxdt[i + d] = -1.0 * (l_[i] + m_[i]) * x[i + d];
    }
    
    for (int j = 0; j < d; ++j) {
      for (int k = 0; k < d; ++k) {
        long double diff_e = x[k] - x[j];
        dxdt[j] +=  q_[j][k] * diff_e;
        
        long double diff_d = x[k + d] - x[j + d];
        dxdt[j + d] += q_[j][k] * diff_d;
      }
    }
    
    return;
  }
  
  double get_l(int index) const {
    return l_[index];
  }
  
  size_t get_d() const {
    return d;
  }
  
private:
  std::vector< double > l_;
  std::vector< double > m_;
  std::vector< std::vector< double >> q_;
  size_t d;
};

class ode_cla {
  // used for normal tree
public:
  ode_cla(const std::vector<std::vector<std::vector<double>>>& l,
          const std::vector<double>& m,
          const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) {
    lambda_sum = std::vector<long double>(d, 0.0);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          lambda_sum[i] += l_[i][j][k];
        }
      }
    }
  }
  
  void operator()(const std::vector< double > &x ,
                std::vector< double > &dxdt,    // NOLINT [runtime/references]
                const double /* t */) const {
    for (int i = 0; i < d; ++i) {
      double Df = 0.0;
      double Ef = 0.0;
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0.0) {    // slightly safer.
            long double ff1 = (x[j] * x[k + d] + x[j + d] * x[k]);
            long double ff2 = (x[j] * x[k]);
            
            Df += l_[i][j][k] * ff1;
            Ef += l_[i][j][k] * ff2;
          }
        }
      }
      
      dxdt[i]     = Ef + m_[i] - (lambda_sum[i] + m_[i]) * x[i];
      dxdt[i + d] = Df + (-lambda_sum[i] - m_[i]) * x[i + d];
      
      for (size_t j = 0; j < d; ++j) {
        // q_[i][j] is always non-zero.
        long double temp1 = (x[j]     - x[i]);
        dxdt[i]     += q_[i][j] * temp1;
        long double temp2 = (x[j + d] - x[i + d]);
        dxdt[i + d] += q_[i][j] * temp2;
      }
    }
    return;
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  }
  
  size_t get_d() const {
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
  std::vector< long double > lambda_sum;
};

class ode_cla_backup {
public:
  ode_cla_backup(const std::vector<std::vector<std::vector<double>>>& l,
                 const std::vector<double>& m,
                 const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) {
    lambda_sum = std::vector<long double>(d, 0.0);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          lambda_sum[i] += l_[i][j][k];
        }
      }
    }
  }
  
  void operator()(const std::vector< double > &x,
                std::vector< double > &dxdt,   // NOLINT [runtime/references]
                const double /* t */) const {
    for (int i = 0; i < d; ++i) {
      long double lamEE = 0.0;
      long double lamDE = 0.0;
      
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0) {
            long double FF1 = x[j] * x[k];
            lamEE += l_[i][j][k] * FF1;
            
            long double FF3 = x[d + j] * x[k];
            long double FF2 = x[d + k] * x[k];
            lamDE += l_[i][j][k] * (FF3 + FF2);
          }
        }
      }
      
      long double FF1 = m_[i] - (lambda_sum[i] + m_[i]) * x[i];
      dxdt[i] = FF1 + lamEE;
      
      long double FF2 = (-lambda_sum[i] - m_[i]) * x[i + d];
      dxdt[i + d] = FF2 + lamDE;
      
      for (int j = 0; j < d; ++j) {
        long double diff = x[j] - x[i];
        dxdt[i] += diff * q_[i][j];
        
        long double diff2 = x[j + d] - x[i + d];
        dxdt[i + d] += diff2 * q_[i][j];
      }
    }
    return;
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  }
  
  size_t get_d() const {
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
  std::vector< long double > lambda_sum;
};

class ode_cla_d {
  // used for complete tree including extinct branches
public:
  ode_cla_d(const std::vector<std::vector<std::vector<double>>>& l,
            const std::vector<double>& m,
            const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) {
    lambda_sum = std::vector<long double>(d, 0.0);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          lambda_sum[i] += l_[i][j][k];
        }
      }
    }
  }
  
  void single_step(const std::vector< double > &x ,
                   std::vector< double > &dxdt) {  // NOLINT [runtime/references]
    for (int i = 0; i < d; ++i) {
      dxdt[i + d] = -1.0 * (lambda_sum[i] + m_[i]) * x[i + d];
      for (int j = 0; j < d; ++j) {
        long double dx = x[j + d] - x[i + d];
        dxdt[i + d] +=  q_[i][j] * dx;
      }
    }
  }
  void operator()(const std::vector< double > &x ,
                std::vector< double > &dxdt,    // NOLINT [runtime/references]
                const double /* t */) const {
    for (int i = 0; i < d; ++i) {
      dxdt[i + d] = -1.0 * (lambda_sum[i] + m_[i]) * x[i + d];
      for (int j = 0; j < d; ++j) {
        long double dx = x[j + d] - x[i + d];
        dxdt[i + d] +=  q_[i][j] * dx;
      }
    }
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  }
  
  size_t get_d() const {
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
  std::vector<long double> lambda_sum;
};

class ode_cla_e {
  // used for ct conditioning.
public:
  ode_cla_e(const std::vector<std::vector<std::vector<double>>>& l,
            const std::vector<double>& m,
            const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) {
  }
  
  void operator()(const std::vector< double > &x ,
                std::vector< double > &dxdt, // NOLINT [runtime/references]
                const double /* t */) const {
    for (int i = 0; i < d; ++i) {
      dxdt[i] = 0.0;
      if (m_[i] != 0.0) {
        dxdt[i] = m_[i] * (1.0 - x[i]);
      }
      for (int j = 0; j < d; ++j) {
        long double diff = (x[j] - x[i]);
        dxdt[i] += q_[i][j] * diff;
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0.0) {
            long double diff2 = (x[j] * x[k] - x[i]);
            dxdt[i] += l_[i][j][k] * diff2;
          }
        }
      }
    }
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  }
  
  size_t get_d() const {
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
};

//////// STORAGE section - these are used for plotting
//////// these versions also store intermediate results!

class ode_standard_store {
public:
  ode_standard_store(const std::vector<double>& l,
                     const std::vector<double>& m,
                     const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q) {
    d = l.size();
  }
  
  ode_standard_store(const Rcpp::NumericVector& l,
                     const Rcpp::NumericVector& m,
                     const Rcpp::NumericMatrix& q) {
    l_ = std::vector<double>(l.begin(), l.end());
    m_ = std::vector<double>(m.begin(), m.end());
    numericmatrix_to_vector(q, &q_);
    d = l_.size();
  }
  
  void operator()(const std::vector< double > &x ,
                std::vector<  double > &dxdt,  // NOLINT [runtime/references]
                const double t) {
    for (size_t i = 0; i < d; ++i) {
      if (l_[i] != 0.0) {
        dxdt[i] = m_[i] - (l_[i] + m_[i]) * x[i]  +
          l_[i] * x[i] * x[i];
        long double FF3 = -1.0 * l_[i] - m_[i] + 2 * l_[i] * x[i];
        dxdt[i + d] = FF3 * x[ i + d];
      } else {
        dxdt[i] = - 1.0 * m_[i] * x[i] + m_[i];
        dxdt[i + d] = -1.0 * m_[i] * x[i + d];
      }
      
      for (size_t j = 0; j < d; ++j) {
        long double diff_e = x[j] - x[i];
        dxdt[i] += diff_e * q_[i][j];
        
        long double diff_d = x[j + d] - x[i + d];
        dxdt[i + d] += diff_d * q_[i][j];
      }
    }
    
    stored_t.push_back(t);
    stored_states.push_back(x);
    return;
  }
  
  double get_l(int index) const {
    return l_[index];
  }
  
  size_t get_d() const {
    return d;
  }
  
  std::vector< std::vector<double >> get_stored_states() {
    return stored_states;
  }
  
  std::vector<double> get_stored_t() {
    return stored_t;
  }
  
private:
  std::vector< double > l_;
  std::vector< double > m_;
  std::vector< std::vector< double >> q_;
  std::vector< std::vector<double >> stored_states;
  std::vector<double> stored_t;
  size_t d;
};

class ode_cla_store {
  // used for normal tree
public:
  ode_cla_store(const std::vector<std::vector<std::vector<double>>>& l,
                const std::vector<double>& m,
                const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) {
    lambda_sum = std::vector<long double>(d, 0.0);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          lambda_sum[i] += l_[i][j][k];
        }
      }
    }
  }
  
  void operator()(const std::vector< double > &x ,
                std::vector< double > &dxdt,    // NOLINT [runtime/references]
                const double t /* t */ )  {
    stored_t.push_back(t);
    stored_states.push_back(x);
    
    for (int i = 0; i < d; ++i) {
      double Df = 0.0;
      double Ef = 0.0;
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0.0) {   // slightly safer.
            long double ff1 = (x[j] * x[k + d] + x[j + d] * x[k]);
            long double ff2 = (x[j] * x[k]);
            
            Df += l_[i][j][k] * ff1;
            Ef += l_[i][j][k] * ff2;
          }
        }
      }
      
      dxdt[i]     = Ef + m_[i] - (lambda_sum[i] + m_[i]) * x[i];
      dxdt[i + d] = Df + (-lambda_sum[i] - m_[i]) * x[i + d];
      
      for (size_t j = 0; j < d; ++j) {
        // q_[i][j] is always non-zero.
        long double temp1 = (x[j]     - x[i]);
        dxdt[i]     += q_[i][j] * temp1;
        
        long double temp2 = (x[j + d] - x[i + d]);
        dxdt[i + d] += q_[i][j] * temp2;
      }
    }
    return;
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  }
  
  size_t get_d() const {
    return d;
  }
  
  std::vector< std::vector<double >> get_stored_states() const {
    return stored_states;
  }
  
  std::vector<double> get_stored_t() const {
    return stored_t;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
  std::vector< long double > lambda_sum;
  std::vector< std::vector<double >> stored_states;
  std::vector<double> stored_t;
};