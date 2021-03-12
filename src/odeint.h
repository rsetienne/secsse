//
//  odeint.h
//  calc_ll_secsse
//
//  Created by Thijs Janzen on 08/03/2021.
//  Copyright © 2021 Thijs Janzen. All rights reserved.
//

#ifndef odeint_h
#define odeint_h

// [[Rcpp::depends(BH)]]
#include "boost/numeric/odeint.hpp"
#include "util.h"

#include "Rcpp.h"


class MyOde{
 public:

  MyOde(const std::vector<double>& l,
        const std::vector<double>& m,
        const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q) {
    d = l.size();
  }

   MyOde(const Rcpp::NumericVector& l,
         const Rcpp::NumericVector& m,
         const Rcpp::NumericMatrix& q) {

     l_ = std::vector<double>(l.begin(), l.end());
     m_ = std::vector<double>(m.begin(), m.end());
     numericmatrix_to_vector(q, q_);

     d = l_.size();
   }


  void operator()( const std::vector< double > &x ,
                         std::vector< double > &dxdt,
                   const double /* t */ ) {
    for (int i = 0; i < d; ++i) {
        double q_mult_one = 0.0;
        double q_ds = 0.0;
        double q_es = 0.0;
        for (int j = 0; j < d; ++j) {
          q_mult_one += q_[i][j];
          q_ds += q_[i][j] * x[j + d];
          q_es += q_[i][j] * x[j];
        }

        dxdt[i] = m_[i] -
                     (l_[i] + m_[i] + q_mult_one) * x[i] +
                     l_[i] * x[i] * x[i] +
                     q_es;

        dxdt[i + d] = (l_[i] + m_[i] + q_mult_one) * x[i + d] * -1 +
                     l_[i] * x[i] * x[i + d] * 2 +
                     q_ds;
    }
    return;
  }

private:
  std::vector< double > l_;
  std::vector< double > m_;
  std::vector< std::vector< double >> q_;
  size_t d;
};

class MyOde_cla{
public:

  MyOde_cla(const Rcpp::List& l,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericMatrix& q) : d(m.size()){

    std::vector< std::vector< double >> entry;
    for (int i = 0; i < l.size(); ++i) {
      Rcpp::NumericMatrix temp = l[i];
      numericmatrix_to_vector(temp, entry);
      l_.push_back(entry);
    }

    m_ = std::vector<double>(m.begin(), m.end());
    numericmatrix_to_vector(q, q_);
  }


  void operator()( const std::vector< double > &x ,
                std::vector< double > &dxdt,
                const double /* t */ ) const {

    assert(x.size() == 2 * d);
    assert(dxdt.size() == 2 * d);
    assert(q_.size() == d);
    assert(q_[0].size() == d);

    for (int i = 0; i < d; ++i) {
      double lambda_sum = 0.0;
      double Q_one = 0.0;
      double Q_Ds = 0.0;
      double Q_Es = 0.0;
      double Df = 0.0;
      double Ef = 0.0;
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] == 0.0)
            continue; // all entries below multiply with 0!

          lambda_sum += l_[i][j][k];
          Df += l_[i][j][k] * (x[j] * x[k + d] + x[j + d] * x[k]);
          Ef += l_[i][j][k] * (x[j] * x[k]);
        }

        Q_one += q_[i][j];
        Q_Es  += q_[i][j] * x[j];
        Q_Ds  += q_[i][j] * x[j + d];
      }

      // dE
      dxdt[i] =     -1 * x[i] *     (lambda_sum + m_[i] + Q_one) + Q_Es + m_[i] + Ef;
      // dD
      dxdt[i + d] = -1 * x[i + d] * (lambda_sum + m_[i] + Q_one) + Q_Ds + Df;
    } // this is correct
    return;
  }

private:
   std::vector< std::vector< std::vector< double > > > l_;
   std::vector< double > m_;
   std::vector< std::vector< double >> q_;
   size_t d;
};

namespace bno = boost::numeric::odeint;


#endif /* odeint_h */
