//
//  odeint.h
//  calc_ll_secsse
//
//  Created by Thijs Janzen on 08/03/2021.
//  Copyright © 2021 Thijs Janzen. All rights reserved.
//

#ifndef odeint_h
#define odeint_h

#include <iostream>

// [[Rcpp::depends(BH)]]
#include "Rcpp.h"
#include "boost/numeric/odeint.hpp"
// #include <boost/multiprecision/mpfr.hpp>
// #include <quadmath.h>
// #include <boost/multiprecision/float128.hpp>
// #include "boost/multiprecision/cpp_bin_float.hpp"

#include <boost/multiprecision/mpfr.hpp>

namespace mp = boost::multiprecision;

// using high_prec_double = mp::cpp_bin_float_quad;   //mp::mpf_float_100;

using high_prec_double = mp::mpf_float_100;

#include "util.h"



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
    numericmatrix_to_vector(q, q_);
    
    d = l_.size();
  }
  
  void operator()( const std::vector< double > &x ,
                std::vector<  double > &dxdt,
                const double ) {
    for (int i = 0; i < d; ++i) {
      double q_mult_one = 0.0;
      double q_ds = 0.0;
      double q_es = 0.0;
      for (int j = 0; j < d; ++j) {
        q_mult_one += q_[i][j];
        q_ds += q_[i][j] * x[j + d];
        q_es += q_[i][j] * x[j];
      }
      
      //dxdt[i] = m_[i] -
      //  (l_[i] + m_[i] + q_mult_one) * x[i] +
      //  l_[i] * x[i] * x[i] +
      //  q_es;
      
      dxdt[i] = m_[i] - x[i] * (m_[i] + q_mult_one + l_[i] * (x[i] + 1)) + q_es;
        
      dxdt[i + d] =  x[i + d] * (l_[i] * (x[i] * 2 - 1) - m_[i] - q_mult_one) + q_ds;
      
    }
    return;
  }
  
  /*
  
  void operator()( const std::vector< double > &x ,
                std::vector<  double > &dxdt,
                const double ) {
    for (size_t i = 0; i < d; ++i) {
      high_prec_double q_mult_one(0.0);
      high_prec_double q_ds(0.0);
      high_prec_double q_es(0.0);
      for (size_t j = 0; j < d; ++j) {
        q_mult_one += q_[i][j];
        q_ds += q_[i][j] * x[j + d];
        q_es += q_[i][j] * x[j];
      }
     
     high_prec_double mult1  = l_[i] * (x[i] + 1);
     high_prec_double sum    = m_[i] + q_mult_one   + mult1;
     high_prec_double Edt = m_[i] - (x[i] * sum) + q_es;
     
     dxdt[i] = Edt.convert_to<double>();
     
     // dxdt[i + d] = (l_[i] + m_[i] + q_mult_one) * x[i + d] * -1 +
    //    l_[i] * x[i] * x[i + d] * 2 +
    //    q_ds;
    //    
     high_prec_double subtr = l_[i] * (x[i] * 2 - 1) - m_[i] - q_mult_one;        // l_[i] * (x[i] * 2) - (l_[i] + m_[i] + q_mult_one);
     high_prec_double mult  = x[i + d] * subtr;
     high_prec_double Ds    = mult + q_ds;
     dxdt[i + d] = Ds.convert_to<double>();   
    }
    return;
  }
   */
  /*
     void operator()( const std::vector< double > &x ,
                   std::vector<  double > &dxdt,
                   const double  ) {
       for (size_t i = 0; i < d; ++i) {
        // dxdt[i] = m_[i] - (l_[i] + m_[i]) * x[i] +  l_[i] * x[i] * x[i];
         
         if (l_[i] != 0.0) {
         
           dxdt[i] = m_[i] - x[i] * (l_[i] * (m_[i] + x[i]));
           
          //dxdt[i + d] = (-1 * (l_[i] + m_[i]) + 2 * l_[i] * x[i]) * x[i + d];
          
           // dxdt[i + d] =  x[i + d] * (2 * l_[i] * x[i] - (l_[i] + m_[i]))
           dxdt[i + d] =  x[i + d] * (l_[i] * (2 * x[i] - 1) - m_[i]);
         } else {
           dxdt[i] = m_[i];
           dxdt[i + d] = -1 * x[i + d] * m_[i];
        }
         
         
         for (int j = 0; j < d; ++j) {
           
           long double dd = (x[j] - x[i]);
           
           dxdt[i] += q_[i][j] * dd;
           
           long double dd2 = (x[j + d] - x[i + d]);
           dxdt[i + d] += q_[i][j] * dd2;
         }
        
       }
       return;
     }
     */    
     
     
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


class ode_cla{
public:
  
  ode_cla(const std::vector<std::vector<std::vector<double>>>& l,
          const std::vector<double>& m,
          const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) {
 //   std::cerr << "ode_cla made\n";
    c_e = 0.0;
    c_d = 0.0;
  }
  
  void operator_base(const std::vector< double > &x ,
                std::vector< double > &dxdt,
                const double /* t */ ) const {
    
    for (int i = 0; i < d; ++i) {
      double lambda_sum(0.0);
      double Df(0.0);
      double Ef(0.0);
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0.0) { // slightly safer.
            lambda_sum += l_[i][j][k];
            Df +=         l_[i][j][k] * (x[j] * x[k + d] + x[j + d] * x[k]);
            Ef +=         l_[i][j][k] * (x[j] * x[k]);
          }
        }
      }
      // dxdt[i] = Ef - x[i] * lambda_sum + m_[i] * (x[i] + 1); // remaining: Q_es - x[i] * Q_one
      dxdt[i]     = Ef + m_[i] + x[i] * (m_[i] - lambda_sum);
      dxdt[i + d] = Df - x[i + d] * (lambda_sum + m_[i]); //remaining: Q_ds - x[i+d] * Q_one
      
      for (size_t j = 0; j < d; ++j) {
        // q_[i][j] is always non-zero.
        long double t1 = (x[j]     - x[i]);
        dxdt[i]     += q_[i][j] * t1;
        
        long double t2 = (x[j + d] - x[i + d]);
        dxdt[i + d] += q_[i][j] * t2;
      }
    }
    return;
  }
  
  void kahan_sum(double& sum, double& c, const double& add) {
    t = sum + add;
    c = (t - sum) - add;
    sum = t;
  //  std::cerr << sum << " " << c << " " << add << "\n";
  }
  
  void operator()(const std::vector< double > &x ,
                      std::vector< double > &dxdt,
                      const double /* t */ ) {
    
    for (int i = 0; i < d; ++i) {
      double lambda_sum = 0.0;
      double Df = 0.0;
      double Ef = 0.0;
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0.0) { // slightly safer.
            lambda_sum += l_[i][j][k];
            Df +=         l_[i][j][k] * (x[j] * x[k + d] + x[j + d] * x[k]);
            Ef +=         l_[i][j][k] * (x[j] * x[k]);
          }
        }
      }
      // dxdt[i] = Ef - x[i] * lambda_sum + m_[i] * (x[i] + 1); // remaining: Q_es - x[i] * Q_one
      dxdt[i]     = Ef + m_[i] + x[i] * (m_[i] - lambda_sum);
      dxdt[i + d] = Df - x[i + d] * (lambda_sum + m_[i]); //remaining: Q_ds - x[i+d] * Q_one
      
      
      for (size_t j = 0; j < d; ++j) {
        kahan_sum(dxdt[i],     c_e, q_[i][j] * (x[j]     - x[i]));
        kahan_sum(dxdt[i + d], c_d, q_[i][j] * (x[j + d] - x[i + d]));
      }
    }
    return;
  }
  
  // code below can be renamed to 'operator()' to function
  // but is incredibly slow.
  void kahan_operator(const std::vector< double > &x ,
                    std::vector< double > &dxdt,
                    const double /* t */ ) {
    
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0; // correction factor for kahan-sum
  
    for (int i = 0; i < d; ++i) {
      double lambda_sum(0.0);
      double Df(0.0);
      double Ef(0.0);
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0.0) { // slightly safer.
            
            kahan_sum(lambda_sum, c1, l_[i][j][k]);
            kahan_sum(Df, c2, l_[i][j][k] * (x[j] * x[k + d] + x[j + d] * x[k]));
            kahan_sum(Ef, c3, l_[i][j][k] * (x[j] * x[k]));
          }
        }
      }
      kahan_sum(dxdt[i], c_e, Ef - x[i] * lambda_sum + m_[i] * (x[i] + 1));
      kahan_sum(dxdt[i + d], c_d, Df - x[i + d] * (lambda_sum + m_[i]));
      
      for (size_t j = 0; j < d; ++j) {
        kahan_sum(dxdt[i], c_e, q_[i][j] * (x[j] - x[i]));
        kahan_sum(dxdt[i + d], c_d, q_[i][j] * (x[j + d] - x[i + d]));
      }
    }
    return;
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    
    std::cerr << i << " " << j << " " << k << " " << l_[i][j][k] << "\n" <<  std::flush;
    
    return l_[i][j][k];
  } 
  
 /* size_t get_ll_size() const {
    return l_.size();
 }*/
  
  size_t get_d() const {
 //   std::cerr << "d: " << d << "\n";
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
  // kahan sum parameters:
  double c_e;
  double c_d;
  double t;
};

namespace bno = boost::numeric::odeint;

namespace odeintcpp {

namespace bno = boost::numeric::odeint;

constexpr char default_stepper_name[] = "odeint::runge_kutta_fehlberg78";
constexpr double default_atol = 1E-10;
constexpr double default_rtol = 1E-10;
constexpr double default_init_steps = 10;

template <
  typename STEPPER,
  typename STATE,
  typename ODE
>
void integrate(STEPPER stepper, ODE ode, STATE& y, double t0, double t1, double dt)
{
  bno::integrate_adaptive(stepper, ode, y, t0, t1, dt);
}


template <
  typename STATE,
  typename ODE
>
void integrate(const std::string& stepper_name, std::unique_ptr<ODE> ode, STATE& y, double t0, double t1, double dt, double atol, double rtol)
{
  if ("odeint::runge_kutta_cash_karp54" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_cash_karp54<STATE>>(atol, rtol),
              std::ref(*ode), y, t0, t1, dt);
  }
  else if ("odeint::runge_kutta_fehlberg78" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_fehlberg78<STATE>>(atol, rtol), std::ref(*ode), y, t0, t1, dt);
  }
  else if ("odeint::runge_kutta_dopri5" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_dopri5<STATE>>(atol, rtol), std::ref(*ode), y, t0, t1, dt);
  }
  else if ("odeint::bulirsch_stoer" == stepper_name) {
    integrate(bno::bulirsch_stoer<STATE>(atol, rtol), std::ref(*ode), y, t0, t1, dt);
  }
  else if ("odeint::runge_kutta4" == stepper_name) {
    integrate(bno::runge_kutta4<STATE>(), std::ref(*ode), y, t0, t1, dt);
  }
 
  else {
    throw std::runtime_error("odeintcpp::integrate: unknown stepper");
  }
}


template <
  typename STATE,
  typename ODE
>
void integrate(const std::string& stepper_name, std::unique_ptr<ODE> ode, STATE& y, double t0, double t1)
{
  integrate(stepper_name, std::move(ode), y, t0, t1, (t1 - t0) / default_init_steps, default_atol, default_rtol);
}


template <
  typename STATE,
  typename ODE
>
void integrate(std::unique_ptr<ODE> ode, STATE& y, double t0, double t1, double dt, double atol, double rtol)
{
  integrate(default_stepper_name, std::move(ode), y, t0, t1, dt, atol, rtol);
}


template <
  typename STATE,
  typename ODE
>
void integrate(std::unique_ptr<ODE> ode, STATE& y, double t0, double t1)
{
  integrate(default_stepper_name, std::move(ode), y, t0, t1);
}
}

#endif /* odeint_h */
