// Copyright 2021 - 2023 Hanno Hildenbrandt
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
#pragma once

#ifdef USE_BULRISCH_STOER_PATCH

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

using bstime_t = boost::units::quantity<boost::units::si::dimensionless, double>;  // NOLINT [whitespace/line_length]

#else  // USE_BULRISCH_STOER_PATCH

// The default. Causes unitialized member m_last_dt in
// boost::odeint::bulrisch_stoer<>, declared in
// boost/numreic/odeint/stepper/bulrisch_stoer.hpp
using bstime_t = double;

#endif  // USE_BULRISCH_STOER_PATCH

#include <iostream>
#include <utility>   // std::move
#include <memory>    // std::unique_ptr
#include <string>
#include <vector>

// [[Rcpp::depends(BH)]]
#include "config.h"                   // NOLINT [build/include_subdir]
#include "boost/numeric/odeint.hpp"   // NOLINT [build/include_subdir]


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
void integrate(STEPPER stepper, ODE ode, STATE* y,
               double t0, double t1, double dt) {
  bno::integrate_adaptive(stepper, ode, (*y), t0, t1, dt);
}

template <
  typename STATE,
  typename ODE
>
void integrate(const std::string& stepper_name,
               std::unique_ptr<ODE> ode,
               STATE* y,
               double t0, double t1,
               double dt, double atol, double rtol) {
  if ("odeint::runge_kutta_cash_karp54" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_cash_karp54<STATE>>(atol,
                                                                        rtol),
              std::ref(*ode), y, t0, t1, dt);
  } else if ("odeint::runge_kutta_fehlberg78" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_fehlberg78<STATE>>(
        atol, rtol), std::ref(*ode), y, t0, t1, dt);
  } else if ("odeint::runge_kutta_dopri5" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_dopri5<STATE>>(atol, rtol),
              std::ref(*ode), y, t0, t1, dt);
  } else if ("odeint::bulirsch_stoer" == stepper_name) {
    integrate(bno::bulirsch_stoer<STATE>(atol, rtol),
              std::ref(*ode), y, bstime_t{t0}, bstime_t{t1}, dt);
  } else if ("odeint::runge_kutta4" == stepper_name) {
    integrate(bno::runge_kutta4<STATE>(), std::ref(*ode), y, t0, t1, dt);
  } else {
    throw std::runtime_error("odeintcpp::integrate: unknown stepper");
  }
}

template <
  typename STATE,
  typename ODE
>
void integrate_full(const std::string& stepper_name,
                    std::unique_ptr<ODE> ode,
                    STATE* y,
                    double t0, double t1, double dt,
                    double atol, double rtol,
                    std::vector< std::vector<double>>* yvals,
                    std::vector<double>* tvals) {
  if ("odeint::runge_kutta_cash_karp54" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_cash_karp54<STATE>>(atol,
                                                                        rtol),
              std::ref(*ode), y, t0, t1, dt);
  } else if ("odeint::runge_kutta_fehlberg78" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_fehlberg78<STATE>>(atol,
                                                                       rtol),
                                                std::ref(*ode), y, t0, t1, dt);
  } else if ("odeint::runge_kutta_dopri5" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_dopri5<STATE>>(atol, rtol),
              std::ref(*ode), y, t0, t1, dt);
  } else if ("odeint::bulirsch_stoer" == stepper_name) {
    integrate(bno::bulirsch_stoer<STATE>(atol, rtol), std::ref(*ode), y,
              t0, t1, dt);
  } else if ("odeint::runge_kutta4" == stepper_name) {
    integrate(bno::runge_kutta4<STATE>(), std::ref(*ode), y, t0, t1, dt);
  } else {
    throw std::runtime_error("odeintcpp::integrate: unknown stepper");
  }

  (*yvals) = (*ode).get_stored_states();
  (*tvals) = (*ode).get_stored_t();
  return;
}

template <
  typename STATE,
  typename ODE
>
void integrate(const std::string& stepper_name,
               std::unique_ptr<ODE> ode, STATE* y, double t0, double t1) {
  integrate(stepper_name, std::move(ode), (*y), t0, t1,
            (t1 - t0) / default_init_steps, default_atol, default_rtol);
}

template <
  typename STATE,
  typename ODE
>
void integrate(std::unique_ptr<ODE> ode, STATE* y,
               double t0, double t1, double dt, double atol, double rtol) {
  integrate(default_stepper_name, std::move(ode), (*y), t0, t1, dt, atol, rtol);
}


template <
  typename STATE,
  typename ODE
>
void integrate(std::unique_ptr<ODE> ode, STATE* y, double t0, double t1) {
  integrate(default_stepper_name, std::move(ode), (*y), t0, t1);
}

}   // namespace odeintcpp
