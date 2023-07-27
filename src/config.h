//
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#ifndef SRC_CONFIG_H_
#define SRC_CONFIG_H_

// Special case to make use of some steppers that would include
// boost/get_pointer.hpp
#ifndef BOOST_NO_AUTO_PTR
# define BOOST_NO_AUTO_PTR
#endif

// Addresses unitialized member variable bulirsch_stoer<>::m_dt_last.
//
// The issue is *not* fixed in BOOST_VERSION 1.81.1.
// We need to check for fixes in upcomming boost (BH) releases.
//
// Uncomment if unitialized member variable bulirsch_stoer::m_dt_last
// is fixed in boost (BH):
#define USE_BULRISCH_STOER_PATCH

// Default initial dt factor for interation stepper.
// The initial dt is calculated as SECSEE_DEFAULT_DTF * (t1 - t0).
// All used steppers are adaptive, thus the value shouldn't really matter
#define SECSSE_DEFAULT_DTF 0.01

// Default initial dt factor for interation stepper in iterative 'store' mode.
// The initial dt is calculated as SECSEE_DEFAULT_EVAL_DTF * (t1 - t0).
// All used steppers are adaptive, thus the value shouldn't really matter
#define SECSSE_DEFAULT_EVAL_DTF 0.1

// Uncomment to enable nested parallelism.
// This feature may improve or may deterioate performance.
// Default is disabled.
//#define SECSSE_NESTED_PARALLELISM

#endif  // SRC_CONFIG_H_
