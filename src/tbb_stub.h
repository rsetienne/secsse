//  Copyright (c) 2026, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cstdlib>
#include <Rcpp.h>
#include <RcppParallel.h>   // pull RCPP_PARALLEL_USE_TBB


#if RCPP_PARALLEL_USE_TBB == 0

#include <algorithm>

namespace tbb {

  namespace task_arena {

    constexpr size_t automatic = size_t(-1);

  } // namespace task_arena
  
  
  class global_control {
  public:
    enum parameter {
        max_allowed_parallelism,
        thread_stack_size,
        terminate_on_exception
    };

    global_control(parameter /*p*/, size_t /*value*/) {}
    ~global_control() {};
    static size_t active_value(parameter param);  // undefined
  };

  
  template<typename InputIterator, typename Body>
  inline void parallel_for_each( InputIterator first, InputIterator last, Body&& body ) {
    std::for_each(first, last, std::forward<Body>(body));
  }

  
  template<typename Index, typename Func>
  inline void parallel_for(Index first, Index last, const Func f) {
    for (; first != last; ++first) {
      f(first);
    }
  }
        

  template<typename Index, typename Func>
  inline void parallel_for(Index first, Index last, Index step, const Func f) {
    for (; first != last; first += step) {
      f(first);
    }
  }
        

} // namespce tbb

#endif


// probably the cleanest way to retrieve RcppParallel's concurrency setting
// set by RcppParallel::setThreadOptions(numThreads)
inline size_t get_rcpp_num_threads() {
  auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  return (nullptr == nt_env) 
    ? tbb::task_arena::automatic  // -1
    : static_cast<size_t>(std::atoi(nt_env));
}
