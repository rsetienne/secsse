## Flawed 'for testing' stuff

Whatever you are testing, this doen't belong to the non-testing
code. And it doesn't look correct.

```R
# secsee_loglik && cla_secsee_loglik.R
 if (num_concealed_states != round(num_concealed_states)) {
    # for testing
    d <- ncol(states) / 2
    new_states <- states[, c(1:sqrt(d), (d + 1):((d + 1) + sqrt(d) - 1))]
    new_states <- new_states[, c(1, 2, 3, 10, 11, 12)]
    states <- new_states
  }
```

## Follow up of the R/C++ indexing nightmare

```R
# cla_secsee_loglik.R
  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    # last row contains safety entry from C++ (all zeros)
    ancestral_states <- states[(num_tips + 1):(nrow(states) - 1), ]
    ancestral_states <-
        ancestral_states[, -1 * (1:(ncol(ancestral_states) / 2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, LL = LL, states states))
  }

# secsee_loglik.R
  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):(nrow(states)), ]
    ancestral_states <-
      ancestral_states[, -1 * (1:(ncol(ancestral_states) / 2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, LL = LL, states = states))
  }
```

Note especially the comment in the `cla` version...
Btw., what about:

```R
  if (see_ancestral_states == TRUE) {
    ancestral_states <- states[(phy$Nnode + 2):nrow(states), (ncol(states) + 1) - ((ncol(states) / 2):1)]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, LL = LL, states = states))
  }
```

## see_ancestral_states

```R
# xxx_loglik.R
return(list(ancestral_states = ancestral_states, LL = LL, states = states))
```

Why both (`ancestral_states` and `states`)? Forces `calc_ll` to return full states.

## Naming conventions

E.g. `q_matrix` in `secsse_loglik.R` vs `Q` in `cla_secsse_loglik,R`

## Differences in ode_cla_x are strange

## Overall very poor rhs performance

## `build_initStates_time` and `build_states` still very slow

Best case: calls `ape::branching.times(phy)` which is slow R-code.<br>
Btw., `build_states` calls `build_initStates_time` (double calculation).
C++ side should accept the `phy` object from the beginning.<br>
As a first step, allow `do_call_ll` to return a list.

## 'Full storage' is not feasable for 'controled' steppers

This **will** blow up memory, recorded states are not in (time) order, duplicates, etc.

## test: missing package `testit`

## `bs_time` leftover

* `cla_secsee_store.cpp::calc_ll_cla_store`

## Misleading comments / argument names

* `secsse_loglik_eval.R`: `ancestral_states = ll$states` not `ancestral states`
* `cla_secsse_eval.R`: `ancestral_states = ll$states` not `ancestral states`

Fix: rename `see_ancestral_states` to `see_states` in C++ code.

## Deja vu

```R
# secsse_plot.R
    calcul <- c()
    ancescpp <- ances - 1
    forTimecpp <- for_time # nolint
    forTimecpp[, c(1, 2)] <- forTimecpp[, c(1, 2)] - 1 # nolint
```

