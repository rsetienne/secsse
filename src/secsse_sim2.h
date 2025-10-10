//
//  Copyright (c) 2022 - 2023, Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <memory>
#include <vector>
#include <array>
#include <random>
#include <tuple>
#include <string>
#include <map>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "mutable_dist.h"  // NOLINT [build/include_subdir]

using num_mat     = std::vector< std::vector<double >>;
using num_mat_mat = std::vector<num_mat>;
using vec_dist    = std::vector< std::discrete_distribution<> >;
using event_dist  = std::vector<mutable_distribution>;

enum event_type {shift, speciation, extinction, max_num};

enum finish_type {done, extinct, overshoot, conditioning, not_run_yet,
                  max_types};

struct species_info {
  species_info(const std::vector<double>& m,
               const std::vector<double>& l,
               const std::vector<double>& s) :
    trait_mu(m), trait_lambda(l), trait_qs(s) {
  }


  double mu(size_t trait) const {
    return trait_mu[trait];
  }

  double lambda(size_t trait) const {
    return trait_lambda[trait];
  }

  double shift(size_t trait) const {
    return trait_qs[trait];
  }

 private:
  const std::vector<double> trait_mu;
  const std::vector<double> trait_lambda;
  const std::vector<double> trait_qs;
};


struct ltab_species {
  enum info_index {time, p_id, self_id, extinct_time, trait_val};

  ltab_species(double brts, int parent, int ID, double death,
               double trait, const species_info& spec,
               int ind, event_dist& draw_dist) :
    index(ind) {
    data_[time] = brts;
    data_[p_id] = static_cast<double>(parent);
    data_[self_id] = static_cast<double>(ID);
    data_[extinct_time] = death;
    data_[trait_val] = trait;

    rates[extinction] = spec.mu(trait);
    rates[speciation] = spec.lambda(trait);
    rates[shift]      = spec.shift(trait);

    draw_dist[extinction].push_back(rates[extinction]);
    draw_dist[speciation].push_back(rates[speciation]);
    draw_dist[shift].push_back(rates[shift]);
  }

  double get_id() const {
    return(data_[self_id]);
  }

  double get_parent() const  {
    return(data_[p_id]);
  }

  double get_trait() const {
    return(data_[trait_val]);
  }

  void set_trait(double new_val,
                 const species_info& spec,
                 event_dist& draw_dist) {
    data_[trait_val] = new_val;
    rates[extinction] = spec.mu(new_val);
    rates[speciation] = spec.lambda(new_val);
    rates[shift] = spec.shift(new_val);
    update_dist(draw_dist);
  }

  void set_death(double d, event_dist& draw_dist) {
    data_[extinct_time] = d;
    rates[extinction] = 0.0;
    rates[speciation] = 0.0;
    rates[shift] = 0.0;
    update_dist(draw_dist);
  }

  bool is_dead() const {
    if (data_[extinct_time] > -1) return true;
    return false;
  }

  std::array<double, 3> get_rates() const {
    return rates;
  }

  ltab_species() : index(-1) {
    data_[time] = -1e6;
    data_[p_id] = -1e6;
    data_[self_id] = -1e6;
    data_[extinct_time] = -1e6;
    data_[trait_val] = -1;
    rates[extinction] = 0.0;
    rates[speciation] = 0.0;
    rates[shift] = 0.0;
  }

  std::array<double, 5> get_data() const {
    return data_;
  }

  double mu() const {
    return rates[extinction];
  }
  double lambda() const {
    return rates[speciation];
  }
  double shift_rate() const {
    return rates[shift];
  }

  void update_dist(event_dist& draw_dist) const {
    draw_dist[extinction].change_val(index, rates[extinction]);
    draw_dist[speciation].change_val(index, rates[speciation]);
    draw_dist[shift].change_val(index, rates[shift]);
  }

 private:
  std::array<double, 5> data_;
  std::array<double, 3> rates;
  int index;
};


struct lambda_dist {
  std::vector<size_t> indices;
  std::discrete_distribution<size_t> d;

  size_t draw_from_dist(std::mt19937_64* rndgen) {
    auto index = d(*rndgen);
    return indices[index];
  }

  lambda_dist(const std::vector<size_t>& i,
              const std::vector<double>& p) : indices(i) {
    d = std::discrete_distribution<size_t>(p.begin(), p.end());
  }
};

struct secsse_sim {
  std::mt19937_64 rndgen_;

  std::vector< ltab_species > L;

  std::vector< lambda_dist > lambda_distributions;

  event_dist draw_dist;

  vec_dist qs_dist;
  const species_info trait_info;

  std::array<int, 2> track_crowns;
  std::array<double, 3> rates;

  // external data:
  const std::vector<double> mus;
  const size_t num_states;
  const double max_t;
  const size_t max_spec;

  const bool non_extinction;
  const bool max_spec_extant;
  const std::vector<double> init_state_probs;
  const bool crown_start;

  finish_type run_info;
  int init_state;
  double t;

  const bool verbose;

  secsse_sim(const std::vector<double>& m,
             const num_mat_mat& l,
             const num_mat& q,
             double mt,
             size_t max_s,
             bool max_s_e,
             const std::vector<double>& init,
             const bool& ne,
             int seed,
             bool start_at_crown,
             bool verb) :
             trait_info(m, update_lambdas(l), update_qs_row_sums(q)),
             mus(m),
             num_states(m.size()), max_t(mt),
             max_spec(max_s),
             non_extinction(ne),
             max_spec_extant(max_s_e),
             init_state_probs(init),
             crown_start(start_at_crown),
             run_info(not_run_yet),
             t(0.0),
             verbose(verb) {
    // randomize randomizer
    rndgen_.seed((seed < 0) ? std::random_device {}() : seed);
    init_state = 0;
    track_crowns = {0, 0};
    rates = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < event_type::max_num; ++i) {
        draw_dist.push_back(mutable_distribution());
    }
  }

  void run() {
    t = 0.0;

    // randomly draw initial trait
    std::discrete_distribution<> init_trait_dist(init_state_probs.begin(),
                                                 init_state_probs.end());
    init_state = init_trait_dist(rndgen_);

    run_info = not_run_yet;

    L.clear();
    for (size_t i = 0; i < event_type::max_num; ++i) {
          draw_dist[i].clear();
    }
    int updateFreq = static_cast<int>(max_t / 20);
    if (updateFreq < 1) updateFreq = 1;

    bool print_bar = false;

    if (crown_start) {
      auto crown_states = root_speciation(init_state);
      L.push_back(ltab_species(0.0,   0, -1, -1,
                               std::get<0>(crown_states),
                               trait_info, 0, draw_dist));
      L.push_back(ltab_species(0.0,  -1,  2, -1,
                               std::get<1>(crown_states),
                               trait_info, 1, draw_dist));
    } else {
      L.push_back(ltab_species(0.0,  0, -1, -1, init_state,
                               trait_info, 0, draw_dist));
      track_crowns = {1, 0};
      evolve_until_crown();
      if (t > max_t) {
        run_info = done;
        return;
      }
      if (track_crowns[0] + track_crowns[1] < 1) {
        run_info = extinct;
        return;
      }
    }

    track_crowns = {1, 1};

    int prev_t = 0;

    while (true) {
      update_rates();
      double dt = draw_dt();
      t += dt;

      if (verbose) {
       auto current_t = static_cast<int>(t);
       if (current_t - prev_t > 0) {
          if (current_t % updateFreq == 0) {
              if (print_bar == false) {
                Rcpp::Rcout << "\n0--------25--------50--------75--------100\n";
                Rcpp::Rcout << "*";
                print_bar = true;
              }
              Rcpp::Rcout << "**";
              prev_t = current_t;
          }
        }
      }

      if (t > max_t)  {
        run_info = done; break;
      }

      event_type event = draw_event();
      apply_event(event);

      if (crown_start) {
        if (track_crowns[0] < 1 || track_crowns[1] < 1) {
          run_info = extinct;
          break;
        }
      }

      if (max_spec_extant) {
        if (track_crowns[0] + track_crowns[1] >= max_spec) {
          run_info = overshoot; break;
        }
      } else {
        if (L.size() >= max_spec) {
          run_info = overshoot; break;
        }
      }
    }
  }

  void evolve_until_crown() {
    while (true) {
        update_rates();
        double dt = draw_dt();
        t += dt;
        if (t > max_t) {
          run_info = done;
          return;
        }
        event_type event = draw_event();
        switch (event) {
          case shift: {
              event_traitshift();
              break;
          }
          case extinction: {
            event_extinction();
            break;
          }
          case speciation: {
              size_t mother = 0;

              auto mother_trait = L[mother].get_trait();

              auto pick_speciation_cell = pick_speciation_id(mother_trait);
              auto trait_to_parent      = calc_y(pick_speciation_cell);
              auto trait_to_daughter    = calc_x(pick_speciation_cell);

              L[mother].set_trait(trait_to_parent, trait_info, draw_dist);

              L.emplace_back(t, L[mother].get_id(), 2, -1,
                             trait_to_daughter, trait_info, 1, draw_dist);
              break;
          }
          default:
              break;
        }

        // stop when we reach two species
        if (L.size() >= 2) break;

        // stop if all are extinct
        if (track_crowns[0] + track_crowns[1] < 1) break;
      }
  }

  void apply_event(const event_type event) {
    switch (event) {
      case shift: {
        event_traitshift();
        break;
      }
      case speciation: {
        event_speciation();
        break;
      }
      case extinction: {
        event_extinction();
        break;
      }
      default: break;
    }
    return;
  }

  void event_extinction() {
    size_t dying = sample_from_pop(event_type::extinction);

    if (L[dying].get_id() < 0) {
          track_crowns[0]--;
    } else {
          track_crowns[1]--;
    }

    L[dying].set_death(t, draw_dist);
  }

  void event_speciation() {
    size_t mother = sample_from_pop(event_type::speciation);

    auto mother_trait = L[mother].get_trait();

    auto pick_speciation_cell = pick_speciation_id(mother_trait);
    auto trait_to_parent      = calc_y(pick_speciation_cell);
    auto trait_to_daughter    = calc_x(pick_speciation_cell);

    L[mother].set_trait(trait_to_parent, trait_info, draw_dist);

    int new_id = static_cast<int>(L.size()) + 1;
    if (L[mother].get_id() < 0) {
      track_crowns[0]++;
      new_id *= -1;
    } else {
      track_crowns[1]++;
    }

    L.emplace_back(t, L[mother].get_id(), new_id, -1,
                   trait_to_daughter, trait_info, L.size(), draw_dist);
  }

  std::tuple<int, int> root_speciation(int root_state) {
    // check of rates:
    size_t num_empty = 0;
    for (const auto& i : lambda_distributions) {
      if (i.indices.empty()) num_empty++;
    }
    if (num_empty == lambda_distributions.size()) {
      Rcpp::stop("all speciation rates are zero - this makes it impossible to create a crown from a root state, \nplease set one rate to a very low non-zero rate");   // NOLINT.
    }

    // calculate if the other crown lineage is the same trait:
    auto mother_trait = root_state;

    auto pick_speciation_cell = pick_speciation_id(mother_trait);
    auto trait_to_parent      = calc_y(pick_speciation_cell);
    auto trait_to_daughter    = calc_x(pick_speciation_cell);

    return {trait_to_parent, trait_to_daughter};
  }

  size_t calc_x(size_t index) {
    return index % num_states;
  }

  size_t calc_y(size_t index) {
    return index / num_states;
  }

  int pick_speciation_id(const size_t& index) {
    return lambda_distributions[index].draw_from_dist(&rndgen_);
  }

  void event_traitshift() {
    size_t index_chosen_species = sample_from_pop(event_type::shift);

    auto trait_chosen_species = L[index_chosen_species].get_trait();

    size_t shift_to = qs_dist[trait_chosen_species](rndgen_);

    L[index_chosen_species].set_trait(shift_to, trait_info, draw_dist);
    return;
  }

  void update_rates() {
    rates[extinction] = draw_dist[extinction].cdf_.back();
    rates[speciation] = draw_dist[speciation].cdf_.back();
    rates[shift] = draw_dist[shift].cdf_.back();
  }

  event_type draw_event() {
    double total_rate = rates[shift] +
                        rates[extinction] +
                        rates[speciation];
    std::uniform_real_distribution<double> unif_dist(0.0, total_rate);
    double r = unif_dist(rndgen_);

    // ordering of rates is:
    // {shift, speciation, extinction, max_num};
     if (r < rates[shift]) return shift;
     if (r < (rates[shift] + rates[speciation])) return speciation;

     return extinction;
  }

  double draw_dt() {
    double total_rate = rates[shift] +
                        rates[extinction] +
                        rates[speciation];
    std::exponential_distribution<double> exp_dist(total_rate);
    return exp_dist(rndgen_);
  }

  std::vector<double> update_lambdas(const num_mat_mat& lambdas) {
    std::vector<double> lambda_sums(lambdas.size(), 0.0);

    std::vector<size_t> indices;
    std::vector<double> probs;

    for (size_t m = 0; m < lambdas.size(); ++m) {
      size_t index = 0;
      indices.clear();
      probs.clear();
      for (const auto& i : lambdas[m]) {
        for (const auto& j : i) {
          lambda_sums[m] += j;
          if (j > 0) {
            indices.push_back(index);
            probs.push_back(j);
          }
          index++;
        }
      }
      lambda_distributions.emplace_back(indices, probs);
    }

    return lambda_sums;
  }

  std::vector<double> update_qs_row_sums(const num_mat& qs) {
    auto qs_row_sums = std::vector<double>(qs.size(), 0.0);
    qs_dist.resize(qs.size());
    for (size_t i = 0; i < qs.size(); ++i) {
      qs_row_sums[i] = std::accumulate(qs[i].begin(), qs[i].end(), 0.0);
      // qs_dist[i] is not accessed if qs_row_sums[i] == 0.0
      qs_dist[i] = std::discrete_distribution<>(qs[i].begin(), qs[i].end());
    }
    return qs_row_sums;
  }

  int sample_from_pop(event_type event) {
      return draw_dist[event](rndgen_);
  }

  void check_states(size_t num_traits,
                    size_t num_concealed_states) {
    auto total_num_traits = num_concealed_states > 0 ?
                              num_traits / num_concealed_states :
                              num_traits;

    std::vector<int> focal_traits;
    for (size_t i = 0; i < total_num_traits; ++i) focal_traits.push_back(0);

    for (const auto& i : L) {
      int trait = static_cast<int>(i.get_trait());
      if (num_concealed_states > 0) trait = trait % num_concealed_states;
      focal_traits[trait]++;
    }

    auto min_val = *std::min_element(focal_traits.begin(),
                                     focal_traits.end());
    if (min_val == 0) {
      run_info = conditioning;
    } else {
      run_info = done;
    }

    return;
  }

  void check_custom_conditioning(const std::vector<double>& condition_vec,
                                 int num_concealed_traits) {
    std::map<int, int> histogram;

    for (const auto& i : L) {
      int trait = static_cast<int>(i.get_trait()) % num_concealed_traits;
        histogram[trait]++;
    }

    for (const auto& c : condition_vec) {
      if (histogram.find(c) == histogram.end()) {
        run_info = conditioning;
        return;
      }
    }

    run_info = done;
    return;
  }

  std::vector<int> get_traits() {
    std::vector<int> traits(L.size() * 2);
    for (size_t i = 0; i < L.size(); ++i) {
      auto index = i * 2;
      traits[index] = L[i].get_trait();
      traits[index + 1] = L[i].get_id();
    }
    return traits;
  }

  void check_conditioning(std::string conditioning_type,
                          size_t num_concealed_states,
                          size_t num_states,
                          const std::vector<double>& condition_vec) {
    if (run_info == extinct) return;

    if (conditioning_type == "none") {
      run_info = done;
    }

    if (conditioning_type == "true_states") {
        check_states(num_states, 0);
    }

    if (conditioning_type == "obs_states") {
        check_states(num_states, num_concealed_states);
    }

    if (conditioning_type == "custom") {
      // do something
      check_custom_conditioning(condition_vec, num_concealed_states);
    }

    return;
  }

  size_t get_initial_state() {
    return init_state;
  }

  size_t num_species() {
    if (max_spec_extant) {
      return track_crowns[0] + track_crowns[1];
    } else {
      return L.size();
    }
  }

  size_t ltable_size() {
    return L.size();
  }

  num_mat extract_ltable() {
    num_mat extracted_ltable(L.size(), std::vector<double>(5));
    for (size_t i = 0; i < L.size(); ++i) {
      auto temp = L[i].get_data();
      std::vector<double> row(temp.begin(), temp.end());
      extracted_ltable[i] = row;
    }
    return extracted_ltable;
  }

  void update_tree_size_hist(int* val) {
    if (run_info == extinct) {
      *val = 0;
      return;
    }

    if (max_spec_extant) {
       *val = (track_crowns[0] + track_crowns[1]);
    } else {
       *val = static_cast<int>(L.size());
    }
    return;
  }
};
