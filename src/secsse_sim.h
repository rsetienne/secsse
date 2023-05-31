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

using num_mat = std::vector< std::vector<double >>;
using num_mat_mat = std::vector<num_mat>;
using vec_dist  = std::vector< std::discrete_distribution<> >;

enum event_type {shift, speciation, extinction, max_num};

enum finish_type {done, extinct, overshoot, conditioning, not_run_yet,
                  max_types};

struct ltab_species {
  enum info_index {time, p_id, self_id, extinct_time};

  ltab_species(double brts, int parent, int ID, double death) {
    data_[time] = brts;
    data_[p_id] = static_cast<double>(parent);
    data_[self_id] = static_cast<double>(ID);
    data_[extinct_time] = death;
  }

  double get_id() const {
    return(data_[self_id]);
  }

  double get_parent() const  {
    return(data_[p_id]);
  }

  void set_death(double d) {
    data_[extinct_time] = d;
  }

  bool is_dead() {
    if (data_[extinct_time] > -1) return true;
    return false;
  }

  ltab_species() {
    data_[time] = -1e6;
    data_[p_id] = -1e6;
    data_[self_id] = -1e6;
    data_[extinct_time] = -1e6;
  }

  std::array<double, 4>& get_data() {
    return data_;
  }

 private:
  std::array<double, 4> data_;
};

struct ltable {
  std::vector< ltab_species > data_;

  ltable() {
    data_.emplace_back(ltab_species(0.0,  0, -1, -1));
    data_.emplace_back(ltab_species(0.0, -1,  2, -1));
  }

  void clear() {
    data_.clear();
  }
};

struct lambda_dist {
  std::vector<size_t> indices;
  std::vector<double> probs;
  std::discrete_distribution<size_t> d;

  size_t draw_from_dist(std::mt19937* rndgen) {
    return indices[d(*rndgen)];
  }

  lambda_dist(const std::vector<size_t>& i,
              const std::vector<double>& p) : indices(i), probs(p) {
    d = std::discrete_distribution<size_t>(probs.begin(), probs.end());
  }
};

struct species_info {
  species_info() {
    max_la = max_mu = max_qs = 0.0;
  }

  species_info(const std::vector<double>& m,
               const std::vector<double>& l,
               const std::vector<double>& s) :
    trait_mu(m), trait_lambda(l), trait_qs(s) {
    max_mu = *std::max_element(trait_mu.begin(), trait_mu.end());
    max_la = *std::max_element(trait_lambda.begin(), trait_lambda.end());
    max_qs = *std::max_element(trait_qs.begin(), trait_qs.end());
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

  double max_ext()const  {
    return max_mu;
  }

  double max_spec() const {
    return max_la;
  }

  double max_shift() const {
    return max_qs;
  }

 private:
  std::vector<double> trait_mu;
  std::vector<double> trait_lambda;
  std::vector<double> trait_qs;
  double max_mu = 0.0;
  double max_la = 0.0;
  double max_qs = 0.0;
};

struct species {
 private:
  size_t trait_;

 public:
  int id_;
  double mu_;
  double lambda_;
  double shiftprob_;
  double sum_rate_;

  species(size_t trait, int ID, const species_info& info) :
    trait_(trait),  id_(ID), mu_(info.mu(trait)), lambda_(info.lambda(trait)),
    shiftprob_(info.shift(trait)) {
    sum_rate_ = mu_ + lambda_ + shiftprob_;
  }

  void change_trait(size_t new_trait, const species_info& info) {
    trait_ = new_trait;
    mu_ = info.mu(new_trait);
    lambda_ = info.lambda(new_trait);
    shiftprob_ = info.shift(new_trait);
    sum_rate_ = mu_ + lambda_ + shiftprob_;
  }

  size_t get_trait() const {
    return trait_;
  }
};


struct population {
  std::vector<species> pop;
  std::array<double, 3> rates;

  population() {
    rates = {0.0, 0.0, 0.0};
  }

  void add(const species& s) {
    rates[shift]      += s.shiftprob_;
    rates[extinction] += s.mu_;
    rates[speciation] += s.lambda_;
    pop.push_back(s);
  }

  void remove(size_t index) {
    rates[shift]      -= pop[index].shiftprob_;
    rates[extinction] -= pop[index].mu_;
    rates[speciation] -= pop[index].lambda_;

    pop[index] = pop.back();
    pop.pop_back();
  }

  void change_trait(size_t index, size_t new_trait, const species_info& info) {
    auto old_mu    = pop[index].mu_;
    auto old_la    = pop[index].lambda_;
    auto old_shift = pop[index].shiftprob_;
    pop[index].change_trait(new_trait, info);

    rates[shift]      += pop[index].shiftprob_ - old_shift;
    rates[extinction] += pop[index].mu_ - old_mu;
    rates[speciation] += pop[index].lambda_ - old_la;
  }

  bool empty() const {
    return pop.empty();
  }

  size_t size() const {
    return pop.size();
  }

  size_t get_trait(size_t index) const {
    return pop[index].get_trait();
  }

  int get_id(size_t index) const {
    return pop[index].id_;
  }

  auto begin() const {
    return pop.begin();
  }

  auto end() const {
    return pop.end();
  }

  void clear() {
    pop.clear();
    rates = {0.0, 0.0, 0.0};
  }
};

struct secsse_sim {
  std::mt19937 rndgen_;

  ltable L;
  double t;

  finish_type run_info;

  int init_state;

  population pop;

  species_info trait_info;

  std::vector< lambda_dist > lambda_distributions;
  vec_dist qs_dist;

  std::array<int, 2> track_crowns;

  // external data:
  const std::vector<double> mus;
  const size_t num_states;
  const double max_t;
  const size_t max_spec;
  const std::vector<double> init_states;
  const bool non_extinction;

  secsse_sim(const std::vector<double>& m,
             const num_mat_mat& l,
             const num_mat& q,
             double mt,
             size_t max_s,
             const std::vector<double>& init,
             const bool& ne) : mus(m),
             num_states(m.size()), max_t(mt),
             max_spec(max_s),
             init_states(init),
             non_extinction(ne) {
    auto l_sums = update_lambdas(l);
    auto q_sums = update_qs_row_sums(q);
    trait_info = species_info(mus, l_sums, q_sums);

    // randomize randomizer
    std::random_device rd;
    std::mt19937 rndgen_t(rd());
    rndgen_ = rndgen_t;
    run_info = not_run_yet;
    t = 0.0;
    init_state = 0;
  }

  void run() {
    t = 0.0;

    // randomly draw initial trait
    std::uniform_int_distribution<size_t> d(0, init_states.size() - 1);
    init_state = init_states[d(rndgen_)];

    run_info = not_run_yet;

    pop.clear();
    L.clear();

    auto crown_states = root_speciation(init_state);

    pop.add(species(std::get<0>(crown_states), -1, trait_info));
    pop.add(species(std::get<1>(crown_states),  2, trait_info));

    track_crowns = {1, 1};

    L = ltable();

    while (true) {
      double dt = draw_dt();
      t += dt;

      if (t > max_t)  {
        run_info = done; break;
      }

      event_type event = draw_event();
      apply_event(event);

      if (track_crowns[0] < 1 || track_crowns[1] < 1) {
        run_info = extinct;
        break;
      }
      if (pop.size() > max_spec) {
        run_info = overshoot; break;
      }
    }
  }

  void check_rates() {
    // for debugging
    std::vector<double> check_rates(3);
    check_rates[shift] = std::accumulate(pop.pop.begin(), pop.pop.end(),
                    0.0,
                    [](double x, const species& s){return x + s.shiftprob_;});
    check_rates[extinction] = std::accumulate(pop.pop.begin(), pop.pop.end(),
                    0.0,
                    [](double x, const species& s){return x + s.mu_;});
    check_rates[speciation] = std::accumulate(pop.pop.begin(), pop.pop.end(),
                    0.0,
                    [](double x, const species& s){return x + s.lambda_;});

    for (int i = shift; i != max_num; ++i) {
      if (std::abs(check_rates[i] - pop.rates[i]) > 1e-3) {
        std::cerr << t << " " << i << " " <<
                     pop.rates[i] << " " << check_rates[i] << "\n";
        exit(0);
      }
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
    size_t dying = 0;
    if (pop.size() > 1) {
      // sample one at randomly following mus
      auto get_val = [](const species& s) { return s.mu_;};
      dying = sample_from_pop(get_val);
    }
    auto dying_id = pop.get_id(dying);

    for (auto& i : L.data_) {
      if (std::abs(i.get_id()) == std::abs(dying_id)) {
        if (i.get_id() < 0) {
          track_crowns[0]--;
        } else {
          track_crowns[1]--;
        }

        i.set_death(t);
        break;
      }
    }
    pop.remove(dying);
  }

  void event_speciation() {
    size_t mother = 0;
    if (pop.size() > 1) {
      // sample one at randomly following lambdas
      auto get_val = [](const species& s) { return s.lambda_;};
      mother = sample_from_pop(get_val);
    }
    auto mother_trait = pop.get_trait(mother);

    auto pick_speciation_cell = pick_speciation_id(mother_trait);
    auto trait_to_parent      = calc_y(pick_speciation_cell);
    auto trait_to_daughter    = calc_x(pick_speciation_cell);

    pop.change_trait(mother, trait_to_parent, trait_info);

    int new_id = static_cast<int>(L.data_.size()) + 1;
    if (pop.get_id(mother) < 0) {
      track_crowns[0]++;
      new_id *= -1;
    } else {
      track_crowns[1]++;
    }

    pop.add(species(trait_to_daughter, new_id, trait_info));
    L.data_.emplace_back(ltab_species(t, pop.get_id(mother), new_id, -1));
  }

  std::tuple<int, int> root_speciation(int root_state) {
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
    size_t index_chosen_species = 0;
    if (pop.size() > 1) {
      // sample one at randomly following shiftprob
      auto get_val = [](const species& s) { return s.shiftprob_;};
      index_chosen_species = sample_from_pop(get_val);
    }

    auto trait_chosen_species = pop.get_trait(index_chosen_species);

    size_t shift_to = qs_dist[trait_chosen_species](rndgen_);
    pop.change_trait(index_chosen_species, shift_to, trait_info);
    return;
  }

  event_type draw_event() {
    double total_rate = pop.rates[shift] +
                        pop.rates[extinction] +
                        pop.rates[speciation];
    std::uniform_real_distribution<double> unif_dist(0.0, total_rate);
    double r = unif_dist(rndgen_);

    // ordering of rates is:
    // {shift, speciation, extinction, max_num};
     if (r < pop.rates[shift]) return shift;
     if (r < pop.rates[shift] + pop.rates[speciation]) return speciation;

     return extinction;
  }

  double draw_dt() {
    double total_rate = pop.rates[shift] +
                        pop.rates[extinction] +
                        pop.rates[speciation];

    std::exponential_distribution<double> exp_dist(total_rate);
    return exp_dist(rndgen_);
  }

  std::vector<double> update_lambdas(const num_mat_mat& lambdas) {
    std::vector<double> lambda_sums(lambdas.size(), 0.0);

    lambda_distributions.clear();
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
      lambda_distributions.emplace_back(lambda_dist(indices, probs));
    }

    return lambda_sums;
  }

  std::vector<double> update_qs_row_sums(const num_mat& qs) {
    auto qs_row_sums = std::vector<double>(qs.size(), 0.0);
    qs_dist.resize(qs.size());
    for (size_t i = 0; i < qs.size(); ++i) {
      qs_row_sums[i] = std::accumulate(qs[i].begin(), qs[i].end(), 0.0);
      qs_dist[i] = std::discrete_distribution<>(qs[i].begin(), qs[i].end());
    }
    return qs_row_sums;
  }

  size_t sample_from_pop(double (*getvalfrom_species)(const species&)) {
    auto max = *std::max_element(pop.pop.begin(), pop.pop.end(),
                         [&](const species& a, const species& b) {
                           return getvalfrom_species(a) < getvalfrom_species(b);
                         });

    std::uniform_int_distribution<> d(0, static_cast<int>(pop.size()) - 1);
    std::uniform_real_distribution<double> r(0.0, 1.0);
    int index;
    double mult = 1.0 / getvalfrom_species(max);
    double ulim = 1.0 - 1e-9;
    while (true) {
      index = d(rndgen_);
      double rel_prob = getvalfrom_species(pop.pop[index]) * mult;
      if (rel_prob > 0.0) {
        if (rel_prob >= (ulim)) break;

        if (r(rndgen_) < rel_prob) {
          break;
        }
      }
    }
    return index;
  }

  size_t get_num_traits() {
    std::vector<int> hist(mus.size(), 0);
    for (size_t i = 0; i < pop.size(); ++i) {
      auto trait = pop.get_trait(i);
      hist[trait]++;
    }
    size_t cnt = 0;
    for (const auto& i : hist) {
      if (i > 0) cnt++;
    }
    return cnt;
  }

  void check_true_states(size_t num_traits) {
      std::vector<int> focal_traits(num_traits);
      std::iota(focal_traits.begin(), focal_traits.end(), 0);
      for (size_t i = 0; i < pop.size(); ++i) {
        auto trait = pop.get_trait(i);
        for (size_t j = 0; j < focal_traits.size(); ++j) {
          if (focal_traits[j] == trait) {
            focal_traits[j] = focal_traits.back();
            focal_traits.pop_back();
            break;
          }
        }
        if (focal_traits.empty()) {
          break;
        }
      }
      if (focal_traits.empty()) {
        run_info = done;
        return;
      } 

      // otherwise, conditioning is a reason to reject:
      run_info = conditioning;

      return;
  }

  void check_obs_states(size_t num_concealed_states,
                        size_t num_observed_states) {

      std::vector<int> focal_traits; //(num_observed_states);
      for (size_t i = 0; i < num_observed_states; ++i) focal_traits.push_back(i);
      //std::iota(focal_traits.begin(), focal_traits.end(), 0);
      for (size_t i = 0; i < pop.size(); ++i) {
        auto trait = pop.get_trait(i) % num_concealed_states;
        for (size_t j = 0; j < focal_traits.size(); ++j) {
          if (focal_traits[j] == trait) {
            focal_traits[j] = focal_traits.back();
            focal_traits.pop_back();
            break;
          }
        }
        if (focal_traits.empty()) {
          break;
        }
      }
      if (focal_traits.empty()) {
        run_info = done;
        return;
      } 

      // otherwise, conditioning is a reason to reject:
      run_info = conditioning;

      return;
  }

  void check_conditioning(std::string conditioning_type,
                          size_t num_concealed_states,
                          size_t num_states) {

    if (run_info == extinct) return;

    if (conditioning_type == "none") {
      run_info = done;
    }

    if (conditioning_type == "true_states") {
        check_true_states(num_states);
    }

    if (conditioning_type == "obs_states") {
      check_obs_states(num_concealed_states, num_states / num_concealed_states); 
    }

    return;
  }



  void check_num_traits(const std::vector<double>& input_traits) {
    std::vector<double> focal_traits = input_traits;
    if (run_info != done) return;

    // check if all focal traits are there
    if (focal_traits.empty()) return;  // no conditioning

    // now check if each trait to be checked is present:
    for (size_t i = 0; i < pop.size(); ++i) {
      auto trait = pop.get_trait(i);
      for (size_t j = 0; j < focal_traits.size(); ++j) {
        if (focal_traits[j] == trait) {
          focal_traits[j] = focal_traits.back();
          focal_traits.pop_back();
          break;
        }
      }
      if (focal_traits.empty()) {
        break;
      }
    }
    // if traits is empty, all traits were found:
    if (focal_traits.empty()) {
      run_info = done;
      return;
    } 

    // otherwise, conditioning is a reason to reject:
    run_info = conditioning;

    return;
  }

  std::vector<int> get_traits() {
    std::vector<int> traits(pop.size() * 2);
    for (size_t i = 0; i < pop.size(); ++i) {
      auto index = i * 2;
      traits[index] = pop.get_trait(i);
      traits[index + 1] = pop.pop[i].id_;
    }
    return traits;
  }

  size_t get_initial_state() {
    return init_state;
  }

  num_mat extract_ltable() {
    num_mat extracted_ltable(L.data_.size(), std::vector<double>(4));
    for (int i = 0; i < L.data_.size(); ++i) {
      auto temp = L.data_[i].get_data();
      std::vector<double> row(temp.begin(), temp.end());
      extracted_ltable[i] = row;
    }
    return extracted_ltable;
  }
};
