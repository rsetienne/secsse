//
//  naive.h
//  secsse_sim
//
//  Created by thijsjanzen on 21/02/2023.
//

#pragma once

#include <memory>
#include <vector>
#include <array>
#include <random>


struct ltab_species {
  ltab_species(double brts, int parent, int ID, double death){
    data_[0] = brts;
    data_[1] = static_cast<double>(parent);
    data_[2] = static_cast<double>(ID);
    data_[3] = death;
  }
  
  double get_id() {
    return(data_[2]);
  }
  double get_parent() {
    return(data_[1]);
  }
  void set_death(double d) {
    data_[3] = d;
  }
  
  
  std::array<double, 4> data_;
};

struct ltable {
  std::vector< ltab_species > data_;
  
  ltable() {
    data_.emplace_back(ltab_species(0.0, 0, -1, -1));
    data_.emplace_back(ltab_species(0.0, -1, 2, -1));
  }
  
  void clear() {
    data_.clear();
  }
};

struct lambda_dist {
  std::vector<size_t> indices;
  std::vector<double> probs;
  std::discrete_distribution<size_t> d;
  
  size_t draw_from_dist(std::mt19937& rndgen) {
    return indices[ d(rndgen) ];
  }
  
  lambda_dist(const std::vector<size_t>& i,
              const std::vector<double>& p) : indices(i), probs(p) {
    d = std::discrete_distribution<size_t>(probs.begin(), probs.end());
  }
};



using num_mat = std::vector< std::vector<double >>;
using num_mat_mat = std::vector<num_mat>;

enum event_type {shift, speciation, extinction, max_num};
struct species_info {
  species_info() {
    max_la = max_mu = max_qs = 0.0;
  }
  
  species_info(const std::vector<double>& m,
               const std::vector<double>& l,
               const std::vector<double>& s) : trait_mu(m), trait_lambda(l), trait_qs(s) {
    
    max_mu = *std::max_element(trait_mu.begin(), trait_mu.end());
    max_la = *std::max_element(trait_lambda.begin(), trait_lambda.end());
    max_qs = *std::max_element(trait_qs.begin(), trait_qs.end());
    
  //  for (const auto& i : trait_lambda) {
  //    std::cerr << i << " ";
  //  } std::cerr << "\n";
    
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
  
  double max_ext() {
    return max_mu;
  }
  
  double max_spec() {
    return max_la;
  }
  
  double max_shift() {
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
  int id_;
  double mu_;
  double lambda_;
  double shiftprob_;
  double sum_rate_;
  
  species(size_t trait, int ID, const species_info& info) :
    trait_(trait),  id_(ID), mu_(info.mu(trait)), lambda_(info.lambda(trait)), shiftprob_(info.shift(trait)) {
    sum_rate_ = mu_ + lambda_ + shiftprob_;
  }
  
  void change_trait(size_t new_trait, const species_info& info) {
    trait_ = new_trait;
    mu_ = info.mu(new_trait);
    lambda_ = info.lambda(new_trait);
    shiftprob_ = info.shift(new_trait);
    sum_rate_ = mu_ + lambda_ + shiftprob_;
  }
  
  size_t get_trait() {
    return trait_;
  }
  
private:
  size_t trait_;
};


struct population {
  std::vector<species> pop;
  std::array<double, 3> rates;
  
  population() {
    rates = {0.0, 0.0, 0.0};
  }
  
  void add(const species& s) {
    pop.emplace_back(s);
    rates[shift] += s.shiftprob_;
    rates[extinction] += s.mu_;
    rates[speciation] += s.lambda_;
  }
  
  void remove(size_t index) {
    rates[shift] -= pop[index].shiftprob_;
    rates[extinction] -= pop[index].mu_;
    rates[speciation] -= pop[index].lambda_;
    
    pop[index] = pop.back();
    pop.pop_back();
  }
  
  void change_trait(size_t index, size_t new_trait, const species_info& info) {
    auto old_mu = pop[index].mu_;
    auto old_la = pop[index].lambda_;
    auto old_shift = pop[index].shiftprob_;
    pop[index].change_trait(new_trait, info);
    
    rates[shift] += pop[index].shiftprob_ - old_shift;
    rates[extinction] += pop[index].mu_ - old_mu;
    rates[speciation] += pop[index].lambda_ - old_la;
  }
  
  bool empty() {
    return pop.empty();
  }
  
  size_t size() {
    return pop.size();
  }
  
  size_t get_trait(size_t index) {
    return pop[index].get_trait();
  }
  
  int get_id(size_t index) {
    return pop[index].id_;
  }
  
  void clear() {
    pop.clear();
    rates = {0.0, 0.0, 0.0};
  }
};

using vec_dist  = std::vector< std::discrete_distribution<> >;


struct secsse_sim {
  
  std::mt19937 rndgen_;
  
  ltable L;
  double t;
  bool extinct;
  int init_state;
  
  population pop;
  
  species_info trait_info;
  
  std::vector< lambda_dist > lambda_distributions;
  // num_mat lambdas_vectorized;
  // vec_dist lambdas_vectorized_distributions;
  vec_dist qs_dist;
  
  // external data:
  const std::vector<double> mus;
  //  const num_mat qs;
  // const num_mat_mat lambdas;
  const size_t num_states;
  const double max_t;
  const size_t max_spec;
  
  
  secsse_sim(const std::vector<double>& m,
             const num_mat_mat& l,
             const num_mat& q,
             double mt,
             size_t max_s) : mus(m),
             num_states(m.size()), max_t(mt),
             max_spec(max_s) {
    auto l_sums = update_lambdas(l);
    auto q_sums = update_qs_row_sums(q);
    trait_info = species_info(mus, l_sums, q_sums);
    
    // randomize randomizer
    std::random_device rd;
    std::mt19937 rndgen_t(rd());
    rndgen_ = rndgen_t;
    extinct = false;
  }
  
  void run(int initial_trait) {
   
    t = 0.0;
    
    if (initial_trait < 0) {
      // randomly draw initial trait
      std::uniform_int_distribution<size_t> d(0, mus.size() - 1);
      initial_trait = d(rndgen_);
    }
    
    init_state = initial_trait;
    extinct = false;
    
    pop.clear();
    L.clear();
    pop.add(species(initial_trait, -1, trait_info));
    pop.add(species(initial_trait,  2, trait_info));
    L = ltable();
    // initial rates:
  //  std::cerr << initial_trait << "\n";
//    for (const auto& i : pop.pop) {
//      std::cerr << i.mu_ << " " << i.lambda_ << " " << i.shiftprob_ << " " << i.sum_rate_ << "\n";  
//    }
    
    while(!pop.empty() &&
          t < max_t) { // these conditions are never exceeded, see below:

      double dt = draw_dt();
      if (t + dt > max_t) break;
      t += dt;
      
      event_type event = draw_event();
      apply_event(event);
      
      if (pop.size() > max_spec) break;
    }
  }
  
  void check_rates() {
    std::vector<double> check_rates(3);
    check_rates[shift] = std::accumulate(pop.pop.begin(), pop.pop.end(), 0.0,
                                         [](double x, const species& s){return x + s.shiftprob_;});
    check_rates[extinction] = std::accumulate(pop.pop.begin(), pop.pop.end(), 0.0,
                                              [](double x, const species& s){return x + s.mu_;});
    check_rates[speciation] = std::accumulate(pop.pop.begin(), pop.pop.end(), 0.0,
                                              [](double x, const species& s){return x + s.lambda_;});
    
    for (int i = shift; i != max_num; ++i) {
      if (std::abs(check_rates[i] - pop.rates[i]) > 1e-3) {
        std::cout << t << " " << i << " " << pop.rates[i] << " " << check_rates[i] << "\n";
        exit(0);
      }
    }
    
  }
  
  void apply_event(const event_type event) {
    switch(event) {
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
      dying = sample_from_pop(get_val, trait_info.max_ext());
    }
    auto dying_id = pop.get_id(dying);
    for (auto& i : L.data_) {
      if (i.get_id() == dying_id) {
        i.set_death(t);
        break;
      }
    }
    
    pop.remove(dying);
  }
  
  void event_speciation() {
    /// this needs some checking with all the indices!
    size_t mother = 0;
    if (pop.size() > 1) {
      // sample one at randomly following lambdas
      auto get_val = [](const species& s) { return s.lambda_;};
      mother = sample_from_pop(get_val, trait_info.max_spec());
    }
    auto mother_trait = pop.get_trait(mother);
    
    auto pick_speciation_cell = pick_speciation_id(mother_trait);
    auto trait_to_parent      = calc_y(pick_speciation_cell);
    auto trait_to_daughter    = calc_x(pick_speciation_cell);
    
    pop.change_trait(mother, trait_to_parent, trait_info);
    
    int new_id = static_cast<int>(L.data_.size()) + 1;
    if (pop.get_id(mother) < 0) new_id *= -1;
    
    pop.add(species(trait_to_daughter, new_id, trait_info));
    L.data_.emplace_back(ltab_species(t, pop.get_id(mother), new_id, -1));
  }
  
  
  size_t calc_x(size_t index) {
    return index % num_states;
  }
  size_t calc_y(size_t index) {
    return index / num_states;
  }
  
  int pick_speciation_id(const size_t& index) {
    //return lambdas_vectorized_distributions[index](rndgen_);
    return lambda_distributions[index].draw_from_dist(rndgen_);
  }
  
  
  void event_traitshift() {
    size_t index_chosen_species = 0;
    if (pop.size() > 1) {
      // sample one at randomly following shiftprob
      auto get_val = [](const species& s) { return s.shiftprob_;};
      index_chosen_species = sample_from_pop(get_val, trait_info.max_spec());
    }
    
    auto trait_chosen_species = pop.get_trait(index_chosen_species);
    
    size_t shift_to = qs_dist[trait_chosen_species](rndgen_);  //d(rndgen_);
    pop.change_trait(index_chosen_species, shift_to, trait_info);
    return;
  }
  
  event_type draw_event() {
    double total_rate = pop.rates[shift] + pop.rates[extinction] + pop.rates[speciation];
    std::uniform_real_distribution<double> unif_dist(0.0, total_rate);
    double r = unif_dist(rndgen_);
    
    // ordering of rates is:
    // {shift, speciation, extinction, max_num};
    
    if (r < pop.rates[shift]) return shift;
    if (r < pop.rates[shift] + pop.rates[speciation]) return speciation;
    
    return extinction;
  }
  
  
  double draw_dt() {
    double total_rate = pop.rates[shift] + pop.rates[extinction] + pop.rates[speciation];

    std::exponential_distribution<double> exp_dist(total_rate);
    return exp_dist(rndgen_);
  }
  
  std::vector<double> update_lambdas(const num_mat_mat& lambdas) {
    std::vector<double> lambda_sums(lambdas.size(), 0.0);
    
    lambda_distributions.clear(); //reserve(lambdas.size());
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
    auto qs_row_sums = std::vector<double>(qs.size(), 0.0) ;
    qs_dist.resize(qs.size());
    for (size_t i = 0; i < qs.size(); ++i) {
      qs_row_sums[i] = std::accumulate(qs[i].begin(), qs[i].end(), 0.0);
      //qs_dist.push_back(std::discrete_distribution<>(qs[i].begin(), qs[i].end()));
      qs_dist[i] = std::discrete_distribution<>(qs[i].begin(), qs[i].end());
    }
    return qs_row_sums;
  }
  
  size_t sample_from_pop(double (*getvalfrom_species)(const species&), double max) {
    // first we get the max
    
    std::uniform_int_distribution<> d(0, static_cast<int>(pop.size()) - 1);
    std::uniform_real_distribution<double> r(0.0, 1.0);
    int index;
    double mult = 1.0 / max;
    double ulim = 1.0 - 1e-9;
    while(true) {
      index = d(rndgen_);
      double rel_prob = getvalfrom_species(pop.pop[index]) *mult;
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
      std::vector<double> row(L.data_[i].data_.size());
      for (int j = 0; j < L.data_[i].data_.size(); ++j) {
        row[j] = L.data_[i].data_[j];
      }
      extracted_ltable[i] = row;
    }
    return extracted_ltable;
  }
};
