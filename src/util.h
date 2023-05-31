//
//  Copyright (c) 2022 - 2023, Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "config.h"
#include "Rcpp.h"
#include <vector>

void force_output();

std::vector<int> find_desNodes(
    const std::vector< std::vector<double>>& phy_edge,
    int focal);

std::vector<int> find_connections(
    const std::vector< std::vector<double>>& phy_edge,
    int focal);

double get_dt(const std::vector< std::vector<double>>& phy_edge,
              int focal);

void find_desNodes(const std::vector< std::vector<double>>& phy_edge,
                   int focal,
                   std::vector<int>* desNodes,
                   std::vector<double>* timeInte);

double get_time_inte(const std::vector< std::vector<double>>& forTime,
                     int focal_node);

void normalize_loglik_node(std::vector<double>* probvec,
                           long double* loglik);

void normalize_loglik(std::vector<double>* probvec,
                      long double* loglik);


void numericmatrix_to_vector(const Rcpp::NumericMatrix& m,
                             std::vector< std::vector< double >>* v);

void vector_to_numericmatrix(const std::vector< std::vector< double >>& v,
                             Rcpp::NumericMatrix* m);

void output_vec(const std::vector<double>& v);

void list_to_vector(const Rcpp::ListOf<Rcpp::NumericMatrix>& l,
                    std::vector< std::vector< std::vector<double >>>* v);

struct data_storage {
  std::vector<double> t;
  std::vector<std::vector<double>> probs;

  void add_entry(double time, std::vector<double> prob) {
    t.push_back(time);
    probs.push_back(prob);
  }
};

struct entry {
  int ances;
  int focal_node;
  data_storage probabilities;

  entry(int a, int fn, const data_storage& probs) :
    ances(a), focal_node(fn), probabilities(probs)
  {};
};

struct storage {
  std::vector< entry > data_;

  void add_entry(int a, int fn, const data_storage& p) {
    data_.push_back(entry(a, fn, p));
  }
};
