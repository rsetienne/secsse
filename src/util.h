// Copyright 2022 - 2023 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//

#pragma once

#include <vector>
#include "Rcpp.h"

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
                   std::vector<int>& desNodes,
                   std::vector<double>& timeInte);

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
