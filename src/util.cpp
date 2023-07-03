//
//  Copyright (c) 2022 - 2023, Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include "config.h"    // NOLINT [build/include_subdir]
#include "util.h"   // NOLINT [build/include_subdir]

std::vector<int> find_desNodes(
    const std::vector< std::vector<double>>& phy_edge,
    int focal) {
  std::vector<int> output;
  for (size_t i = 0; i < phy_edge.size(); ++i) {
    if (phy_edge[i][0] == focal) {
      output.push_back(phy_edge[i][1]);
    }
  }
  return(output);
}

double get_dt(const std::vector< std::vector<double>>& phy_edge,
              int focal) {
  for (size_t i = 0; i < phy_edge.size(); ++i) {
    if (phy_edge[i][1] == focal) {
      return phy_edge[i][2];
    }
  }
  return 0.0;
}

void find_desNodes(const std::vector< std::vector<double>>& phy_edge,
                   int focal,
                   std::vector<int>* desNodes,
                   std::vector<double>* timeInte) {
  (*desNodes).resize(2);
  (*timeInte).resize(2);
  size_t cnt = 0;
  for (size_t i = 0; i < phy_edge.size(); ++i) {
    if (phy_edge[i][0] == focal) {
      (*desNodes)[cnt] = phy_edge[i][1];
      (*timeInte)[cnt] = phy_edge[i][2];
      cnt++;
    }
    if (cnt > 1) break;
  }
}

std::vector<int> find_connections(
    const std::vector< std::vector<double>>& phy_edge,
    int focal) {
  std::vector<int> output(2);
  int cnt = 0;
  for (size_t i = 0; i < phy_edge.size(); ++i) {
    if (phy_edge[i][0] == focal) {
     output[cnt] = phy_edge[i][1];
     cnt++;
    }
    if (cnt >= 2) break;
  }
  return output;
}


double get_time_inte(const std::vector< std::vector<double>>& forTime,
                     int focal_node) {
  //  R code:   timeInte <- forTime[which(forTime[,2] == desNodes[desIndex]), 3]
  for (const auto& i : forTime) {
    if (i[1] == focal_node) {
      return(i[2]);
    }
  }
  return 0.0;
}

void normalize_loglik_node(std::vector<double>* probvec,
                           long double* loglik) {
  size_t d = (*probvec).size() / 2;

  double sumabsprobs(0.0);
  for (size_t i = d; i < (d + d); ++i) {
    sumabsprobs += std::abs((*probvec)[i]);
  }
  for (size_t i = d; i < (d + d); ++i) {
    (*probvec)[i] *= 1.0 / sumabsprobs;
  }
  (*loglik) += log(sumabsprobs);
  return;
}

void normalize_loglik(std::vector<double>* probvec,
                      long double* loglik) {
  static const auto abssum = [] (auto x, auto y) {return x + std::abs(y);};

  double sumabsprobs = std::accumulate((*probvec).begin(), (*probvec).end(),
                                       0.0,
                                       abssum);

  if (sumabsprobs > 0.0) {
      for (auto& i : (*probvec)) {
        i *= 1.0 / sumabsprobs;
      }
      (*loglik) += log(sumabsprobs);
  }
  return;
}

void numericmatrix_to_vector(const Rcpp::NumericMatrix& m,
                             std::vector< std::vector< double >>* v) {
  (*v) = std::vector< std::vector< double> >(m.nrow(),
            std::vector<double>(m.ncol(), 0.0));
  for (int i = 0; i < m.nrow(); ++i) {
    std::vector<double> row(m.ncol(), 0.0);
    for (int j = 0; j < m.ncol(); ++j) {
      row[j] = m(i, j);
    }
    (*v)[i] = row;
  }
  return;
}

std::vector< std::vector< double> >
    num_mat_to_vec(const Rcpp::NumericMatrix& m) {
  auto v = std::vector< std::vector< double> >(m.nrow(),
                                               std::vector<double>(m.ncol(),
                                                                   0.0));
  for (int i = 0; i < m.nrow(); ++i) {
    std::vector<double> row(m.ncol(), 0.0);
    for (int j = 0; j < m.ncol(); ++j) {
      row[j] = m(i, j);
    }
    v[i] = row;
  }
  return v;
}

std::vector< std::vector< std::vector<double >>>
  list_to_vector(const Rcpp::ListOf<Rcpp::NumericMatrix>& ll) {
  std::vector< std::vector< std::vector< double > >> ll_cpp;
  for (size_t i = 0; i < ll.size(); ++i) {
    Rcpp::NumericMatrix temp = ll[i];
    std::vector< std::vector< double >> temp2;
    for (size_t j = 0; j < temp.nrow(); ++j) {
      std::vector<double> row;
      for (size_t k = 0; k < temp.ncol(); ++k) {
        row.push_back(temp(j, k));
      }
      temp2.push_back(row);
    }
    ll_cpp.push_back(temp2);
  }

  return ll_cpp;
}

void vector_to_numericmatrix(const std::vector< std::vector< double >>& v,
                             Rcpp::NumericMatrix* m) {
  size_t n_rows = v.size();
  size_t n_cols = v[0].size();
  (*m) = Rcpp::NumericMatrix(n_rows, n_cols);
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      (*m)(i, j) = v[i][j];
    }
  }
  return;
}

void output_vec(const std::vector<double>& v) {
  // std::cerr << "vec: ";
  //  for (size_t i = 0; i < v.size(); ++i) {
  //    std::cerr << v[i] << " ";
  //} std::cerr << "\n";
}
