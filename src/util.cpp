#include <thread>
#include <chrono>

#include "util.h"

void force_output() {
  //  std::this_thread::sleep_for(std::chrono::nanoseconds(100));
  std::this_thread::sleep_for(std::chrono::milliseconds(30));
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}


std::vector<int> find_desNodes(const std::vector< std::vector<double>>& phy_edge,
                               int focal) {
  std::vector<int> output;
  for(int i = 0; i < phy_edge.size(); ++i) {
    if (phy_edge[i][0] == focal) {
      output.push_back( phy_edge[i][1]);
    }
  }
  return(output);
}

void find_desNodes(const std::vector< std::vector<double>>& phy_edge,
                   int focal,
                   std::vector<int>& desNodes,
                   std::vector<double>& timeInte) {
  std::vector<int> output;
  for(int i = 0; i < phy_edge.size(); ++i) {
    if (phy_edge[i][0] == focal) {
      desNodes.push_back( phy_edge[i][1]);
      timeInte.push_back( phy_edge[i][2]);
    }
  }
}

double get_time_inte(const std::vector< std::vector<double>>& forTime,
                     int focal_node) {
  //     timeInte <- forTime[which(forTime[,2] == desNodes[desIndex]), 3]

  for (const auto& i : forTime) {
    if (i[1] == focal_node) {
      //      double output = ;
      return(i[2]);
    }
  }
  return 0.0;
  // stop("could not find time inte");
}

std::vector<double> normalize_loglik_node(std::vector<double>& probvec,
                                          double& loglik,
                                          size_t d) {
  // 4 5 6 7
  // // d = 4
  double sumabsprobs = 0.0;
  for(size_t i = d; i < (d + d); ++i) {
    sumabsprobs += std::abs(probvec[i]);
  }
  for(size_t i = d; i < (d + d); ++i) {
    probvec[i] *= 1.0 / sumabsprobs;
  }
  loglik = loglik + log(sumabsprobs);
  return probvec;
}

void normalize_loglik(std::vector<double>& probvec,
                      double& loglik) {
  static const auto abssum = [] (auto x, auto y) {return x + std::abs(y);};

  double sumabsprobs = std::accumulate(probvec.begin(), probvec.end(), 0.0,
                                       abssum);
  for (auto& i : probvec) {
    i *= 1.0 / sumabsprobs;
  }
  loglik = loglik + log(sumabsprobs);
  return;
}

void numericmatrix_to_vector(const Rcpp::NumericMatrix& m,
                             std::vector< std::vector< double >>& v) {

  v = std::vector< std::vector< double>>(m.nrow());
  for (int i = 0; i < m.nrow(); ++i) {
    std::vector<double> row(m.ncol());
    for (int j = 0; j < m.ncol(); ++j) {
      row[j] = m(i, j);
    }
    v[i] = row;
  }
  return;
}

void vector_to_numericmatrix(const std::vector< std::vector< double >>& v,
                             Rcpp::NumericMatrix& m) {

  int n_rows = v.size();
  int n_cols = v[0].size();
  m = Rcpp::NumericMatrix(n_rows, n_cols);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      m(i, j) = v[i][j];
    }
  }
  return;
}
