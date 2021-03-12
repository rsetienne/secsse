#ifndef util_h
#define util_h

#include <vector>
#include "Rcpp.h"

void force_output();

std::vector<int> find_desNodes(const std::vector< std::vector<double>>& phy_edge,
                               int focal);

void find_desNodes(const std::vector< std::vector<double>>& phy_edge,
                   int focal,
                   std::vector<int>& desNodes,
                   std::vector<double>& timeInte);

double get_time_inte(const std::vector< std::vector<double>>& forTime,
                     int focal_node);

std::vector<double> normalize_loglik_node(std::vector<double>& probvec,
                                          double& loglik,
                                          size_t d);

void normalize_loglik(std::vector<double>& probvec,
                      double& loglik);

void numericmatrix_to_vector(const Rcpp::NumericMatrix& m,
                             std::vector< std::vector< double >>& v);

void vector_to_numericmatrix(const std::vector< std::vector< double >>& v,
                             Rcpp::NumericMatrix& m);

#endif /* util_h */