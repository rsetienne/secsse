//
//  mutable_dist.h
//  test_speed_secssesim
//
//  Created by thijsjanzen on 02/10/2025.
//

#ifndef mutable_dist_h
#define mutable_dist_h

#include <vector>
#include <random>

struct mutable_distribution {
    
    mutable_distribution() {
    }
    
    void push_back(const double& val) {
        probs_.push_back(val);
        
        auto x = cdf_.empty() ? 0.0 : cdf_.back();
        
        cdf_.push_back(x + val);
    }
    
    void clear() {
        probs_.clear();
        cdf_.clear();
    }
    
    void change_val(size_t index, double new_val) {
        if (new_val == probs_[index]) return; // no change in value
        probs_[index] = new_val;
        double sum(index > 0 ? cdf_[index - 1] : 0);
        for (size_t i = index; i < cdf_.size(); ++i) {
            cdf_[i] = sum += probs_[i];
        }
    }
    
    template <typename Reng>
    int operator()(Reng& reng) {
      auto p = cdf_.back() * unif_dist(reng);
      return static_cast<int>(std::distance(cdf_.cbegin(), std::lower_bound(cdf_.cbegin(), cdf_.cend(), p)));
    }
    
    std::vector<double> cdf_;
    std::vector<double> probs_;
    std::uniform_real_distribution<> unif_dist = std::uniform_real_distribution<>(0, 1.0);
};

#endif /* mutable_dist_h */
