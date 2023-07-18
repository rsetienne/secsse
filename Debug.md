## secsse_sim.h

No need to improve performance but still...

* Missing includes

```C++
#include <algotithm>
#include <numeric>
```

* Poor seeding

```C++
// randomize randomizer
std::random_device rd;
if (seed < 0) seed = rd();
std::mt19937 rndgen_t(seed);
rndgen_ = rndgen_t;
```

```C++
// shorter version, still poor seeding
rndgen_.seed((seed < 0) ? std::random_device{}() : seed);
```

* Use `std::mt19937_64` for drawing `doubles` (64 Bit).
* Setter/getter pattern is rather useless for structs.
* `draw_event` is `std::discrete_distribution` re-invented ;)
* Are all instances of std::discrete_distribution<> hardened againts 'all-zero' and 'empty' (both UB)?<br>
i.e. `secsse_sim::qs_dist` are copied from 'user-space' -- and users are evil.
* `size_t sample_from_pop(double (*getvalfrom_species)(const species&))` prevents inlining of the lambda.
* `lambda_dist::probs` is unused member.
* `L.data_.emplace_back(ltab_species(0.0, 0, -1, -1, pop.get_trait(0)));` undermines the purpose of `emplace_back`.
* `lambda_distributions.emplace_back(lambda_dist(indices, probs));` undermines the purpose of `emplace_back`.
* `lambda_dist::lambda_dist()` use move-construction.
* `secsse_sim::secsse_sim()` and other ctors: move member initialization out of body.
* `secsse_sim::apply_event()` could be replaced by a call table:

```C++
// declare member-function-pointer-table somewhere.
// this is a std::initializer_list<correct mbrfn-pointer-type nobody wants to figure>{...}.
const auto event_handler = { &secsse_sim::event_traitshift, &secsse_sim::event_speciation, &secsse_sim::event_extinction };
// ...
while (true) {
    // call the hanlder through the mbrfn ptrs.
    (this->*std::data(event_handler)[draw_event()])();
}
```
Ja, I know - ugly. But only syntax-wise ;)

* Don't build vectors iterative if you know it's size, e.g. `secsse_sim::check_states`.
* `std::min_element` is an exhausting search, `std::any_of` shortcuts. Thus:

```C++
  // orig
  void check_states(size_t num_traits,
                    size_t num_concealed_states) {
    
    auto total_num_traits = num_concealed_states > 0 ? num_traits / num_concealed_states : num_traits;

    std::vector<int> focal_traits;
    for (size_t i = 0; i < total_num_traits; ++i) focal_traits.push_back(0);    // see above, *might* be optimized by the compiler

    for (const auto& i : L.data_) {
      int trait = static_cast<int>(i.get_trait());
      if (num_concealed_states > 0) trait %= num_concealed_states;
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
```

```C++
  // proposed fiddled version
  void check_states(size_t num_traits,
                    size_t num_concealed_states) {
      auto total_num_traits = num_concealed_states > 0 ? num_traits / num_concealed_states : num_traits;
      std::vector<bool> missing_traits(total_num_traits, true);     // consider as member 
      for (const auto& i : L.data_) {
        auto trait = static_cast<size_t>(i.get_trait()) % total_num_traits;     // eq. signed, branchless
        missing_traits[trait] = false;
      }
      run_info = std::any_of(missing_traits.cbegin(), missing_traits.cend(), [](const auto& missing) { return missing; })
               ? conditioning
               : done;
  }
```

* `population::add(const species& s)` shall accept RValues and shall use `pop.push_back(std::move(...))`.
* RValue/move semantic underused in general.
