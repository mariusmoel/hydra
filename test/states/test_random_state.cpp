#include "../catch.hpp"

#include "../blocks/electron/testcases_electron.h"
#include <hydra/algebra/algebra.h>
#include <hydra/states/random_state.h>
#include <hydra/utils/close.h>

#include <iostream>
#include <set>

using namespace hydra;

TEST_CASE("random_state", "[states]") {
  using namespace hydra::testcases::electron;

  Log.out("random state Spinhalf distinction test");
  double first_r = 0.;
  complex first_c = 0.;

  // Test whether random states from different blocks are different
  for (int n_sites = 6; n_sites <= 8; ++n_sites) {
    Log("N={}", n_sites);
    auto [group, irreps] = get_cyclic_group_irreps(n_sites);

    for (int nup = 0; nup <= n_sites; ++nup) {
      for (auto irrep : irreps) {
        auto block = Spinhalf(n_sites, nup, group, irrep);

        if (block.size() > 3) {
          // HydraPrint(irrep);
          // HydraPrint(block);
          auto state_real = StateReal(block, RandomState());
          auto state_cplx = StateCplx(block, RandomState());
          // HydraPrint(state_real);
          // HydraPrint(state_real.vector());
          // HydraPrint(state_real.vector());
          // HydraPrint(state_cplx.vector());
          if (first_r == 0.) {
            first_r = state_real(0);
          } else {
            REQUIRE(std::abs(state_real(0) - first_r) > 1e-12);
          }
          if (first_c == 0.) {
            first_c = state_cplx(0);
          } else {
            REQUIRE(std::abs(state_cplx(0) - first_c) > 1e-12);
          }
        }
      }
    }
  }
#ifdef _OPENMP

  // Check whether result with multiple threads is the same as with a single
  // thread
  Log.out("random state Spinhalf omp test");

  auto block = Spinhalf(4);
  for (int seed = 0; seed < 10; ++seed) {
    auto state = StateReal(block, RandomState(seed));
    auto state_cplx = StateCplx(block, RandomState(seed));

    omp_set_num_threads(1);
    auto state2 = StateReal(block, RandomState(seed));
    auto state_cplx2 = StateCplx(block, RandomState(seed));
    REQUIRE(arma::norm(state.vector() - state2.vector()) < 1e-12);
    REQUIRE(arma::norm(state_cplx.vector() - state_cplx2.vector()) < 1e-12);
  }
#endif
}
