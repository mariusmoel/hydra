#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/utils/bitops.h>

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra::terms::electron_symmetric_simple {

template <class bit_t, class GroupAction, class Filler>
void do_ising_symmetric(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetricSimple<bit_t, GroupAction> const &block, Filler &&fill) {

  using bitops::gbit;
  using bitops::popcnt;

  auto ising = bonds.bonds_of_type("HEISENBERG") +
               bonds.bonds_of_type("ISING") + bonds.bonds_of_type("HB");
  auto ising_tj = bonds.bonds_of_type("TJHEISENBERG") +
                  bonds.bonds_of_type("TJISING") + bonds.bonds_of_type("TJHB");

  for (auto bond : ising + ising_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing tJ Ising: "
                    "bond must have exactly two sites defined");

    if (!utils::coupling_is_zero(bond, couplings)) {

      double J = lila::real(couplings[bond.coupling()]);

      // Set values for same/diff (tJ block definition)
      std::string type = bond.type();
      double val_same, val_diff;
      if ((type == "HEISENBERG") || (type == "ISING") || (type == "HB")) {
        val_same = J / 4.;
        val_diff = -J / 4.;
      } else if ((type == "TJHEISENBERG") || (type == "TJISING") ||
                 (type == "TJHB")) {
        val_same = 0.;
        val_diff = -J / 2.;
      }

      int s1 = bond.site(0);
      int s2 = bond.site(1);
      bit_t s1mask = (bit_t)1 << s1;
      bit_t s2mask = (bit_t)1 << s2;
      bit_t mask = s1mask | s2mask;

      for (auto [up, lower_upper] : block.ups_lower_upper_) {
        idx_t lower = lower_upper.first;
        idx_t upper = lower_upper.second;

        if ((up & mask) == mask) { // both spins pointing up
          for (idx_t idx = lower; idx < upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if (!(dn & mask))
              fill(idx, idx, val_same);
          }
        } else if (up & s1mask) { // s1 is pointing up
          for (idx_t idx = lower; idx < upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if ((dn & mask) == s2mask)
              fill(idx, idx, val_diff);
          }
        } else if (up & s2mask) { // s2 is pointing up
          for (idx_t idx = lower; idx < upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if ((dn & mask) == s1mask)
              fill(idx, idx, val_diff);
          }
        } else { // no upspins
          for (idx_t idx = lower; idx < upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if ((dn & mask) == mask)
              fill(idx, idx, val_same);
          }
        }
      }
    }
  }
}

} // namespace hydra::terms::electron_symmetric_simple