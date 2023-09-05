#pragma once

#include <functional>

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/fermion/terms/apply_term_offdiag_no_sym.h>
#include <hydra/blocks/fermion/terms/apply_term_offdiag_sym.h>

namespace hydra::fermion {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_hopping(Bond const &bond, IndexingIn &&indexing_in,
                    IndexingOut &&indexing_out, Fill &&fill) {

  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "HOP"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  int s1 = bond[0];
  int s2 = bond[1];
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int l = std::min(s1, s2);
  int u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);
  coeff_t t = bond.coupling<coeff_t>();

  auto non_zero_term = [&flipmask](bit_t const &spins) -> bool {
    return bitops::popcnt(spins & flipmask) & 1;
  };

  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bool fermi = bitops::popcnt(spins & fermimask) & 1;
    spins ^= flipmask;
    if constexpr (is_complex<coeff_t>()) {
      coeff_t tt = (bitops::gbit(spins, s1)) ? t : hydra::conj(t);
      return {spins, fermi ? tt : -tt};
    } else {
      return {spins, fermi ? t : -t};
    }
  };


  // Dispatch either symmetric of unsymmetric term application
  if constexpr (symmetric) {
    fermion::apply_term_offdiag_sym<bit_t, coeff_t>(
        indexing_in, indexing_out, non_zero_term, term_action, fill);
  } else {
    fermion::apply_term_offdiag_no_sym<bit_t, coeff_t>(
        indexing_in, indexing_out, non_zero_term, term_action, fill);
  }
}

} // namespace hydra::fermion
