#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_symmetric_hopping(Bond const &bond, Indexing &&indexing,
                             Filler &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert((bond.type() == "HOPUP") || (bond.type() == "HOPDN"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  // Get group/irrep info
  auto const &group_action = indexing.group_action();
  auto const &irrep = indexing.irrep();
  std::vector<coeff_t> bloch_factors;
  if constexpr (is_complex<coeff_t>()) {
    bloch_factors = irrep.characters();
  } else {
    bloch_factors = irrep.characters_real();
  }
  std::vector<coeff_t> prefacs(irrep.size());

  assert(group_action.n_symmetries() == (int)bloch_factors.size());
  int n_sites = group_action.n_sites();

  coeff_t t = bond.coupling<coeff_t>();
  int s1 = bond[0];
  int s2 = bond[1];

  // Prepare bitmasks
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int l = std::min(s1, s2);
  int u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  // // DEBUG PRINT states
  // Log("\n\n\n");
  // Log("s1: {} s2: {}", s1, s2);
  // idx_t idx = 0;
  // for (idx_t idx_up = 0; idx_up < indexing.n_reps_up(); ++idx_up) {
  //   bit_t ups = indexing.rep_up(idx_up);
  //   auto syms = indexing.syms_up(ups);
  //   auto const &dnss = indexing.dns_for_up_rep(ups);
  //   auto const &norms = indexing.norms_for_up_rep(ups);

  //   if (syms.size() == 1) {
  //     bit_t not_ups = (~ups) & sitesmask;

  //     idx_t idx_dn = 0;
  //     for (bit_t dnsc : dnss) {
  //       bit_t dns = bitops::deposit(dnsc, not_ups);
  //       Log("{}: {};{} {}", idx, BSTR(ups), BSTR(dns),
  //       norms[idx_dn]);
  //       ++idx_dn;
  //       ++idx;
  //     }
  //   } else {
  //     idx_t idx_dn = 0;
  //     for (bit_t dns : dnss) {
  //       Log("{}: {};{} {}", idx, BSTR(ups), BSTR(dns),
  //       norms[idx_dn]);

  //       ++idx_dn;
  //       ++idx;
  //     }
  //   }
  // }

  // Apply hoppings on dnspins
  if (bond.type() == "HOPDN") {

    for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
      bit_t ups = indexing.rep_ups(idx_ups);

      // skip if no hopping for dns works
      if (bitops::popcnt(ups & flipmask) != 0) {
        continue;
      }

      idx_t up_offset = indexing.ups_offset(idx_ups);
      auto dnss = indexing.dns_for_ups_rep(ups);
      auto syms = indexing.syms_ups(ups);

      // trivial stabilizer of ups -> dns have to be deposited
      if (syms.size() == 1) {
        bit_t not_ups = (~ups) & sitesmask;
        idx_t idx_in = up_offset;
        for (bit_t dnsc : dnss) {
          bit_t dns = bitops::deposit(dnsc, not_ups);

          // If hopping is possible ...
          if (bitops::popcnt(dns & flipmask) == 1) {
            bit_t dns_flip = dns ^ flipmask;
            bit_t dns_flip_c = bitops::extract(dns_flip, not_ups);
            idx_t idx_out = up_offset + indexing.dnsc_index(dns_flip_c);

            // Complex conjugate t if necessary
            coeff_t val;
            if constexpr (is_complex<coeff_t>()) {
              val = -(bitops::gbit(dns, s1) ? t : hydra::conj(t));
            } else {
              val = -t;
            } // Comment: norm is always 1.0 for trivial stabilizers

            // fill with correct fermi sign
            bool fermi = bitops::popcnt(dns & fermimask) & 1;

            // Log("-----------------------------------------");
            // Log("HOPDN");
            // Log("CASE: stab-ups-FALSE");
            // Log("from: {};{}", BSTR(ups), BSTR(dns));
            // Log("mask: {} {}", BSTR(0), BSTR(flipmask));
            // Log("to  : {};{}", BSTR(ups), BSTR(dns_flip));
            // Log("fermi: {}", fermi);
            // Log("idx_in: {}, idx_out: {}", idx_in, idx_out);
            // Log("val: {}", real(val));
            // Log("fill: {}", real((fermi) ? -val : val));

            fill(idx_out, idx_in, (fermi) ? -val : val);
          }

          ++idx_in;
        }
      }
      // non-trivial stabilizer of ups -> dns don't have to be deposited
      else {
        auto norms = indexing.norms_for_ups_rep(ups);

        idx_t idx_in = up_offset;
        idx_t idx_dns = 0;
        for (bit_t dns : dnss) {

          if (bitops::popcnt(dns & flipmask) == 1) {
            bit_t dns_flip = dns ^ flipmask;
            auto [idx_dns_flip, fermi_dn, sym] =
                indexing.index_dns_fermi_sym(dns_flip, syms, dnss, fermimask);

            if (idx_dns_flip != invalid_index) {

              idx_t idx_out = up_offset + idx_dns_flip;
              bool fermi_up = indexing.fermi_bool_ups(sym, ups);

              // Complex conjugate t if necessary
              coeff_t val;
              if constexpr (is_complex<coeff_t>()) {
                val = -(bitops::gbit(dns, s1) ? t : hydra::conj(t)) *
                      bloch_factors[sym] * norms[idx_dns_flip] / norms[idx_dns];
              } else {
                // Log("sz: {}, sym: {}", bloch_factors.size(), sym);
                val = -t * bloch_factors[sym] * norms[idx_dns_flip] /
                      norms[idx_dns];
              }

              // Log("-----------------------------------------");
              // Log("HOPDN");
              // Log("CASE: stab-ups-TRUE");
              // Log("from: {};{}", BSTR(ups), BSTR(dns));
              // Log("mask: {} {}", BSTR(0), BSTR(flipmask));
              // Log("to  : {};{}", BSTR(ups), BSTR(dns_flip));
              // Log("sym : {}", sym);
              // bit_t dns_flip_rep = group_action.apply(sym, dns_flip);
              // Log("rep : {};{}", BSTR(ups), BSTR(dns_flip_rep));
              // Log("fermi_up: {}, fermi_dn: {}", fermi_up, fermi_dn);
              // Log("idx_in: {}, idx_out: {}", idx_in, idx_out);
              // Log("val: {}", real(val));
              // Log("fill: {}",
              //           real((fermi_up ^ fermi_dn) ? -val : val));

              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            }
          }

          ++idx_in;
          ++idx_dns;
        } // non-trivial stabilizer of ups

      } // for (auto [up, lower_upper]
    }
  } //     if ((hop.type() == "HOP") || (hop.type() == "HOPDN")) {

  // Apply hoppings on upspins
  if (bond.type() == "HOPUP") {

    for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
      bit_t ups = indexing.rep_ups(idx_ups);

      // continue if hopping not possible
      if (bitops::popcnt(ups & flipmask) != 1) {
        continue;
      }

      bit_t ups_flip = ups ^ flipmask;
      idx_t idx_ups_flip = indexing.index_ups(ups_flip);
      bit_t ups_flip_rep = indexing.rep_ups(idx_ups_flip);
      bit_t not_ups_flip_rep = (~ups_flip_rep) & sitesmask;

      // Get limits, syms, and dns for ingoing ups
      idx_t ups_offset_in = indexing.ups_offset(idx_ups);
      auto syms_ups_in = indexing.syms_ups(ups);
      auto dnss_in = indexing.dns_for_ups_rep(ups);

      // Get limits, syms, and dns for outgoing ups
      idx_t ups_offset_out = indexing.ups_offset(idx_ups_flip);
      auto syms_ups_out = indexing.syms_ups(ups_flip);
      auto dnss_out = indexing.dns_for_ups_rep(ups_flip_rep);

      // Trivial stabilizer of target ups
      if (syms_ups_out.size() == 1) {
        int sym = syms_ups_out.front();

        // Complex conjugate t if necessary
        coeff_t prefac;
        if constexpr (is_complex<coeff_t>()) {
          prefac = -(bitops::gbit(ups, s1) ? t : hydra::conj(t)) *
                   bloch_factors[sym];
        } else {
          prefac = -t * bloch_factors[sym];
        }

        // Fermi-sign of up spins
        bool fermi_up = (bitops::popcnt(ups & fermimask) & 1);
        fermi_up ^= indexing.fermi_bool_ups(sym, ups_flip);

        // Origin ups trivial stabilizer -> dns need to be deposited
        if (syms_ups_in.size() == 1) {
          idx_t idx_in = ups_offset_in;
          bit_t not_ups = (~ups) & sitesmask;
          for (bit_t dnsc : dnss_in) {
            bit_t dns = bitops::deposit(dnsc, not_ups);
            if (bitops::popcnt(dns & flipmask) == 0) {
              bit_t dns_rep = group_action.apply(sym, dns);
              bit_t dns_rep_c = bitops::extract(dns_rep, not_ups_flip_rep);
              idx_t idx_out = ups_offset_out + indexing.dnsc_index(dns_rep_c);
              bool fermi_dn = indexing.fermi_bool_dns(sym, dns);

              // Log("-----------------------------------------");
              // Log("HOPUP");
              // Log("CASE: stab-origin-FALSE stab-target-FALSE");
              // Log("from: {};{}", BSTR(ups), BSTR(dns));
              // Log("mask: {} {}", BSTR(flipmask), BSTR(0));
              // Log("to  : {};{}", BSTR(ups_flip), BSTR(dns));
              // Log("rep : {};{}", BSTR(ups_flip_rep), BSTR(dns_rep));
              // Log("idx_in: {}, idx_out: {}", idx_in, idx_out);
              // Log("val: {}", real(prefac));
              // Log("fill: {}",
              //           real((fermi_up ^ fermi_dn) ? -prefac :
              //           prefac));

              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -prefac : prefac);
            }
            ++idx_in;
          }
        }
        // Origin ups have stabilizer -> dns DONT need to be deposited
        else {
          auto norms_in = indexing.norms_for_ups_rep(ups);
          idx_t idx_dn = 0;
          idx_t idx_in = ups_offset_in;
          for (bit_t dns : dnss_in) {
            if (bitops::popcnt(dns & flipmask) == 0) {

              auto [idx_dn_out, fermi_dn] =
                  indexing.index_dns_fermi(dns, sym, not_ups_flip_rep);
              coeff_t val = prefac / norms_in[idx_dn];
              idx_t idx_out = ups_offset_out + idx_dn_out;

              // Log("-----------------------------------------");
              // Log("HOPUP");
              // Log("CASE: stab-origin-TRUE stab-target-FALSE");
              // Log("from: {};{}", BSTR(ups), BSTR(dns));
              // Log("mask: {} {}", BSTR(flipmask), BSTR(0));
              // Log("to  : {};{}", BSTR(ups_flip), BSTR(dns));
              // Log("sym : {}", sym);
              // bit_t dns_rep = group_action.apply(sym, dns);
              // Log("rep : {};{}", BSTR(ups_flip_rep), BSTR(dns_rep));
              // Log("fermi_up: {}, fermi_dn: {}", fermi_up, fermi_dn);
              // Log("idx_in: {}, idx_out: {}", idx_in, idx_out);
              // Log("val: {}", real(val));
              // Log("fill: {}",
              //           real((fermi_up ^ fermi_dn) ? -val : val));

              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            }
            ++idx_dn;
            ++idx_in;
          }
        }
      }
      // Target ups have non-trivial stabilizer
      else {
        auto norms_out = indexing.norms_for_ups_rep(ups_flip_rep);
        auto syms = syms_ups_out;

        // Fix the bloch/prefactors
        if constexpr (is_complex<coeff_t>()) {
          for (int i = 0; i < (int)irrep.size(); ++i) {
            prefacs[i] = -(bitops::gbit(ups, s1) ? t : hydra::conj(t)) *
                         irrep.character(i);
          }
        } else {
          for (int i = 0; i < (int)irrep.size(); ++i) {
            prefacs[i] = -t * real(irrep.character(i));
          }
        }

        bool fermi_up_hop = (bitops::popcnt(ups & fermimask) & 1);

        // Origin ups trivial stabilizer -> dns need to be deposited
        if (syms_ups_in.size() == 1) {
          idx_t idx_in = ups_offset_in;
          bit_t not_ups = (~ups) & sitesmask;
          for (bit_t dnsc : dnss_in) {
            bit_t dns = bitops::deposit(dnsc, not_ups);
            if (bitops::popcnt(dns & flipmask) == 0) {

              auto [idx_dn_out, fermi_dn, sym] =
                  indexing.index_dns_fermi_sym(dns, syms, dnss_out);

              if (idx_dn_out != invalid_index) {
                idx_t idx_out = ups_offset_out + idx_dn_out;
                bool fermi_up =
                    fermi_up_hop ^ indexing.fermi_bool_ups(sym, ups_flip);
                coeff_t val = prefacs[sym] * norms_out[idx_dn_out];

                // Log("-----------------------------------------");
                // Log("HOPUP");
                // Log("CASE: stab-origin-FALSE stab-target-TRUE");
                // Log("from: {};{}", BSTR(ups), BSTR(dns));
                // Log("mask: {} {}", BSTR(flipmask), BSTR(0));
                // Log("to  : {};{}", BSTR(ups_flip), BSTR(dns));
                // Log("sym : {}", sym);
                // bit_t dns_rep = group_action.apply(sym, dns);
                // Log("repl: {};{}", BSTR(ups_flip_rep),
                // BSTR(dns_rep)); Log("fermi_up: {}, fermi_dn: {}",
                // fermi_up, fermi_dn); Log("idx_in: {}, idx_out: {}",
                // idx_in, idx_out); Log("val: {}", real(val));
                // Log("fill: {}",
                //           real((fermi_up ^ fermi_dn) ? -val : val));

                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
            }
            ++idx_in;
          }
        }
        // Origin ups non-trivial stabilizer -> dns DONT need to be deposited
        else {
          auto norms_in = indexing.norms_for_ups_rep(ups);
          idx_t idx_in = ups_offset_in;
          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {
            if (bitops::popcnt(dns & flipmask) == 0) {

              auto [idx_dn_out, fermi_dn, sym] =
                  indexing.index_dns_fermi_sym(dns, syms, dnss_out);
              if (idx_dn_out != invalid_index) {
                idx_t idx_out = ups_offset_out + idx_dn_out;
                bool fermi_up =
                    fermi_up_hop ^ indexing.fermi_bool_ups(sym, ups_flip);
                coeff_t val =
                    prefacs[sym] * norms_out[idx_dn_out] / norms_in[idx_dn];

                // Log("-----------------------------------------");
                // Log("HOPUP");
                // Log("CASE: stab-origin-TRUE stab-target-TRUE");
                // Log("from: {};{}", BSTR(ups), BSTR(dns));
                // Log("mask: {} {}", BSTR(flipmask), BSTR(0));
                // Log("to  : {};{}", BSTR(ups_flip), BSTR(dns));
                // Log("sym : {}", sym);
                // bit_t dns_rep = group_action.apply(sym, dns);
                // Log("rep : {};{}", BSTR(ups), BSTR(dns_rep));
                // Log("fermi_up: {}, fermi_dn: {}", fermi_up,
                // fermi_dn); Log("idx_in: {}, idx_out: {}", idx_in,
                // idx_out); Log("val: {}", real(val));
                // Log("fill: {}",
                //           real((fermi_up ^ fermi_dn) ? -val : val));

                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
            }
            ++idx_dn;
            ++idx_in;
          }
        }
      } // Target ups have non-trivial stabilizer
    }   // loop over ups
  }     // type == "HOPUP"
}

} // namespace hydra::tj
