#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/fermion/terms/apply_hopping.h>
#include <hydra/blocks/fermion/terms/apply_number.h>


#include <hydra/common.h>

#include <hydra/utils/print_macro.h>
#include <hydra/utils/timing.h>

namespace hydra::fermion {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_terms(BondList const &bonds, IndexingIn const &indexing_in,
                 IndexingOut const &indexing_out, Fill &fill) {
  for (auto bond : bonds) {

    //if (bond.type_defined()) {
      if (bond.type() == "HOP") {
        fermion::apply_hopping<bit_t, coeff_t, symmetric>(bond, indexing_in,
                                                            indexing_out, fill);
      /*} else if (bond.type() == "C") {
        fermion::apply_c<bit_t, coeff_t, symmetric>(bond, indexing_in,
                                                         indexing_out, fill);
      } else if (bond.type() == "CDAG") {
        fermion::apply_cdag<bit_t, coeff_t, symmetric>(bond, indexing_in,
                                                      indexing_out, fill);
      */
      } else if (bond.type() == "NUMBER") {
        fermion::apply_number<bit_t, coeff_t, symmetric>(bond, indexing_in,
                                                        indexing_out, fill);
      } 
    //} else {
    //  fermion::apply_non_branching<bit_t, coeff_t, symmetric>(
    //      bond, indexing_in, indexing_out, fill);
    //}
  }
}

} // namespace hydra::fermion
