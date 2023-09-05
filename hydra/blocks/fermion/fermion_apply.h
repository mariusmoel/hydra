#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/fermion/fermion.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, Fermion const &block_in,
           arma::Col<coeff_t> const &vec_in, Fermion const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace hydra
