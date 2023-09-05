#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/fermion/fermion.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

arma::Mat<double> matrix_real(BondList const &bonds, Fermion const &block_in,
                              Fermion const &block_out);

arma::Mat<complex> matrix_cplx(BondList const &bonds, Fermion const &block_in,
                               Fermion const &block_out);

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds, Fermion const &block_in,
                              Fermion const &block_out);

} // namespace hydra
