#include "lanczos_eigenvector.h"

#include <hydra/random/hash_functions.h>
#include <hydra/random/hashes.h>
#include <hydra/random/random_utils.h>
#include <hydra/states/random_state.h>

namespace hydra {

std::pair<Tmatrix, arma::Col<double>> lanczos_eigenvector_real(
    BondList const &bonds, Block const &block, int num_eigenvector,
    double precision, int max_iterations, uint64_t seed, double deflation_tol) {
  return lanczos_eigenvector<double>(bonds, block, num_eigenvector, precision,
                                     max_iterations, seed, deflation_tol);
}

std::pair<Tmatrix, arma::Col<complex>> lanczos_eigenvector_cplx(
    BondList const &bonds, Block const &block, int num_eigenvector,
    double precision, int max_iterations, uint64_t seed, double deflation_tol) {
  return lanczos_eigenvector<complex>(bonds, block, num_eigenvector, precision,
                                      max_iterations, seed, deflation_tol);
}

template <class coeff_t>
std::pair<Tmatrix, arma::Col<coeff_t>>
lanczos_eigenvector(BondList const &bonds, Block const &block,
                    int num_eigenvector, double precision, int max_iterations,
                    uint64_t seed, double deflation_tol) {

  // Create random starting vector with normal distributed entries
  auto set_v0 = [&seed, &block](arma::Col<coeff_t> &v0) {
    uint32_t seed_modified = random::hash_combine(seed, block.hash());
    random::fill_random_normal_vector(v0, seed_modified);
  };

  // Run Lanczos algorithm
  return lanczos_eigenvector<coeff_t>(bonds, block, set_v0, num_eigenvector,
                                      precision, max_iterations, deflation_tol);
}

template <class coeff_t>
std::pair<Tmatrix, arma::Col<coeff_t>>
lanczos_eigenvector(BondList const &bonds, State<coeff_t> const &v0,
                    int num_eigenvector, double precision, int max_iterations,
                    double deflation_tol) {
  auto set_v0 = [&v0](arma::Col<coeff_t> &v0_copy) { v0_copy = v0.vector(); };
  return lanczos_eigenvector<coeff_t>(bonds, v0.block(), set_v0,
                                      num_eigenvector, precision,
                                      max_iterations, deflation_tol);
}

template std::pair<Tmatrix, arma::Col<double>>
lanczos_eigenvector(BondList const &, State<double> const &, int, double, int,
                    double);
template std::pair<Tmatrix, arma::Col<complex>>
lanczos_eigenvector(BondList const &, State<complex> const &, int, double, int,
                    double);

} // namespace hydra
