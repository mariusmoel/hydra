#include "fermion_apply.h"

#include <hydra/blocks/fermion/terms/apply_terms_dispatch.h>
#include <hydra/blocks/fermion/terms/compile.h>
#include <hydra/operators/compiler.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, Fermion const &block_in,
           arma::Col<coeff_t> const &vec_in, Fermion const &block_out,
           arma::Col<coeff_t> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == (idx_t)vec_in.size());
  assert(block_out.size() == (idx_t)vec_out.size());

  // get Bondlist out of generic/special bonds -> see compile.cpp
  BondList bonds_c = fermion::compile(bonds, 1e-12);
  operators::check_bonds_in_range(bonds, block_in.n_sites());

  if ((is_real<coeff_t>()) && (bonds_c.is_complex())) {
    Log.err("Error in matrix_gen: trying to create a real matrix from an "
            "intrisically complex BondList");
  }

  vec_out.zeros();
  //smart lambda
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, coeff_t val) {
#ifdef _OPENMP
    if constexpr (is_real<coeff_t>()) {
      coeff_t x = val * vec_in(idx_in);
      coeff_t *pos = vec_out.memptr();
#pragma omp atomic update
      pos[idx_out] += x;
    } else {
      complex x = val * vec_in(idx_in);
      complex *pos = vec_out.memptr();
      double *r = &reinterpret_cast<double(&)[2]>(pos[idx_out])[0];
      double *i = &reinterpret_cast<double(&)[2]>(pos[idx_out])[1];
#pragma omp atomic update
      *r += x.real();
#pragma omp atomic update
      *i += x.imag();
    }
#else
   // Matrix-vector Multiplication 
    //std::cout << vec_out(idx_out)  << " " << vec_in(idx_in) << std::endl;
    vec_out(idx_out) += val * vec_in(idx_in);
    //std::cout << " " << idx_in <<" "<< idx_out << " " << vec_out(idx_out) << " " << val << " " << vec_in(idx_in) << std::endl;
#endif
  };

  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();
  fermion::apply_terms_dispatch<coeff_t>(bonds_c, indexing_in, indexing_out,
                                          fill);
}

template void apply<double>(BondList const &, Fermion const &,
                            arma::Col<double> const &, Fermion const &,
                            arma::Col<double> &);

template void apply<complex>(BondList const &, Fermion const &,
                             arma::Col<complex> const &, Fermion const &,
                             arma::Col<complex> &);

} // namespace hydra
