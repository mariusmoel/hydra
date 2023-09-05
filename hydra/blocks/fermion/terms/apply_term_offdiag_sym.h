#pragma once

#include <hydra/common.h>

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::fermion {

template <typename bit_t, typename coeff_t, class IndexingIn, class IndexingOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_offdiag_sym_to_fermions(bit_t fermions_in, idx_t idx_in,
                                     std::vector<coeff_t> const &characters,
                                     IndexingIn &&indexing_in,
                                     IndexingOut &&indexing_out,
                                     NonZeroTerm &&non_zero_term,
                                     TermAction &&term_action, Fill &&fill) {

  if (non_zero_term(fermions_in)) {
    auto [fermions_out, coeff] = term_action(fermions_in);
    auto [idx_out, sym] = indexing_out.index_sym(fermions_out);
    
    if (idx_out != invalid_index) {
      double norm_out = indexing_out.norm(idx_out);
      double norm_in = indexing_in.norm(idx_in);
      coeff_t bloch = characters[sym];
      coeff_t val = coeff * bloch * norm_out / norm_in;
      fill(idx_out, idx_in, val);
    }
  }
}

template <typename bit_t, typename coeff_t, class IndexingIn, class IndexingOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_offdiag_sym(IndexingIn &&indexing_in,
                            IndexingOut &&indexing_out,
                            NonZeroTerm &&non_zero_term,
                            TermAction &&term_action, Fill &&fill) {
  Representation irrep_out = indexing_out.irrep();
  std::vector<coeff_t> characters;
  if constexpr (is_complex<coeff_t>()) {
    characters = irrep_out.characters();
  } else {
    characters = irrep_out.characters_real();
  }

#ifdef _OPENMP
  idx_t size = indexing_in.size();

#pragma omp parallel for schedule(guided)
  for (idx_t idx_in = 0; idx_in < size; ++idx_in) {
    bit_t fermions_in = indexing_in.state(idx_in);
    apply_term_offdiag_sym_to_fermions(fermions_in, idx_in, characters, indexing_in,
                                    indexing_out, non_zero_term, term_action,
                                    fill);
  }

#else
  // Run through all states and apply term
  for (auto [fermions_in, idx_in] : indexing_in) {
    apply_term_offdiag_sym_to_fermions(fermions_in, idx_in, characters, indexing_in,
                                    indexing_out, non_zero_term, term_action,
                                    fill);
  }
#endif
}

} // namespace hydra::fermion
