#include "indexing_np.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations_index.h>

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::indexing::fermion {

template <typename bit_t>
IndexingNp<bit_t>::IndexingNp(int n_sites, int n)
    : n_sites_(n_sites), n_(n), lintable_(n_sites, n),
      states_(combinatorics::binomial(n_sites, n)), size_(states_.size()),
      begin_(n_sites, n, 0), end_(n_sites, n, size_) {
  utils::check_n_fermion(n_sites, n, "Fermion");

#ifdef _OPENMP
  for (auto [state, idx] :
       combinatorics::CombinationsIndexThread<bit_t>(n_sites, n)) {
    states_[idx] = state;
  }
#else
  for (auto [state, idx] :
       combinatorics::CombinationsIndex<bit_t>(n_sites, n)) {
    states_[idx] = state;
  }
#endif
}

template <typename bit_t>
bool IndexingNp<bit_t>::operator==(IndexingNp<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_ == rhs.n_);
}

template <typename bit_t>
bool IndexingNp<bit_t>::operator!=(IndexingNp<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class IndexingNp<uint16_t>;
template class IndexingNp<uint32_t>;
template class IndexingNp<uint64_t>;

} // namespace hydra::indexing::fermion
