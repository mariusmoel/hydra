#include "indexing_no_np.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing::fermion {

template <typename bit_t>
IndexingNoNp<bit_t>::IndexingNoNp(int n_sites)
    : n_sites_(n_sites), size_(pow(2, n_sites)), begin_(0), end_(size_) {
  assert(n_sites_ >= 0);
}

template <typename bit_t>
bool IndexingNoNp<bit_t>::operator==(IndexingNoNp<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_);
}

template <typename bit_t>
bool IndexingNoNp<bit_t>::operator!=(IndexingNoNp<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class IndexingNoNp<uint16_t>;
template class IndexingNoNp<uint32_t>;
template class IndexingNoNp<uint64_t>;

} // namespace hydra::indexing::fermion
