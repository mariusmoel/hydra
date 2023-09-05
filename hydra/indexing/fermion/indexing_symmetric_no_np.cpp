#include "indexing_symmetric_no_np.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/symmetries/operations/representative_list.h>
#include <hydra/symmetries/operations/fermi_sign.h>

namespace hydra::indexing::fermion {

template <class bit_t>
IndexingSymmetricNoNp<bit_t>::IndexingSymmetricNoNp(
    int n_sites, PermutationGroup permutation_group, Representation irrep)
    : n_sites_(n_sites),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep),
      subsets_indexing_(n_sites), fermi_table_(n_sites_, permutation_group) {

  assert(n_sites >= 0);
  utils::check_n_sites(n_sites, permutation_group);

  std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
      symmetries::representatives_indices_symmetries_limits_norms<bit_t>(
          subsets_indexing_, group_action_, irrep);

  size_ = (idx_t)reps_.size();
  begin_ = iterator_t(reps_, 0);
  end_ = iterator_t(reps_, size_);
}

template <typename bit_t>
bool IndexingSymmetricNoNp<bit_t>::operator==(
    IndexingSymmetricNoNp<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (group_action_ == rhs.group_action_) &&
         (irrep_ == rhs.irrep_);
}

template <typename bit_t>
bool IndexingSymmetricNoNp<bit_t>::operator!=(
    IndexingSymmetricNoNp<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class IndexingSymmetricNoNp<uint16_t>;
template class IndexingSymmetricNoNp<uint32_t>;
template class IndexingSymmetricNoNp<uint64_t>;

} // namespace hydra::indexing::fermion
