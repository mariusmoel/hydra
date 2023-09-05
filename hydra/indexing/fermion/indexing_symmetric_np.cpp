#include "indexing_symmetric_np.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/symmetries/operations/representative_list.h>
#include <hydra/utils/logger.h>
#include <hydra/indexing/fermi_table.h>

namespace hydra::indexing::fermion {

template <class bit_t>
IndexingSymmetricNp<bit_t>::IndexingSymmetricNp(
    int n_sites, int n, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), n_(n),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep),
      combinations_indexing_(n_sites, n),  fermi_table_fermions_(n_sites, n, allowed_subgroup(permutation_group, irrep)) {

  utils::check_n_fermion(n_sites, n, "IndexingSymmetricNp");
  utils::check_n_sites(n_sites, permutation_group);

  std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
      symmetries::representatives_indices_symmetries_limits_norms<bit_t>( // CHANGE HERE!!!
          combinations_indexing_, group_action_, irrep);

  size_ = (idx_t)reps_.size();
  begin_ = iterator_t(reps_, 0);
  end_ = iterator_t(reps_, size_);
}

template <typename bit_t>
bool IndexingSymmetricNp<bit_t>::operator==(
    IndexingSymmetricNp<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_ == rhs.n_) &&
         (group_action_ == rhs.group_action_) && (irrep_ == rhs.irrep_);
}

template <typename bit_t>
bool IndexingSymmetricNp<bit_t>::operator!=(
    IndexingSymmetricNp<bit_t> const &rhs) const {
  return !operator==(rhs);
}   

template class IndexingSymmetricNp<uint16_t>;
template class IndexingSymmetricNp<uint32_t>;
template class IndexingSymmetricNp<uint64_t>;

} // namespace hydra::indexing::fermion
