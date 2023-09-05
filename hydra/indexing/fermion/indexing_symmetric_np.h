#pragma once

#include <extern/gsl/span>
#include <utility>
#include <vector>

#include <hydra/common.h>
#include <hydra/indexing/fermi_table.h>
#include <hydra/indexing/combinations_indexing.h>
#include <hydra/indexing/fermion/symmetric_iterator.h>
#include <hydra/symmetries/operations/symmetry_operations.h>

#include <hydra/symmetries/group_action/group_action_lookup.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::indexing::fermion {

template <typename bit_t> class IndexingSymmetricNp {
public:
  using iterator_t = SymmetricIterator<bit_t>;
  using span_size_t = gsl::span<int const>::size_type;

  IndexingSymmetricNp() = default;
  IndexingSymmetricNp(int n_sites, int n, PermutationGroup permutation_group,
                      Representation irrep);
  inline int n_sites() const { return n_sites_; }
  inline int n() const { return n_; }
  GroupActionLookup<bit_t> const &group_action() const { return group_action_; }
  Representation const &irrep() const { return irrep_; }
  inline idx_t size() const { return size_; }

  inline bit_t state(idx_t idx) const { return reps_[idx]; }
  inline double norm(idx_t idx) const { return norms_[idx]; }

  inline idx_t index(bit_t state) const {
    return index_for_rep_[combinations_indexing_.index(state)];
  }
  inline bit_t representative(bit_t raw_state) const {
    return reps_[index(raw_state)];
  }
  inline std::pair<idx_t, int> index_sym(bit_t raw_state) const {
    idx_t raw_idx = combinations_indexing_.index(raw_state);
    idx_t index = index_for_rep_[raw_idx];
    if (index == invalid_index) {
      return {invalid_index, 0};
    }
    idx_t start = sym_limits_for_rep_[raw_idx].first;
    return {index, syms_[start]};
  }
  inline std::pair<idx_t, gsl::span<int const>>
  index_syms(bit_t raw_state) const {
    idx_t raw_idx = combinations_indexing_.index(raw_state);
    idx_t index = index_for_rep_[raw_idx];
    auto [start, length] = sym_limits_for_rep_[raw_idx];
    return {index, {syms_.data() + start, length}};
  }

// -----------------------------------------------------------------------------------

  std::pair<idx_t, bool> index_fermi(bit_t fermions, int sym) const {
    bit_t fermions_rep = group_action_.apply(sym, fermions);
    idx_t idx_fermions_rep = combinations_indexing_.index(fermions_rep);
    bool fermi_fermions = fermi_table_fermions_.sign(sym, fermions);
    return {idx_fermions_rep, fermi_fermions};
  }

  std::pair<idx_t, bool> index_fermions_fermi(bit_t fermions, int sym,
                                         bit_t fermimask) const {
    bit_t fermions_rep = group_action_.apply(sym, fermions);
    idx_t idx_fermions_rep = combinations_indexing_.index(fermions_rep);
    bool fermi_fermions = (bitops::popcnt(fermions & fermimask) & 1);
    fermi_fermions ^= fermi_table_fermions_.sign(sym, fermions);
    return {idx_fermions_rep, fermi_fermions};
  }

  // Fermi sign when applying sym on states
  bool fermi_bool_fermions(int sym, bit_t ups) const {
    return fermi_table_fermions_.sign(sym, ups);
  }


// -----------------------------------------------------------------------------------

  inline iterator_t begin() const { return begin_; }
  inline iterator_t end() const { return end_; }

  bool operator==(IndexingSymmetricNp const &rhs) const;
  bool operator!=(IndexingSymmetricNp const &rhs) const;
  
private:
  int n_sites_;
  int n_;
  GroupActionLookup<bit_t> group_action_;
  Representation irrep_;
  CombinationsIndexing<bit_t> combinations_indexing_;

  FermiTableCombinations<bit_t> fermi_table_fermions_;

  std::vector<bit_t> reps_;
  std::vector<idx_t> index_for_rep_;
  std::vector<int> syms_;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits_for_rep_;
  std::vector<double> norms_;

  idx_t size_;
  iterator_t begin_, end_;
};


} // namespace hydra::indexing::fermion
