#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/combinations_index.h>

#include <hydra/indexing/lin_table.h>

namespace hydra::indexing::fermion {
template <typename bit_t> class IndexingNp {
public:
  using iterator_t = combinatorics::CombinationsIndexIterator<bit_t>;

  IndexingNp() = default;
  IndexingNp(int n_sites, int n);

  inline int n_sites() const { return n_sites_; }
  inline int n() const { return n_; }
  inline idx_t size() const { return size_; }
  inline idx_t index(bit_t fermions) const { return lintable_.index(fermions); }
  inline bit_t state(idx_t index) const { return states_[index]; }

  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

  bool operator==(IndexingNp const &rhs) const;
  bool operator!=(IndexingNp const &rhs) const;
  
private:
  int n_sites_;
  int n_;
  indexing::LinTable<bit_t> lintable_;
  std::vector<bit_t> states_;
  idx_t size_;
  iterator_t begin_, end_;
};

} // namespace hydra::indexing::fermion
