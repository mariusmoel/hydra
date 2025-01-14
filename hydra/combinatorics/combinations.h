#pragma once

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/bit_patterns.h>
#include <hydra/common.h>

namespace hydra::combinatorics {

template <typename bit_t> class CombinationsIterator;

// Combinations
template <typename bit_t = std_bit_t> class Combinations {
public:
  using iterator_t = CombinationsIterator<bit_t>;

  Combinations() = default;
  Combinations(int n, int k);

  inline int n() const { return n_; }
  inline int k() const { return k_; }
  inline idx_t size() const { return size_; };
  inline iterator_t begin() const { return begin_; }
  inline iterator_t end() const { return end_; }

private:
  int n_, k_;
  idx_t size_;
  iterator_t begin_, end_;
};

// CombinationsIterator
template <typename bit_t> class CombinationsIterator {
public:
  CombinationsIterator() = default;
  CombinationsIterator(int n, int k, idx_t idx);

  inline bool operator==(CombinationsIterator<bit_t> const &rhs) const {
    return idx_ == rhs.idx_;
  }
  inline bool operator!=(CombinationsIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }
  inline CombinationsIterator &operator++() {
    current_ = combinatorics::get_next_pattern(current_);
    ++idx_;
    return *this;
  }
  inline bit_t operator*() const { return current_; }

private:
  bit_t current_;
  idx_t idx_;
};

#ifdef _OPENMP
// CombinationsThread
template <typename bit_t = std_bit_t> class CombinationsThread {
public:
  using iterator_t = CombinationsIterator<bit_t>;

  CombinationsThread() = default;
  CombinationsThread(int n, int k);

  int n() const { return n_; }
  int k() const { return k_; }
  idx_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int n_, k_;
  idx_t size_;
  iterator_t begin_, end_;
};

template <typename bit_t>
inline CombinationsThread<bit_t> ThreadStates(Combinations<bit_t> const &si) {
  return CombinationsThread<bit_t>(si.n(), si.k());
}

#endif

} // namespace hydra::combinatorics
