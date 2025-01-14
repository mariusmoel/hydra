#include "non_branching_bonds.h"

#include <tuple>
#include <vector>

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

namespace hydra::operators {

BondList non_branching_bonds(Bond const &bond, double precision) {
  if (bond.type_defined()) {
    return BondList({bond});
  }

  arma::cx_mat mat = bond.matrix();
  int m = (int)mat.n_rows;
  int n = (int)mat.n_cols;
  if (m != n) {
    Log.err("Error: bond matrix is not square");
  }

  std::vector<std::tuple<int, int, complex>> all_entries;

  // Get diagonal elements
  for (int i = 0; i < n; ++i) {
    if (std::abs(mat(i, i)) > precision) {
      all_entries.push_back({i, i, mat(i, i)});
    }
  }

  // Get offidagonal elements
  for (int n_diag = 1; n_diag < n; ++n_diag) {
    for (int i = 0; i < n - n_diag; ++i) {
      if (std::abs(mat(i, i + n_diag)) > precision) {
        all_entries.push_back({i, i + n_diag, mat(i, i + n_diag)});
      }
      if (std::abs(mat(i + n_diag, i)) > precision) {
        all_entries.push_back({i + n_diag, i, mat(i + n_diag, i)});
      }
    }
  }

  // Reduce to minimal number of non-branching terms
  std::vector<arma::cx_mat> mats_nb;

  while (all_entries.size() != 0) {
    std::vector<int> forbidden_columns;
    std::vector<int> forbidden_rows;
    std::vector<std::tuple<int, int, complex>> current_entries;
    std::vector<int> delete_entries;
    int i = 0;
    for (auto [row, column, coeff] : all_entries) {

      if ((std::find(forbidden_rows.begin(), forbidden_rows.end(), row) ==
           forbidden_rows.end()) &&
          (std::find(forbidden_columns.begin(), forbidden_columns.end(),
                     column) == forbidden_columns.end())) {
        current_entries.push_back({row, column, coeff});
        forbidden_rows.push_back(row);
        forbidden_columns.push_back(column);
        delete_entries.push_back(i);
      }
      ++i;
    }

    for (int i = delete_entries.size() - 1; i >= 0; --i)
      all_entries.erase(all_entries.begin() + delete_entries[i]);

    // Create non-branching matrix
    arma::cx_mat mat_nb(m, n, arma::fill::zeros);
    for (auto [i, j, coeff] : current_entries) {
      mat_nb(i, j) = coeff;
    }
    mats_nb.push_back(mat_nb);
  }

  BondList bonds_nb;
  for (arma::cx_mat mat_nb : mats_nb) {
    if (bond.coupling_defined()) {
      bonds_nb << Bond(mat_nb, bond.coupling(), bond.sites());
    } else {
      bonds_nb << Bond(mat_nb, bond.coupling_name(), bond.sites());
    }
  }
  return bonds_nb;
}

BondList non_branching_bonds(BondList const &bonds, double precision) {
  BondList bonds_nb;
  for (Bond bond : bonds) {
    bonds_nb = bonds_nb + non_branching_bonds(bond, precision);
  }
  return bonds_nb;
}

bool is_non_branching_bond(Bond const &bond, double precision) {
  if (bond.type_defined()) {
    return false;
  }
  arma::cx_mat mat = bond.matrix();
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    int non_zero_in_row = 0;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      if (std::abs(mat(i, j)) > precision) {
        ++non_zero_in_row;
      }
    }
    if (non_zero_in_row > 1) {
      return false;
    }
  }
  return true;
}

template <typename bit_t, typename coeff_t>
NonBranchingBond<bit_t, coeff_t>::NonBranchingBond(Bond const &bond,
                                                   double precision)
    : sites_(bond.sites()), dim_(1 << sites_.size()), mask_(0) {
  if (!is_non_branching_bond(bond, precision)) {
    Log.err("Error: trying to create a NonBranchingBond from a Bond which is "
            "branching");
  }

  if (!bond.coupling_defined()) {
    Log.err("Error: cannot create Nonbranching bond from Bond without having "
            "its coupling defined.");
  }
  coeff_t cpl = bond.coupling<coeff_t>();

  for (auto s : bond.sites()) {
    mask_ |= ((bit_t)1 << s);
  }
  // int n_sites = 2;
  // Log("X: {}", BSTR(mask_));

  mask_ = ~mask_;
  // Log("Y: {}", BSTR(mask_));

  arma::cx_mat matrix_ = bond.matrix();

  // Matrix dimension is 2**(no. sites of bond)
  if ((matrix_.n_cols != dim_) || (matrix_.n_rows != dim_)) {
    Log.err("Error: invalid matrix dimension for non-branching bond matrix.");
  }

  non_zero_term_ = std::vector<bool>(dim_, false);
  state_applied_ = std::vector<bit_t>(dim_, 0);
  coeff_ = std::vector<coeff_t>(dim_, 0.);

  for (bit_t in = 0; in < dim_; ++in) {
    // int non_zero_in_row = 0;

    for (bit_t out = 0; out < dim_; ++out) {
      if (std::abs(matrix_(out, in)) > precision) {
        non_zero_term_[in] = true;
        state_applied_[in] = out;
        if constexpr (is_real<coeff_t>()) {
          if (std::abs(imag(matrix_(out, in))) > precision) {
            Log.err("Error: trying to create a real NonBranchingBond, but "
                    "found a truly complex matrix entry");
          }
          coeff_[in] = cpl * real(matrix_(out, in));
        } else {
          coeff_[in] = cpl * matrix_(out, in);
        }
        // ++non_zero_in_row;
      }
    }

    // // security check
    // if (non_zero_term_[in]) {
    //   assert(non_zero_in_row == 1);
    // } else {
    //   assert(non_zero_in_row == 0);
    // }
  }

  // int n_sites = 1;
  // for (bit_t in = 0; in < dim_; ++in) {
  //   std::cout << "a " << BSTR(in) << " -> " << BSTR(state_applied_[in]) << " "
  //             << coeff_[in] << " " << non_zero_term_[in] << std::endl;
  // }
}

template <typename bit_t, typename coeff_t>
bool NonBranchingBond<bit_t, coeff_t>::is_diagonal() const {
  for (bit_t i = 0; i < dim_; ++i) {
    if ((non_zero_term_[i]) && (state_applied_[i] != i)) {
      return false;
    }
  }
  return true;
}

template <typename bit_t, typename coeff_t>
bool NonBranchingBond<bit_t, coeff_t>::non_zero_term(bit_t local_state) const {
  return non_zero_term_[local_state];
}
template <typename bit_t, typename coeff_t>
coeff_t NonBranchingBond<bit_t, coeff_t>::coeff(bit_t local_state) const {
  return coeff_[local_state];
}

template <typename bit_t, typename coeff_t>
std::pair<bit_t, coeff_t>
NonBranchingBond<bit_t, coeff_t>::state_coeff(bit_t local_state) const {
  return {state_applied_[local_state], coeff_[local_state]};
}

template <typename bit_t, typename coeff_t>
bit_t NonBranchingBond<bit_t, coeff_t>::extract_local_state(bit_t state) const {
  bit_t local_state = 0;
  for (int i = 0; i < (int)sites_.size(); ++i) {
    local_state |= bitops::gbit(state, sites_[i]) << i;
  }
  return local_state;
}

template <typename bit_t, typename coeff_t>
bit_t NonBranchingBond<bit_t, coeff_t>::deposit_local_state(bit_t local_state,
                                                            bit_t state) const {
  // int n_sites = 2;
  // Log("a: {}", BSTR(state));
  // Log("mask: {}", BSTR(mask_));

  state &= mask_; // clear bits on site
  // Log("b: {}", BSTR(state));
  for (int i = 0; i < (int)sites_.size(); ++i) {
    state |= bitops::gbit(local_state, i) << sites_[i];
    // Log("c: {}", BSTR(state));
  }
  return state;
}

template <typename bit_t, typename coeff_t>
int NonBranchingBond<bit_t, coeff_t>::number_difference() const {
  int diff = 0;
  bool first_diff = true;
  for (bit_t state = 0; state < dim_; ++state) {
    if (non_zero_term_[state]) {
      int diff_state =
          bitops::popcnt(state_applied_[state]) - bitops::popcnt(state);
      if (first_diff) {
        diff = diff_state;
        first_diff = false;
      } else {
        if (diff_state != diff) {
          return undefined_qn;
        }
      }
    }
  }
  return diff;
}

template class NonBranchingBond<uint16_t, double>;
template class NonBranchingBond<uint32_t, double>;
template class NonBranchingBond<uint64_t, double>;
template class NonBranchingBond<uint16_t, complex>;
template class NonBranchingBond<uint32_t, complex>;
template class NonBranchingBond<uint64_t, complex>;

} // namespace hydra::operators
