#include "../../catch.hpp"

#include <iostream>

#include "testcases_fermion.h"
#include <hydra/blocks/fermion/fermion_apply.h>
#include <hydra/algebra/algebra.h>
#include <hydra/algebra/matrix.h>
#include <hydra/algorithms/sparse_diag.h>
#include <hydra/utils/close.h>

using namespace hydra;

void test_apply(BondList bonds) {
  int N = bonds.n_sites();
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Fermion(N, nup);
    auto H = matrix_real(bonds, block, block);
    REQUIRE(H.is_hermitian(1e-8));

    arma::vec v(block.size(), arma::fill::randn);
    arma::vec w1 = H * v;
    arma::vec w2(block.size(), arma::fill::zeros);
    apply(bonds, block, v, block, w2);
    REQUIRE(close(w1, w2));

    arma::vec evals_mat;
    arma::eig_sym(evals_mat, H);

    double e0_mat = evals_mat(0);
    double e0_app = eig0_real(bonds, block);
    REQUIRE(close(e0_mat, e0_app));
  }
}

TEST_CASE("fermion_apply", "[fermion]") {
  using namespace hydra::testcases::fermion;

  Log.out("fermion_apply: Chain apply test, J=1.0, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto bonds = Chain(N, 1.0);
    test_apply(bonds);
  }

  Log.out("fermion_apply: Chain alltoall apply test, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto bonds = HB_alltoall(N);
    test_apply(bonds);
  }

  Log.out("fermion_apply: Chain all-to-all Sz <-> NoSz comparison");
  for (int n_sites = 2; n_sites <= 6; ++n_sites) {
    auto bonds = HB_alltoall(n_sites);
    auto block_no_sz = Fermion(n_sites);
    auto e0_no_sz = eig0_real(bonds, block_no_sz);
    auto e0s_sz = std::vector<double>();

    for (int nup = 0; nup <= n_sites; ++nup) {
      auto block_sz = Fermion(n_sites, nup);
      auto e0_sz = eig0_real(bonds, block_sz);
      e0s_sz.push_back(e0_sz);
    }
    auto e0_sz = *std::min_element(e0s_sz.begin(), e0s_sz.end());
    REQUIRE(close(e0_sz, e0_no_sz));
  }

  {
    Log("fermion_matrix: Triangular J1J2Jchi N=12");
    std::string lfile =
        HYDRA_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.lat";

    auto bondlist = read_bondlist(lfile);
    bondlist["J1"] = 1.00;
    bondlist["J2"] = 0.15;
    bondlist["Jchi"] = 0.09;

    int n_sites = 12;
    int n_up = 6;
    auto fermion = Fermion(n_sites, n_up);
    auto e0 = eig0_cplx(bondlist, fermion);
    double energy = -6.9456000700824329641;

    // Log("{:.18f} {:.18f}", e0, energy);

    REQUIRE(close(e0, energy));
  }
}
