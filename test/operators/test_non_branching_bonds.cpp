#include "../catch.hpp"

#include <iostream>
#include <hydra/utils/close.h>
#include <hydra/operators/non_branching_bonds.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/algebra/algebra.h>
#include <hydra/algebra/matrix.h>
#include <hydra/algorithms/sparse_diag.h>

TEST_CASE("non_branching_bonds", "[operators]") {
  using namespace hydra;
  using namespace arma;
  Log("testing non_branching_bonds");

  cx_mat sx(mat({{0., 0.5}, {0.5, 0.}}), mat({{0., 0.}, {0., 0.}}));
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
  cx_mat sz(mat({{0.5, 0.0}, {0.0, -0.5}}), mat({{0., 0.}, {0., 0.0}}));

  cx_mat sp(mat({{0.0, 1.0}, {0.0, 0.0}}), mat({{0., 0.}, {0., 0.0}}));
  cx_mat sm(mat({{0.0, 0.0}, {1.0, 0.0}}), mat({{0., 0.}, {0., 0.0}}));
  cx_mat ones(mat({{1.0, 1.0}, {1.0, 1.0}}), mat({{1.0, 1.0}, {1.0, 1.0}}));

  for (auto ss : {ones}) {
    auto bond = Bond(ss, 0);
    auto block = Spinhalf(1);
    auto h = matrix(bond, block);
    REQUIRE(close(h, ss));
  }

  // Check ground state energy of TFI model
  int N = 14;
  double J = 1.0;
  double H = 1.0;
  BondList bonds;

  for (int i = 0; i < N - 1; ++i) {
    bonds << Bond("ISING", J, {i, (i + 1) % N});
  }
  for (int i = 0; i < N; ++i) {
    bonds << Bond(sx, H, i);
  }

  auto block = Spinhalf(N);
  double e = eig0(bonds, block);
  double e_dmrg = -7.411918598647893;
  REQUIRE(std::abs(e - e_dmrg) < 1e-8);

  // Check whether random bonds are created correctly
  for (int r = 0; r < 5; ++r) {

    for (int k = 1; k <= 6; ++k) {
      auto block = Spinhalf(k);
      int p2 = pow(2, k);

      auto mr = mat(p2, p2, fill::randn);
      std::vector<int> sites(k);
      std::iota(sites.begin(), sites.end(), 0);

      auto bondr = Bond(mr, sites);
      auto hr = matrix_real(bondr, block);
      REQUIRE(close(hr, mr));

      auto mc = cx_mat(p2, p2, fill::randn);
      auto bondc = Bond(mc, sites);
      auto hc = matrix_cplx(bondc, block);
      REQUIRE(close(hc, mc));
    }
  }

  // compare scalar chirality
  {
    auto block6 = Spinhalf(6);

    BondList bonds1;
    bonds1 << Bond("SCALARCHIRALITY", "Jchi", {0, 1, 2});
    bonds1 << Bond("SCALARCHIRALITY", "Jchi", {1, 2, 3});
    bonds1 << Bond("SCALARCHIRALITY", "Jchi", {2, 3, 4});
    bonds1 << Bond("SCALARCHIRALITY", "Jchi", {3, 4, 5});
    bonds1["Jchi"] = 1.0;
    auto H1 = matrix(bonds1, block6);

    cx_mat jchimat = kron(sx, kron(sy, sz) - kron(sz, sy)) +
                     kron(sy, kron(sz, sx) - kron(sx, sz)) +
                     kron(sz, kron(sx, sy) - kron(sy, sx));

    BondList bonds2;
    bonds2 << Bond(jchimat, "Jchi", {0, 1, 2});
    bonds2 << Bond(jchimat, "Jchi", {1, 2, 3});
    bonds2 << Bond(jchimat, "Jchi", {2, 3, 4});
    bonds2 << Bond(jchimat, "Jchi", {3, 4, 5});
    bonds2["Jchi"] = 1.0;
    auto H2 = matrix(bonds2, block6);

    REQUIRE(close(H1, H2));
  }

  // compare scalar chirality (12 site kagome)
  auto block12 = Spinhalf(12);

  BondList bonds1;
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {0, 4, 6});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {3, 1, 9});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {9, 7, 4});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {4, 2, 10});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {10, 8, 5});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {6, 10, 1});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {1, 5, 7});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {7, 11, 2});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {2, 3, 8});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {8, 9, 0});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {5, 0, 11});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {11, 6, 3});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {4, 10, 6});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {1, 7, 9});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {7, 2, 4});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {2, 8, 10});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {8, 0, 5});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {10, 5, 1});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {5, 11, 7});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {11, 3, 2});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {3, 9, 8});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {9, 4, 0});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {0, 6, 11});
  bonds1 << Bond("SCALARCHIRALITY", "Jchi", {6, 1, 3});
  bonds1["Jchi"] = 1.0;
  auto H1 = matrix(bonds1, block12);

  cx_mat jchimat = kron(kron(sx, sy), sz) + kron(kron(sy, sz), sx) +
                   kron(kron(sz, sx), sy) - kron(kron(sx, sz), sy) -
                   kron(kron(sy, sx), sz) - kron(kron(sz, sy), sx);

  BondList bonds2;
  bonds2 << Bond(jchimat, "Jchi", {0, 4, 6});
  bonds2 << Bond(jchimat, "Jchi", {3, 1, 9});
  bonds2 << Bond(jchimat, "Jchi", {9, 7, 4});
  bonds2 << Bond(jchimat, "Jchi", {4, 2, 10});
  bonds2 << Bond(jchimat, "Jchi", {10, 8, 5});
  bonds2 << Bond(jchimat, "Jchi", {6, 10, 1});
  bonds2 << Bond(jchimat, "Jchi", {1, 5, 7});
  bonds2 << Bond(jchimat, "Jchi", {7, 11, 2});
  bonds2 << Bond(jchimat, "Jchi", {2, 3, 8});
  bonds2 << Bond(jchimat, "Jchi", {8, 9, 0});
  bonds2 << Bond(jchimat, "Jchi", {5, 0, 11});
  bonds2 << Bond(jchimat, "Jchi", {11, 6, 3});
  bonds2 << Bond(jchimat, "Jchi", {4, 10, 6});
  bonds2 << Bond(jchimat, "Jchi", {1, 7, 9});
  bonds2 << Bond(jchimat, "Jchi", {7, 2, 4});
  bonds2 << Bond(jchimat, "Jchi", {2, 8, 10});
  bonds2 << Bond(jchimat, "Jchi", {8, 0, 5});
  bonds2 << Bond(jchimat, "Jchi", {10, 5, 1});
  bonds2 << Bond(jchimat, "Jchi", {5, 11, 7});
  bonds2 << Bond(jchimat, "Jchi", {11, 3, 2});
  bonds2 << Bond(jchimat, "Jchi", {3, 9, 8});
  bonds2 << Bond(jchimat, "Jchi", {9, 4, 0});
  bonds2 << Bond(jchimat, "Jchi", {0, 6, 11});
  bonds2 << Bond(jchimat, "Jchi", {6, 1, 3});
  bonds2["Jchi"] = 1.0;
  auto H2 = matrix(bonds2, block12);

  // HydraPrint(norm(H1));
  // HydraPrint(norm(H2));

  REQUIRE(close(H1, H2));
}
