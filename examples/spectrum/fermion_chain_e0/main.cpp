#include <hydra/all.h>

int main() {
  using namespace hydra;

  int n_sites = 4;
  int nup = 2;

  // Define the Hilbert space block
  auto block = Fermion(n_sites,nup);

  // Define the nearest-neighbor Hopping Hamiltonian
  BondList bonds;
  for (int i = 0; i < n_sites; ++i) {
    bonds << Bond("HOP", "J", {i, (i + 1) % n_sites});
  }

  // Set the coupling constant "J" to one
  bonds["J"] = -2.234;
  // Compute and print the ground state energy
  double e0 = eig0(bonds, block);
  HydraPrint(e0);
  
  return EXIT_SUCCESS;
}
