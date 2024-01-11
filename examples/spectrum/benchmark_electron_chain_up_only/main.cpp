#include <hydra/all.h>

int main() {
  using namespace hydra;

  int n_sites = 12;
  int nup = 4;
  int ndn = 0;

  // Define the Hilbert space block
  auto block = Electron(n_sites,nup, ndn);

  // Define the nearest-neighbor Hopping Hamiltonian
  BondList bonds;
  for (int i = 0; i < n_sites; ++i) {
    bonds << Bond("HOP", "J", {i, (i + 1) % n_sites});
    //bonds << Bond("HOP", "J", {(i+1)% n_sites, i});
  }

  // Set the coupling constant "J" to one
  bonds["J"] = -2.234;
  // Compute and print the ground state energy
  double e0 = eig0(bonds, block);
  HydraPrint(e0);
  
  return EXIT_SUCCESS;
}
