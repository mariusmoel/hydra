#include "hashes.h"
#include <complex>
#include <hydra/random/hash_functions.h>

namespace hydra::random {

uint64_t hash(Permutation const &perm) {
  uint64_t h = 0;
  for (int i = 0; i < perm.size(); ++i) {
    h = hash_combine(h, hash_fnv1((uint64_t)perm[i]));
  }
  return h;
}

uint64_t hash(PermutationGroup const &group) {
  uint64_t h = 0;
  for (auto perm : group) {
    h = hash_combine(h, hash(perm));
  }
  return h;
}

uint64_t hash(Representation const &irrep) {
  uint64_t h = 0;
  for (auto sym : irrep.allowed_symmetries()) {
    h = hash_combine(h, hash_fnv1((uint64_t)sym));
  }
  for (auto ch : irrep.characters()) {
    double phase = std::arg(ch) + M_PI;
    double phase_scattered = std::log(
        std::abs(std::sin(1.7365529164217 * phase + 2.56381234623457)) +
        0.1234567);
    uint64_t angle = uint64_t(std::abs(phase_scattered * 1234567.0));
    h = hash_combine(h, hash_fnv1(angle));
  }
  return h;
}

uint64_t hash(Spinhalf const &spinhalf) {
  uint64_t h =
      spinhalf.n_sites() == 0 ? 0 : hash_fnv1((uint64_t)spinhalf.n_sites());
  if (spinhalf.sz_conserved()) {
    h = hash_combine(h, hash_fnv1((uint64_t)spinhalf.n_up()));
  }
  if (spinhalf.symmetric()) {
    h = hash_combine(h, hash(spinhalf.permutation_group()));
    h = hash_combine(h, hash(spinhalf.irrep()));
  }
  return h;
}

uint64_t hash(Fermion const &fermion) {
  uint64_t h =
      fermion.n_sites() == 0 ? 0 : hash_fnv1((uint64_t)fermion.n_sites());
  if (fermion.np_conserved()) {
    h = hash_combine(h, hash_fnv1((uint64_t)fermion.n()));
  }
  if (fermion.symmetric()) {
    h = hash_combine(h, hash(fermion.permutation_group()));
    h = hash_combine(h, hash(fermion.irrep()));
  }
  return h;
}

uint64_t hash(tJ const &tj) {
  uint64_t h = tj.n_sites() == 0 ? 0 : hash_fnv1((uint64_t)tj.n_sites());
  if (tj.charge_conserved() && tj.sz_conserved()) {
    h = hash_combine(h, hash_fnv1((uint64_t)tj.n_up()));
    h = hash_combine(h, hash_fnv1((uint64_t)tj.n_dn()));
  }
  if (tj.symmetric()) {
    h = hash_combine(h, hash(tj.permutation_group()));
    h = hash_combine(h, hash(tj.irrep()));
  }
  return h;
}

uint64_t hash(Electron const &electron) {
  uint64_t h =
      electron.n_sites() == 0 ? 0 : hash_fnv1((uint64_t)electron.n_sites());
  if (electron.charge_conserved() && electron.sz_conserved()) {
    h = hash_combine(h, hash_fnv1((uint64_t)electron.n_up()));
    h = hash_combine(h, hash_fnv1((uint64_t)electron.n_dn()));
  }
  if (electron.symmetric()) {
    h = hash_combine(h, hash(electron.permutation_group()));
    h = hash_combine(h, hash(electron.irrep()));
  }
  return h;
}

} // namespace hydra::random
