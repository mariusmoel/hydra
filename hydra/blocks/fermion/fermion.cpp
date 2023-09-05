#include "fermion.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/utils/logger.h>

namespace hydra {

using namespace indexing;

Fermion::Fermion(int n_sites)
    : n_sites_(n_sites), //charge_conserved_(false), charge_(undefined_qn),
      //sz_conserved_(false), sz_(undefined_qn), 
      n_(undefined_qn), symmetric_(false), permutation_group_(), irrep_() {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Fermion: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingNoNp<uint16_t>(n_sites));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingNoNp<uint32_t>(n_sites));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingNoNp<uint64_t>(n_sites));
  } else {
    Log.err("Error creating Fermion: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);
  assert(n_sites >= 0);
}

Fermion::Fermion(int n_sites, int n)
    : n_sites_(n_sites), np_conserved_(true), //charge_conserved_(true), charge_(nup + ndn),
      //sz_conserved_(true), sz_(nup - ndn), 
      n_(n), symmetric_(false), permutation_group_(), irrep_() {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Fermion: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingNp<uint16_t>(n_sites, n));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingNp<uint32_t>(n_sites, n ));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingNp<uint64_t>(n_sites, n ));
  } else {
    Log.err("Error creating Fermion: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);
  utils::check_n_fermion(n_sites, n , "Fermion");
}

Fermion::Fermion(int n_sites, PermutationGroup group, Representation irrep)
    : n_sites_(n_sites), //charge_conserved_(false), charge_(undefined_qn),
      //sz_conserved_(false), sz_(undefined_qn),
      n_(undefined_qn), symmetric_(true),
      permutation_group_(allowed_subgroup(group, irrep)), irrep_(irrep) {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Fermion: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingSymmetricNoNp<uint16_t>(n_sites, group, irrep));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingSymmetricNoNp<uint32_t>(n_sites, group, irrep));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingSymmetricNoNp<uint64_t>(n_sites, group, irrep));
  } else {
    Log.err("Error creating Fermion: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);

  utils::check_n_sites(n_sites, group);
}

Fermion::Fermion(int n_sites, int n, PermutationGroup group,
                   Representation irrep)
    : n_sites_(n_sites), np_conserved_(true), //charge_conserved_(true), charge_(nup + ndn),
      //sz_conserved_(true), sz_(nup - ndn), n_up_(nup), 
      n_(n), symmetric_(true), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Fermion: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingSymmetricNp<uint16_t>(n_sites, n , group,
                                                irrep));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingSymmetricNp<uint32_t>(n_sites, n , group,
                                                irrep));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<FermionIndexing>(
        fermion::IndexingSymmetricNp<uint64_t>(n_sites, n , group,
                                                irrep));
  } else {
    Log.err("Error creating Fermion: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);

  utils::check_n_fermion(n_sites, n , "Fermion");
  utils::check_n_sites(n_sites, group);
}

bool Fermion::operator==(Fermion const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         //(charge_conserved_ == rhs.charge_conserved_) &&
         //(charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         //(sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (n_ == rhs.n_) && (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
bool Fermion::operator!=(Fermion const &rhs) const {
  return !operator==(rhs);
}

indexing::FermionIndexing const &Fermion::indexing() const {
  return *indexing_;
}

} // namespace hydra
