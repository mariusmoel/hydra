#pragma once

#ifdef HYDRA_ENABLE_MPI
#include <mpi.h>
#endif

#include <hydra/operators/bondlist.h>


namespace hydra::testcases::fermion {

BondList Chain(int n_sites, double J1, double J2 = 0);
std::tuple<BondList, arma::vec> Chain_fullspectrum_nup(int n_sites, int nup);
BondList HB_alltoall(int n_sites);
std::tuple<BondList, double> triangular_12_complex(int nup, double eta);

} // namespace hydra::testcases::fermion
