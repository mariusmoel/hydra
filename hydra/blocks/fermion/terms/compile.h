#pragma once

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra::fermion {

const std::vector<std::string> special_bond_types = {
    "HOP"}; //,  "NUMBER"}; // "C", "CDAG",

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace hydra::fermion
