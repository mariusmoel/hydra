#include "compile.h"
#include <hydra/operators/compiler.h>
#include <hydra/operators/non_branching_bonds.h>
#include <hydra/utils/logger.h>
#include <hydra/utils/print_macro.h>
#include <hydra/utils/timing.h>

namespace hydra::fermion {

BondList compile(BondList const &bonds, double precision) {

  BondList bonds_explicit =
      operators::compile_explicit(bonds, precision, "keep");
  BondList bonds_special;
  BondList bonds_generic;
  
  for (auto bond : bonds_explicit) {
    if (bond.type_defined()) {
      std::string type = bond.type();
      if (std::find(special_bond_types.begin(), special_bond_types.end(),
                    type) == special_bond_types.end()) {
        Log.err("Error compiling BondList: invalid or undefined type found: {}",
                type);
      } else {

        if ((type == "HOP") || (type == "NUMBER")) {
          bonds_special << Bond("HOP", bond.coupling(), bond.sites());
          bonds_special << Bond("NUMBER", bond.coupling(), bond.sites());
        } else {
          bonds_special << bond;
        }
      }
    } else {
      BondList bonds_nb = operators::non_branching_bonds(bond, precision);
      for (auto b : bonds_nb){
	bonds_generic << b;
      }
    }
  }
  return bonds_special + bonds_generic;
}

} // namespace hydra::fermion
