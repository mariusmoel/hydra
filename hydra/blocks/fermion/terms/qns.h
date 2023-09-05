#include <map>
#include <string>

#include <hydra/operators/bondlist.h>

namespace hydra::fermion {

const std::map<std::string, int> special_bonds_nup = {
    {"HOP", 0}, {"NUMBER", 0}}; //{"C", -1}, {"CDAG", +1},
    //{"SZ", 0}, {"S+", 1},         {"S-", -1},      {"SCALARCHIRALITY", 0}};

int nup(BondList bonds, double precision = 1e-12);

} // namespace hydra::fermion
