#include "Laser.hpp"

namespace FHHG::Laser {
    /** 
     * converts e E_0 / (hbar omega_L) to 1 / pm
     * if E_0 is given in MV / cm and (hbar omega_L) in meV
     */
    constexpr h_float field_conversion = 1e-1; 

    Laser::Laser(h_float photon_energy, h_float E_0, int max_frequency, int n_subdivisions)
        : momentum_amplitude{field_conversion * E_0 / (photon_energy)}, 
          max_frequency{max_frequency},
          n_subdivisions{n_subdivisions} 
    {}
}