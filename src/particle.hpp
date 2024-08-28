#pragma once

#include "falcon.hpp"
//#include "diagnostics.h"

namespace falcon::particles {
    struct particle {
        real_t pos[3];
        real_t mass;
        real_t vel[3];
        real_t w[3]; /*needed for Auxiliary Velocity Algorithm with PN terms. See Hellstrom and Mikkola (2010) for more details */
        real_t t_last;
        real_t timestep;
        real_t acc[3];
        real_t acc_pn[3];
        real_t jerk[3];
        size_t id;
        real_t pot;
        bool use_w;
    };
    
    struct particle_system {
        size_t n;
        struct particle *part;
        struct particle *last;
    };

    


}
