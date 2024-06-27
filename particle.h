#pragma once

#include "falcon.h"
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
    };
    
    struct particle_system {
        size_t n;
        struct particle *part;
        struct particle *last;
    };

    //struct particle_system mainsys;
    //struct sys copysys;
    

    /*void printSystemEnergyDiagnostics (struct particle_system s) {
            fflush(stdout);
            real_t kinetic = falcon::diagnostics::system_kinetic_energy(s);
            real_t pot = falcon::diagnostics::system_potential_energy(s);
            real_t e0 = -pot+kinetic;
            std::cout<<"Energy: "<<e0<<std::endl;
    }*/
    /*void printSystemCOMDiagnostics (struct particle_system s) {
        long double cmpos[3], cmvel[3];
        falcon::diagnostics::system_center_of_mass(mainsys, cmpos, cmvel);
        printf("com=%18.12Lg %18.12Lg %18.12Lg\n",cmpos[0], cmpos[1], cmpos[2]);
        fflush(stdout);

    }*/

    /*void readFile(struct particle_system* system, const size_t restart_val) {
    if(restart_val == 0)
    {
        snapnum = 0;
        t_now = 0;

        bool gadget_file = false;
        std::vector < double >array;

        //	  start("reading bodies");
        size_t number_of_lines = 0;
        std::string line;
        std::ifstream file(input_fname);
        while(std::getline(file, line))
            ++number_of_lines;
        numBodies = number_of_lines;
        file.close();

        std::ifstream fin(input_fname);
        array.resize(numBodies * 7);
        array.assign(std::istream_iterator < double >(fin), std::istream_iterator < double >());
        fin.close();
        //	  stop("reading bodies");

          
        

        mainsys.n = numBodies;
        mainsys.part = (struct particle *) malloc(numBodies * sizeof(struct particle));
        mainsys.last = &mainsys.part[numBodies - 1];

        for(size_t b = 0; b < numBodies; b++)
        {
            mainsys.part[b].id = b;
            mainsys.part[b].mass = array[b * 7+0];
            mainsys.part[b].pos[0] = array[b * 7+1];
            mainsys.part[b].pos[1] = array[b * 7+2];
            mainsys.part[b].pos[2] = array[b * 7+3];
            mainsys.part[b].vel[0] = array[b * 7+4];
            mainsys.part[b].vel[1] = array[b * 7+5];
            mainsys.part[b].vel[2] = array[b * 7+6];
            // for pn terms
            mainsys.part[b].acc_pn[0] = 0.0;
            mainsys.part[b].acc_pn[1] = 0.0;
            mainsys.part[b].acc_pn[2] = 0.0;
            mainsys.part[b].w[0] = mainsys.part[b].vel[0];
            mainsys.part[b].w[1] = mainsys.part[b].vel[1];
            mainsys.part[b].w[2] = mainsys.part[b].vel[2];
        }
    
        
    }
    else
    {
        falcon::io::read_hdf5_snapshot(snapnum);
        for(size_t b=0; b<mainsys.n;b++) {
        //for pn terms
            mainsys.part[b].acc_pn[0] = 0.0;
            mainsys.part[b].acc_pn[1] = 0.0;
            mainsys.part[b].acc_pn[2] = 0.0;
            mainsys.part[b].w[0] = mainsys.part[b].vel[0];
            mainsys.part[b].w[1] = mainsys.part[b].vel[1];
            mainsys.part[b].w[2] = mainsys.part[b].vel[2];

        }

    }

    }*/

}
