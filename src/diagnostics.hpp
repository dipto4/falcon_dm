/* --------------------------------------------------------------------------------
 * Filename: diagnostics.h
 * Purpose: Return useful diagnostics about the N-body system
 *
 * Main functions:
 * real_t system_kinetic_energy(const struct falcon::particles::particle_system *s)
 *
 * real_t system_potential_energy(struct falcon::particles::particle_system *s)
 *
 * void system_center_of_mass(const struct falcon::particles::particle_system *s, real_t *cmpos, real_t *cmvel)
 *
 * --------------------------------------------------------------------------------*/

#pragma once
#include <math.h>
#include <omp.h>
#include <iomanip>
#include "falcon.hpp"
#include "particle.hpp"

namespace falcon::diagnostics {

    /* --------------------------------------------------------------------------
     * Function: real_t system_kinetic_energy(const struct falcon::particles::particle_system *s)
     * 
     * Input: falcon::particles::particle_system * (pointer to particle system)
     * Output: Kinetic energy (of type real_t)
     * -------------------------------------------------------------------------- */


    real_t system_kinetic_energy(const struct falcon::particles::particle_system *s) {
        real_t kin_e = 0;
#pragma omp parallel for reduction(+:kin_e)
        for(size_t i=0;i<s->n;i++) {
            kin_e += 0.5 * s->part[i].mass * (s->part[i].vel[0] * s->part[i].vel[0] + s->part[i].vel[1] * s->part[i].vel[1] + s->part[i].vel[2] * s->part[i].vel[2]);
        }
        return kin_e;

    }
    /*--------------------------------------------------------------------------
     * Function: real_t system_potential_energy(struct falcon::particles::particle_system *s)
     *
     * Input: falcon::particles::particle_system * (pointer to particle system)
     * Output: Potential energy (of type real_t)
     ---------------------------------------------------------------------------*/

    real_t system_potential_energy(struct falcon::particles::particle_system *s) {
        real_t pot_e = 0;
        for(size_t i=0;i<s->n;i++) {
            real_t pot_i = 0;
            //#pragma omp parallel for if(s.part[i].id == 0) reduction(+:pot_i)
            for(size_t j=0;j<s->n;j++) {
                if((s->part[i].id == s->part[j].id) || (s->part[j].id >= 1 && s->part[i].id >= 1))
                    continue;

                real_t dr[3];
                for(size_t d=0;d<3;d++)
                    dr[d] = s->part[i].pos[d] - s->part[j].pos[d];

                real_t r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
                //printf("part %i %21.16f\n",s.part[i].id,r);
                pot_i += s->part[j].mass / r;
            }
            s->part[i].pot = pot_i;
        }
#pragma omp parallel for reduction(+:pot_e)
        for(size_t i=0; i<s->n;i++) {

            pot_e += 0.5 * s->part[i].mass * s->part[i].pot;

        }

        return pot_e;

    }
    
    /*-------------------------------------------------------------------
     * Function: system_center_of_mass(const struct falcon::particles::particle_system *s, real_t *cmpos, real_t *cmvel)
     * 
     * Inputs: falcon::particles::particle_system* (pointer to particle system), real_t * (center of mass position array), real_t *(center of mass velocity array)
     * Output: void (center of mass position and vel are written to input arrays)
     ------------------------------------------------------------------*/

    void system_center_of_mass(const struct falcon::particles::particle_system *s, real_t *cmpos, real_t *cmvel) {
        real_t mass = 0., pos[3] = { 0., 0., 0. }, vel[3] = {0., 0., 0.};
#pragma omp parallel for reduction(+:pos,vel,mass)
        for(size_t p = 0; p < s->n; p++)
        {
            for(int i = 0; i < 3; i++)
            {
                pos[i] += s->part[p].mass*s->part[p].pos[i];
                vel[i] +=  s->part[p].mass*s->part[p].vel[i];
            }
            mass +=  s->part[p].mass;
        }
        for(int i = 0; i < 3; i++)
        {
            cmpos[i] = pos[i] / mass;
            cmvel[i] = vel[i] / mass;
        }

    }

    void print_system_energy_diagnostics (struct falcon::particles::particle_system *s) {
            fflush(stdout);
            real_t kinetic = falcon::diagnostics::system_kinetic_energy(s);
            real_t pot = falcon::diagnostics::system_potential_energy(s);
            real_t e0 = -pot+kinetic;
            std::cout<<std::setprecision(15)<<std::scientific<<e0<<std::endl;
    }
    void print_system_com_diagnostics (struct falcon::particles::particle_system *s) {
        real_t cmpos[3], cmvel[3];
        falcon::diagnostics::system_center_of_mass(s, cmpos, cmvel);
        std::cout<<std::fixed<<std::setprecision(15)<<cmpos[0]<<" "<<cmpos[1]<<" "<<cmpos[2]<<"\t";
        fflush(stdout);

    }



}

