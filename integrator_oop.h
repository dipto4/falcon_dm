/*------------------------------------------------------------------------------------------------
 * Filename: integrator_oop.h
 *
 * Purpose: Defines the integrator and associated helper functions 
 *
 * Main classes:
 * class GenericHHSIntegrator (Abstract class for a Hamiltonian splitting integrator)
 * 
 * class HOLD_DKD (Derived class of GenericHHSIntegrator for the 2nd order hamiltonian splitting based HOLD 
 *                  scheme. Uses the auxiliary velocity algorithm (Hellstrom and Mikkola 2010) for the in-
 *                  tegration of non-separable post-Newtonian forces)
 *
 * Main functions:
 * void split(real_t dt, struct falcon::particles::particle_system s,struct falcon::particles::particle_system *slow, struct falcon::particles::particle_system *fast)
 * void kick_sf(struct falcon::particles::particle_system sinks, struct falcon::particles::particle_system sources, real_t dt, bool *includes_bh, size_t *parti_pn, size_t *partj_pn)
 * void kick_self(struct falcon::particles::particle_system sinks, bool *includes_bh, size_t *parti_pn, size_t *partj_pn)
 * void drift(struct falcon::particles::particle_system s, real_t dt)
 * void kick_sync(struct falcon::particles::particle_system s, real_t dt)
 * void kick_sync_w(struct falcon::particles::particle_system s, real_t dt)
 * void findtimesteps(struct falcon::particles::particle_system s)
 * virtual void step(size_t clevel, falcon::particles::particle_system total, real_t stime, real_t etime, real_t dt, bool calc_timestep)
 -------------------------------------------------------------------------------------------------*/
#pragma once
#include<cmath>
#include <limits>
#include "falcon.h"
#include "forces.h"
#include "particle.h"

/*some legacy macros for integrators*/
#define ENDRUN(fmt, ...) {                \
    printf("ENDRUN at %s:%d ", __FILE__, __LINE__);    \
    printf(fmt, ## __VA_ARGS__);            \
    fflush(stdout);                    \
    exit(-1);                        \
}


#define SWAP(a,b,c) {c t;t=(a);(a)=(b);(b)=t;}
#define ABS(X) (((X) >= 0) ? (X) :-(X))
#define LOG(fmt, ...) {				\
    printf("%s:%d\t", __FILE__, __LINE__);	\
    printf(fmt, ## __VA_ARGS__);		\
}


namespace falcon::integrator {

    class GenericHHSIntegrator {
        protected:
            /*it should take an Nbodysystem containing the global params,
             * and a particle system containing the particle system
             * ideally one would create a template that can handle two different kinds of particle systems
             * one that consists of w (auxiliary velocity) terms and one that does not*/
            falcon::Nbodysystem *system;
            falcon::particles::particle_system *partsys;
            falcon::forces::GenericForce *f;
            bool handle_pn {};
            /*some functions common to HHS integrators follow*/
            falcon::particles::particle_system join(struct falcon::particles::particle_system sinks, 
                    struct falcon::particles::particle_system sources); 


            void split(real_t dt, struct falcon::particles::particle_system s,
                    struct falcon::particles::particle_system *slow, struct falcon::particles::particle_system *fast); 

            /*kick and drift functions next*/

            void kick_sf(struct falcon::particles::particle_system sinks, struct falcon::particles::particle_system sources,
                    real_t dt, bool *includes_bh, size_t *parti_pn, size_t *partj_pn);

            /*function can be eliminated or refractored into the previous one*/
            void kick_self(struct falcon::particles::particle_system sinks, bool *includes_bh, size_t *parti_pn, size_t *partj_pn);

            void drift(struct falcon::particles::particle_system s, real_t dt);

            /*setup in this form to handle some legacy code*/
            void kick_sync(struct falcon::particles::particle_system s, real_t dt);
            void kick_sync_w(struct falcon::particles::particle_system s, real_t dt);
            void kick_pn(struct falcon::particles::particle_system s1, 
                    struct falcon::particles::particle_system s2, size_t pi, size_t pj, bool use_w);

            void findtimesteps(struct falcon::particles::particle_system s);

            //falcon::particles::ParticleSystem sys
            virtual void step(size_t clevel, falcon::particles::particle_system total, 
                    real_t stime, real_t etime, real_t dt, bool calc_timestep) = 0; 
        public:

            GenericHHSIntegrator() = default;
            ~GenericHHSIntegrator() = default;

            GenericHHSIntegrator(falcon::Nbodysystem *init_system, falcon::particles::particle_system *init_partsys) {
                system = init_system;
                partsys = init_partsys;
            }
            void integrate() {
                if(system->get_dt()==0)
                    return;

                size_t clevel = 0;
                real_t dt = system->get_dt();
                real_t simtime = system->get_simtime();
                findtimesteps(*partsys);
                step(clevel,*partsys, simtime, simtime+dt, dt,true);
            }



    };



    class HOLD_DKD : public GenericHHSIntegrator {
        protected:
            virtual void step(size_t clevel, falcon::particles::particle_system total, 
                    real_t stime, real_t etime, real_t dt, bool calc_timestep);


        public:
            HOLD_DKD(falcon::Nbodysystem *init_system, falcon::particles::particle_system *init_partsys) : GenericHHSIntegrator(init_system, init_partsys) {
                handle_pn = true;
            }


    };

    inline falcon::particles::particle_system GenericHHSIntegrator::join(struct falcon::particles::particle_system sinks, 
                                                                         struct falcon::particles::particle_system sources) {
        struct falcon::particles::particle_system s = {0, nullptr, nullptr};
        if(sinks.n == 0)
            return sources;
        if(sources.n == 0)
            return sinks;
        s.n = sinks.n + sources.n;
        if(sinks.part + sinks.n == sources.part)
        {
            s.part = sinks.part;
            s.last = sources.last;
        }
        else
        {
            if(sources.part + sources.n == sinks.part)
            {
                s.part = sources.part;
                s.last = sinks.last;
            }
            else
                ENDRUN("join error 1");
        }
        if(s.last-s.part + 1 != s.n)
            ENDRUN("join error 2");
        return s;

    }

    inline void GenericHHSIntegrator::split(real_t dt, 
                                            struct falcon::particles::particle_system s,
                                            falcon::particles::particle_system *slow, struct falcon::particles::particle_system *fast) {
        size_t i = 0;
        struct falcon::particles::particle *left, *right;
        left = s.part;
        right = s.last;
        dt = fabs(dt);
        while(1)
        {
            if(i >= s.n){
                printf("i:%d s.n:%d dt=%g left=%g right=%g \n", i, s.n, dt, left->timestep, right->timestep);
                ENDRUN("split error 1\n");
            }
            i++;
            while(left->timestep < dt && left < right)
                left++;
            while(right->timestep >= dt && left < right)
                right--;
            if(left < right)
            {
                SWAP(*left, *right, struct falcon::particles::particle);
            }
            else
                break;
        }
        if(left->timestep < dt)
            left++;
        slow->n = s.last-left + 1;
        fast->n = (left-s.part);

        if(fast->n == 1){
            fast->n = 0;
            slow->n = s.n;
        }

        if(slow->n > 0){
            slow->part = s.part + fast->n;
            slow->last = s.last;
        }

        if(fast->n > 0)
        {
            fast->part = s.part;
            fast->last = s.part + fast->n-1;}

        if(fast->n + slow->n != s.n)
            ENDRUN("split error 2\n");

    }

    inline void GenericHHSIntegrator::kick_self(struct falcon::particles::particle_system sinks, 
                                                bool *includes_bh, size_t *parti_pn, size_t *partj_pn) {

        f = new falcon::forces::NewtonianForce();

        for(size_t i=0; i<sinks.n; i++) {
            sinks.part[i].acc[0] = 0.0; sinks.part[i].acc[1] = 0.0; sinks.part[i].acc[2] = 0.0;
            for(size_t j=0; j<sinks.n; j++) {
                if((sinks.part[i].id == sinks.part[j].id) || (sinks.part[j].id > 1 && sinks.part[i].id > 1))
                    continue;

                f->calculateForce(sinks.part[i], sinks.part[j], *system);

                if((sinks.part[i].id==0 && sinks.part[j].id==1) || (sinks.part[i].id==1 && sinks.part[j].id==0)) {
                    *includes_bh = true;
                    *parti_pn = i;
                    *partj_pn = j;
                }
            }


        }

        delete f;


    }


    inline void GenericHHSIntegrator::kick_sf(struct falcon::particles::particle_system sinks, struct falcon::particles::particle_system sources,
                                              real_t dt, 
                                              bool *includes_bh, size_t *parti_pn, size_t *partj_pn) {
        // sources -> sinks first
        // sinks -> sources second

        f = new falcon::forces::NewtonianForce();

        for(size_t i=0; i<sinks.n ; i++) {
            sinks.part[i].acc[0] = 0.0; sinks.part[i].acc[1] = 0.0; sinks.part[i].acc[2] = 0.0;

            for(size_t j=0; j<sources.n; j++) {
                if(sinks.part[i].id > 1 && sources.part[j].id > 1)
                    continue;

                f->calculateForce(sinks.part[i],sources.part[j], *system);
                if((sinks.part[i].id==0 && sources.part[j].id==1) || (sinks.part[i].id==1 && sources.part[j].id==0)) {
                    *includes_bh = true;
                    *parti_pn = i;
                    *partj_pn = j;
                }




            }


        } 

        delete f;

    }

    inline void GenericHHSIntegrator::drift(struct falcon::particles::particle_system s, real_t dt) {
        for(size_t i=0;i<s.n;i++) {
            for(size_t d=0;d<3;d++) {
                s.part[i].pos[d] += dt * s.part[i].vel[d];
            }
        }

    }

    inline void GenericHHSIntegrator::kick_sync(struct falcon::particles::particle_system s, real_t dt) {
        for(size_t i=0;i<s.n;i++) {
            for(size_t d=0;d<3;d++) {     // add the newtonian and post newtonian part separately
                s.part[i].vel[d] += dt * (s.part[i].acc[d]+s.part[i].acc_pn[d]);
            }
        }
    }

    inline void GenericHHSIntegrator::kick_sync_w(struct falcon::particles::particle_system s, real_t dt) {
        for(size_t i=0;i<s.n;i++) {
            for(size_t d=0;d<3;d++) {   // add the newtonian and post newtonian part separately
                s.part[i].w[d] += dt * (s.part[i].acc[d] + s.part[i].acc_pn[d]);
            }
        }

    }

    inline void GenericHHSIntegrator::kick_pn(struct falcon::particles::particle_system s1, 
            struct falcon::particles::particle_system s2, size_t pi, size_t pj, bool use_w) {

        s1.part[pi].acc_pn[0] = 0.0; s1.part[pi].acc_pn[1] = 0.0; s1.part[pi].acc_pn[2] = 0.0;   
        s2.part[pj].acc_pn[0] = 0.0; s2.part[pj].acc_pn[1] = 0.0; s2.part[pj].acc_pn[2] = 0.0;   

        s1.part[pi].use_w = use_w;
        s2.part[pj].use_w = use_w;

        if (system->use_precession()) {
            f = new falcon::forces::PN1Force();
            f->calculateForce(s1.part[pi], s2.part[pj], *system);
            delete f;

            f = new falcon::forces::PN2Force();
            f->calculateForce(s1.part[pi], s2.part[pj], *system);
            delete f;
        }

        if (system->use_radiation()) {
            f = new falcon::forces::PN25Force();
            f->calculateForce(s1.part[pi], s2.part[pj], *system);
            delete f;
        }
    }

    inline void GenericHHSIntegrator::findtimesteps(struct falcon::particles::particle_system s) {
        /*one would want the imbh and the compact object to have the saame timestep 
         * for correct pn evolution. It reduces the code efficiency a bit but produces the 
         * correct results*/
        size_t part0=0;
        size_t part1=0;

        for(size_t i=0;i<s.n;i++) {
            real_t stepi = std::numeric_limits<real_t>::max();
            //#pragma omp parallel for if(s.part[i].id==0 && s.n > 100) reduction(min:stepi)
            for(size_t j=0;j<s.n;j++) {
                if((s.part[i].id == s.part[j].id) || (s.part[j].id > 1 && s.part[i].id > 1)) 
                    continue;

                real_t dr[3], dv[3];
                //real_t dr2, vdotdr2;

                for(size_t d=0;d<3;d++) {
                    dr[d] = s.part[i].pos[d] - s.part[j].pos[d];
                    dv[d] = s.part[i].vel[d] - s.part[j].vel[d];

                }

                real_t dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

                /*only soften interactions between the IMBH and DM particles */
                if((s.part[i].id > 1 || s.part[j].id > 1)) {
                    if((s.part[i].id ==0) || (s.part[j].id == 0))
                        dr2 += system->get_soft() * system->get_soft();//soft * soft;
                }


                real_t v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];

                real_t r = sqrt(dr2);
                real_t vdotdr2 = (dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2])/dr2;

                real_t tau = system->get_eta_sym() / M_SQRT2 * r * sqrt(r/(s.part[i].mass+ s.part[j].mass));
                real_t dtau = 3 * tau * vdotdr2 / 2;

                dtau = (dtau < 1) ? dtau : 1;
                tau /= (1 - dtau / 2);

                stepi = std::min(tau, stepi);

                tau = system->get_eta_sym() * r / sqrt(v2);
                dtau = tau * vdotdr2 * (1 + (s.part[i].mass+s.part[j].mass) / (v2 * r));
                dtau = (dtau < 1)? dtau: 1;

                tau /= (1 - dtau / 2);
                stepi = std::min(tau, stepi);

            }
            if(s.part[i].id==0) part0=i;
            if(s.part[i].id==1) part1=i;

            s.part[i].timestep = stepi;
        }

        real_t min_step = std::min(s.part[part0].timestep, s.part[part1].timestep);
        s.part[part0].timestep = min_step;
        s.part[part1].timestep = min_step;
        s.part[part1].timestep=s.part[part0].timestep;

    }



    inline void HOLD_DKD::step(size_t clevel, falcon::particles::particle_system total,
            real_t stime, real_t etime, real_t dt, bool calc_timestep) {

        falcon::particles::particle_system slow = {0,nullptr,nullptr}, fast = {0,nullptr,nullptr};
        bool includes_bh_self = false, includes_bh_sf = false; //BH-BH interaction either present in S or F
        size_t parti_pn, partj_pn;
        if(calc_timestep) {
            findtimesteps(total);
        }

        split(dt, total, &slow, &fast);

        if(fast.n == 0)
        {
            system->set_simtime( system->get_simtime()+ dt);
#if DEBUG
            printf("level=%i, t=%g s=%d \n",clevel, diag->simtime, total.n);
            fflush(stdout);
#endif
        }

        //hold for fast system
        if(fast.n > 0)
            step(clevel+1, fast, stime, stime+dt/2, dt/2, false);


        if(slow.n > 0)
            drift(slow, dt / 2);

        if(slow.n > 0) { //kicksf in between
            kick_self(slow,  &includes_bh_self, &parti_pn, &partj_pn);
            //add pn terms here if both are in slow
            if(fast.n > 0) {
                kick_sf( slow, fast, dt, &includes_bh_sf, &parti_pn, &partj_pn);
                kick_sf(fast, slow, dt, &includes_bh_sf, &parti_pn, &partj_pn);
            }

            // Check Hellstrom and Mikkola (2010) for the integration scheme
            if((includes_bh_self || includes_bh_sf )&& system->use_pn()) {
                if(includes_bh_self)
                    kick_pn(slow, slow, parti_pn, partj_pn,false);
                else if(includes_bh_sf)
                    kick_pn(fast,slow,parti_pn,partj_pn,false);

            }
            kick_sync_w(total, dt/2.);

            if((includes_bh_self || includes_bh_sf )&& system->use_pn()) {
                if(includes_bh_self)
                    kick_pn(slow, slow, parti_pn, partj_pn,true);
                else if(includes_bh_sf)
                    kick_pn(fast,slow,parti_pn,partj_pn,true);
            }

            kick_sync(total, dt);

            if((includes_bh_self || includes_bh_sf )&& system->use_pn()) {
                if(includes_bh_self)
                    kick_pn(slow, slow, parti_pn, partj_pn,false);
                else
                    kick_pn(fast,slow,parti_pn,partj_pn,false);

            }
            kick_sync_w(total,dt/2.);

        }

        if(slow.n > 0)
            drift( slow,  dt / 2);
        //hold for fast system
        if(fast.n > 0)
            step(clevel+1, fast, stime+dt/2, etime, dt/2, true);


    }


}
