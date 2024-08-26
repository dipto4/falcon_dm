#pragma once
#include<cmath>
#include "falcon.h"
#include "particle.h"
namespace falcon::forces {

    class GenericForce {
        protected:
            const int type {} ;
            GenericForce(int init_type): type {init_type}
            {}
        public:
            GenericForce() = default;
            virtual ~GenericForce() = default;
            virtual void calculateForce(struct falcon::particles::particle &pi,
                    struct falcon::particles::particle &pj,  falcon::Nbodysystem &system) = 0;
    };

    class NewtonianForce : public GenericForce {

        public:
            NewtonianForce() : GenericForce(0)
        {}
            ~NewtonianForce() = default;
            virtual void calculateForce(struct falcon::particles::particle &pi,
                    struct falcon::particles::particle &pj, falcon::Nbodysystem &system);
    };

    class PN1Force : public GenericForce {
        public:
            PN1Force() : GenericForce(2)
        {}
            virtual void calculateForce(struct falcon::particles::particle &pi,
                    struct falcon::particles::particle &pj, falcon::Nbodysystem &system);
    };

    class PN2Force : public GenericForce {
        public: 
            PN2Force() : GenericForce(4)
        {}
            virtual void calculateForce(struct falcon::particles::particle &pi,
                    struct falcon::particles::particle &pj, falcon::Nbodysystem &system);

    };

    class PN25Force : public GenericForce {
        public:
            PN25Force() : GenericForce(5)
        {}
            virtual void calculateForce(struct falcon::particles::particle &pi,
                    struct falcon::particles::particle &pj, falcon::Nbodysystem &system);
    };

    inline void NewtonianForce::calculateForce(struct falcon::particles::particle &p1,struct falcon::particles::particle &p2,  
            falcon::Nbodysystem &system) {
        real_t dr3, dr2, dr;

        real_t dx[3];

        for(size_t d=0;d<3;d++) {
            dx[d] = (p1.pos[d]-p2.pos[d]);

        }
        dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        //only soften the ecnounter between IMBH + DM particle
        if((p1.id > 1) || (p2.id > 1)) {
            if((p1.id == 0) || (p2.id == 0) )
                dr2 += system.get_soft()*system.get_soft();
        }
        //dr2 += EPS*EPS;

        dr = sqrt(dr2);
        dr3 = dr * dr2;
        dr = p2.mass / dr3;

        p1.acc[0] += -dx[0] * dr;
        p1.acc[1] += -dx[1] * dr;
        p1.acc[2] += -dx[2] * dr;


    } 
    inline void PN1Force::calculateForce(struct falcon::particles::particle &p1,struct falcon::particles::particle &p2,  falcon::Nbodysystem &system) {
        real_t C_FACTOR = 1./pow(system.get_C(),type);

        real_t v1[3], v2[3];
        bool use_w = p1.use_w || p2.use_w;
        if(use_w) {
            for(int d=0;d<3;d++) {
                v1[d] = p1.w[d];
                v2[d] = p2.w[d];
            }
        } else {
            for(int d=0; d<3;d++) {
                v1[d] = p1.vel[d];
                v2[d] = p2.vel[d];
            }
        }


        real_t dr[3] = {p2.pos[0]-p1.pos[0],p2.pos[1]-p1.pos[1],p2.pos[2]-p1.pos[2]};
        real_t dv[3] = {v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]};

        real_t dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        real_t r2 = dr2;
        real_t r = sqrt(dr2);
        real_t n[3] = {-dr[0]/r,-dr[1]/r,-dr[2]/r};

        real_t v1s = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
        //real_t v2s = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
        real_t v2s = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
        real_t v12 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

        real_t nv1 = n[0] * v1[0] + n[1] * v1[1] + n[2] * v1[2];
        //real_t nv1 = n[0] * dv[0] + n[1] * dv[1] + n[2] * dv[2];
        real_t nv2 = n[0] * v2[0] + n[1] * v2[1] + n[2] * v2[2];
        real_t gmr1 = p1.mass / r;
        real_t gmr2 = p2.mass / r;

        real_t Ai = -v1s - 2 * v2s + 4 * v12 + 1.5 * nv2 * nv2 + 5 * gmr1 + 4 * gmr2;
        real_t Aj = -v2s - 2 * v1s + 4 * v12 + 1.5 * nv1 * nv1 + 5 * gmr2 + 4 * gmr1;

        real_t Bi = 4 * nv1 - 3 * nv2;
        real_t Bj = -4 * nv2 + 3 * nv1;

        real_t coef = (1.0 / r2) * C_FACTOR;
        real_t a_com[3] = {0.0,0.0,0.0};
        for(int d=0;d<3;d++) {
            p1.acc_pn[d] += (coef * p2.mass) * (Ai * n[d] - Bi * dv[d]);
            p1.acc_pn[d] -= (coef * p1.mass) * (Aj * n[d] - Bj * dv[d]);

            a_com[d] = (p1.mass*p1.acc_pn[d] + p2.mass*p2.acc_pn[d])/(p1.mass + p2.mass);

        }

        for(size_t d=0;d<3;d++) {
            p1.acc_pn[d] -= a_com[d];//coef * (Ai * n[d] - Bi * dv[d]);
            p2.acc_pn[d] -= a_com[d];//coef * (Aj * n[d] - Bj * dv[d]);

            //a_com[d] = (parti.mass*acc_pn_i[d] + partj.mass*acc_pn_j[d])/(parti.mass + partj.mass);
        }



    }
    inline void PN2Force::calculateForce(struct falcon::particles::particle &p1,struct falcon::particles::particle &p2,  falcon::Nbodysystem &system) {
        real_t C_FACTOR = 1./pow(system.get_C(),type);

        real_t v1[3], v2[3];
        bool use_w = p1.use_w || p2.use_w;
        if(use_w) {
            for(int d=0;d<3;d++) {
                v1[d] = p1.w[d];
                v2[d] = p2.w[d];
            }
        } else {
            for(int d=0; d<3;d++) {
                v1[d] = p1.vel[d];
                v2[d] = p2.vel[d];
            }
        }
        real_t dr[3] = {p2.pos[0]-p1.pos[0],p2.pos[1]-p1.pos[1],p2.pos[2]-p1.pos[2]};
        real_t dv[3] = {v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]};

        //real_t dv[3] = {partj.vel[0]-parti.vel[0],partj.vel[1]-parti.vel[1],partj.vel[2]-parti.vel[2]};

        real_t dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        real_t r2 = dr2;
        real_t r = sqrt(dr2);
        real_t n[3] = {-dr[0]/r,-dr[1]/r,-dr[2]/r};

        real_t v1s = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
        real_t v1q = v1s * v1s;
        //real_t v2s = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
        real_t v2s = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
        real_t v2q = v2s * v2s;
        real_t v12 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

        real_t nv1 = n[0] * v1[0] + n[1] * v1[1] + n[2] * v1[2];
        real_t nv1s = nv1 * nv1;
        //real_t nv1 = n[0] * dv[0] + n[1] * dv[1] + n[2] * dv[2];
        real_t nv2 = n[0] * v2[0] + n[1] * v2[1] + n[2] * v2[2];
        real_t nv2s = nv2 * nv2;
        real_t gmr1 = p1.mass / r;
        real_t gmr2 = p2.mass / r;

        real_t m1s = p1.mass * p1.mass;
        real_t m2s = p2.mass * p2.mass;
        real_t m12 = p1.mass * p2.mass;

        real_t Ai = -2 * v2q + 4 * v2s * v12 - 2 * v12 * v12 +
            nv2s * (1.5 * v1s + 4.5 * v2s - 6 * v12 - 1.875 * nv2s) +
            gmr1 * (-3.75 * v1s + 1.25 * v2s - 2.5 * v12 + 19.5 * nv1s - 39 * nv1 * nv2 + 8.5 * nv2s) +
            gmr2 * (4 * v2s - 8 * v12 + 2 * nv1s - 4 * nv1 * nv2 - 6 * nv2s) +
            1.0 / dr2 * (-14.25 * m1s - 9 * m2s - 34.5 * m12);

        real_t Aj = -2 * v1q + 4 * v1s * v12 - 2 * v12 * v12 +
            nv1s * (1.5 * v2s + 4.5 * v1s - 6 * v12 - 1.875 * nv1s) +
            gmr2 * (-3.75 * v2s + 1.25 * v1s - 2.5 * v12 + 19.5 * nv2s - 39 * nv1 * nv2 + 8.5 * nv1s) +
            gmr1 * (4 * v1s - 8 * v12 + 2 * nv2s - 4 * nv1 * nv2 - 6 * nv1s) +
            1.0 / r2 * (-14.25 * m2s - 9 * m1s - 34.5 * m12);

        real_t Bi = v1s * nv2 + 4 * v2s * nv1 - 5 * v2s * nv2 - 4 * v12 * nv1 + 4 * v12 * nv2 - 6 * nv1 * nv2s +
            4.5 * nv2 * nv2s + gmr1 * (-15.75 * nv1 + 13.75 * nv2) + gmr2 * (-2 * nv1 - 2 * nv2);

        real_t Bj = -v2s * nv1 - 4 * v1s * nv2 + 5 * v1s * nv2 + 4 * v12 * nv2 - 4 * v12 * nv1 + 6 * nv2 * nv1s -
            4.5 * nv1 * nv1s + gmr2 * (15.75 * nv2 - 13.75 * nv1) + gmr1 * (2 * nv2 + 2 * nv1);

        real_t coef = 1.0 / dr2 * C_FACTOR;
        real_t a_com[3] = {0.0,0.0,0.0};
        for(int d=0;d<3;d++) {
            p1.acc_pn[d] += (coef * p2.mass) * (Ai * n[d] - Bi * dv[d]);
            p2.acc_pn[d] -= (coef * p1.mass) * (Aj * n[d] - Bj * dv[d]);

            a_com[d] = (p1.mass*p1.acc_pn[d] + p2.mass*p2.acc_pn[d])/(p1.mass + p2.mass);

        }

        for(size_t d=0;d<3;d++) {
            p1.acc_pn[d] -= a_com[d];//coef * (Ai * n[d] - Bi * dv[d]);
            p2.acc_pn[d] -= a_com[d];//coef * (Aj * n[d] - Bj * dv[d]);

        }


    }
    inline void PN25Force::calculateForce(struct falcon::particles::particle &p1,struct falcon::particles::particle &p2,  falcon::Nbodysystem &system) {
        real_t C_FACTOR = 1./pow(system.get_C(),type);

        real_t v1[3], v2[3];
        /*needed for the auxiliary velocity algorithm*/
        bool use_w = p1.use_w || p2.use_w;
        if(use_w) {
            for(int d=0;d<3;d++) {
                v1[d] = p1.w[d];
                v2[d] = p2.w[d];
            }
        } else {
            for(int d=0; d<3;d++) {
                v1[d] = p1.vel[d];
                v2[d] = p2.vel[d];
            }
        }

        real_t dr[3] = {p2.pos[0]-p1.pos[0],p2.pos[1]-p1.pos[1],p2.pos[2]-p1.pos[2]};
        real_t dv[3] = {v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]};

        real_t r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        real_t r = sqrt(r2);
        real_t dv2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
        real_t n[3] = {-dr[0]/r, -dr[1]/r, -dr[2]/r};

        real_t n_dot_v = -(n[0]*dv[0] + n[1]*dv[1]+n[2]*dv[2]);

        real_t gmr1 = p1.mass / r;
        real_t gmr2 = p2.mass / r;

        real_t Ai = n_dot_v * (3 * dv2 - 6 * gmr1 + 52.0 / 3 * gmr2);
        real_t Aj = n_dot_v * (3 * dv2 - 6 * gmr2 + 52.0 / 3 * gmr1);

        real_t Bi = -dv2 + 2 * gmr1 - 8 * gmr2;
        real_t Bj = -dv2 + 2 * gmr2 - 8 * gmr1;

        real_t coef = (4.0/5.0)  * p1.mass * p2.mass / (r2 * r) * C_FACTOR;
        real_t a_com[3] = {0.0,0.0,0.0};
        for(size_t d=0;d<3;d++) {
            p1.acc_pn[d] += coef * (Ai * n[d] - Bi * dv[d]);
            p2.acc_pn[d] -= coef * (Aj * n[d] - Bj * dv[d]);

            a_com[d] = (p1.mass*p1.acc_pn[d] + p2.mass*p2.acc_pn[d])/(p1.mass + p2.mass);
        }
        // This is a quick fix for centering the object at the CoM every step. Performed tests to ensure that results do not
        // change due to this fix
        for(size_t d=0;d<3;d++) {
            p1.acc_pn[d] -= a_com[d];
            p2.acc_pn[d] -= a_com[d];//

        }


    }
}
