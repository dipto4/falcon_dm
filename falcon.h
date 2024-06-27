#ifndef falcon_h
#define falcon_h
//#include <complex>
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<iostream>


#define HUGE 1.e30
#define ETA_SYM 0.025
#define ETA_EXT 0.0125
using namespace std;

namespace falcon {
    typedef double real_t;
     
    /*struct particle {
        real_t pos[3];
        real_t mass;
        real_t vel[3];
        real_t w[3]; //needed for Auxiliary Velocity Algorithm with PN terms. See Hellstrom and Mikkola (2010) for more details 
        real_t t_last;
        real_t timestep;
        real_t acc[3];
        real_t acc_pn[3];
        real_t jerk[3];
        size_t id;
        real_t pot;
    };

    struct sys {
        size_t n;
        struct particle *part;
        struct particle *last;
    };

    struct sys mainsys;
    struct sys copysys;

    static char input_fname[200];
    int snapnum;
    real_t t_now;
    size_t numBodies;


    struct external {
        std::string pot_type;
        std::string geometry_type;
        bool include_df;
        bool include_pn;
        bool include_pot;
        real_t mass;
        
        real_t pos[3];
        
        //needed for GW
        real_t vel[3];
        real_t acc[3];

        

        //for NFW
        real_t rs;
        real_t c;

        for Dehnen
        real_t gamma;
        real_t lambda;
    };*/




    class Nbodysystem {
        real_t eta_sym {};
        real_t soft {};
        real_t t {0};
        real_t t_end {};
        real_t simtime {};
        real_t dt {};
        //real_t eta_ext {}; 
        size_t n {};
        std::string filename;
        bool include_pn {};
        bool include_precession {}; //use PN1 & PN2 terms
        bool include_radiation {}; //use PN2.5 terms
        public:
            ~Nbodysystem () = default ;
            Nbodysystem () : 
                eta_sym {0.025},
                t_end {-1},
                dt {0.01},        //eta_ext {0.0125},
                n {0},
                soft {0},
                filename {"none"},
                include_pn {true},
                include_precession {true},
                include_radiation {true}
            {}

            Nbodysystem (const real_t eta_sym_init, const real_t soft_init,//const real_t eta_ext_init,
                    const size_t n_init, 
                    const real_t t_end_init, const real_t dt_init,
                    const std::string filename_init, bool pn_init, bool precession_init, bool rad_init) : 
                        eta_sym {eta_sym_init},
                        soft {soft_init},
                        //eta_ext {eta_ext_init},
                        t_end {t_end_init},
                        simtime {0.0},
                        dt {dt_init},
                        n {n_init}, 
                        filename {filename_init},
                        include_pn {pn_init},
                        include_precession {precession_init},
                        include_radiation {rad_init}
                        {
                            /* print some diagnostics */
                            std::cout<<"Falcon starting with the following values:"<<std::endl;
                            std::cout<<"eta: "<<eta_sym<<std::endl;
                            std::cout<<"soft: "<<soft<<std::endl;
                            std::cout<<"t_end: "<<t_end<<std::endl;
                            std::cout<<"dt: "<<dt<<std::endl;
                            std::cout<<"using PN: "<<include_pn<<std::endl;
                            std::cout<<"precession (PN1 + PN2): "<<include_precession<<std::endl;
                            std::cout<<"radiation (PN2.5): "<<include_radiation<<std::endl;
                            std::cout<<"filename: "<<filename<<std::endl;
                        }
            void set_n(const size_t n_set) {n = n_set;}
            
            size_t get_n() {return n;}
            real_t get_eta_sym() {return eta_sym;}
            real_t get_soft() {return soft;}
            real_t get_dt() {return dt;} 
            bool use_pn() {return include_pn;}
            bool use_precession() {return include_precession;}
            bool use_radiation() {return include_radiation;}
            real_t get_simtime() {return simtime;}
            void set_simtime(real_t update) {simtime = update;}
            //void get_n();
            

    };
    
}

#endif
