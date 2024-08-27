#ifndef falcon_h
#define falcon_h
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<iostream>


#define HUGE 1.e30
using namespace std;

namespace falcon {
    typedef double real_t;
     
    class Nbodysystem {
        real_t eta_sym {};
        real_t soft {};
        real_t t {0};
        real_t t_end {};
        real_t simtime {};
        real_t dt {};
        real_t C {};
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
                    const std::string filename_init, 
                    bool pn_init, bool precession_init, bool rad_init, const real_t C_init) : 
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
                        include_radiation {rad_init},
                        C {C_init}
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
                            std::cout<<"speed of light: "<<C<<std::endl;
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
            real_t get_C() {return C;}
            void set_simtime(real_t update) {simtime = update;}
            //void get_n();
            

    };
    
}

#endif
