#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <cstring>
#include <cstdlib>

/*to read ini files*/

#include "falcon.h"
//#include "integrator.h"
#include "io.h"
//#include "external.h"
#include "integrator_oop.h"
#include "diagnostics.h"
#include "particle.h"
using namespace std;
//using namespace falcon;




int main()
{
    /*
    
    */
    /*global nbodysystem handler, handles global variables like file name, units, softening, etc*/
    auto [input_fname, restart_val, eta, t_end, dt, use_pn, use_precession, use_radiation, c, soft] = falcon::io::read_config(); 
    falcon::Nbodysystem system {eta, soft,0,t_end, dt, input_fname, use_pn, use_precession, use_radiation};
    
    struct falcon::particles::particle_system mainsys; 
    
    auto [snapnum, t_now] = falcon::io::load_data(&mainsys, restart_val, input_fname);
    

    falcon::diagnostics::print_system_energy_diagnostics(&mainsys);
    falcon::diagnostics::print_system_com_diagnostics(&mainsys);

    falcon::integrator::GenericHHSIntegrator *integrator = new falcon::integrator::HOLD_DKD(&system,&mainsys);  
     
    if(t_end > 0)
    {
        while(t_end > t_now)
        {
            
            //double start_time = omp_get_wtime();
            fflush(stdout);
            falcon::diagnostics::print_system_energy_diagnostics(&mainsys);
            falcon::diagnostics::print_system_com_diagnostics(&mainsys);

            falcon::io::write_hdf5_snapshot(snapnum, t_now,mainsys);
            snapnum++;

            double start_time = omp_get_wtime();
         
            if(mainsys.n > 0)
                integrator->integrate();

            //update_step=true;//false;
            double total_time = omp_get_wtime()-start_time;
            printf("time=%21.16f\n",total_time);

            t_now += dt;
        }
    }


    return 0;
}

