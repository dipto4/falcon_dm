/* --------------------------------------------------------------------------------------------------------
 * Falcon-DM: N-body code to simulate Intermediate Mass Ratio Inspirals in DM spikes
 * DOI: 10.1093/mnras/stae1989 
 * 
 * Created by: Diptajyoti Mukherjee
 * Date: Aug 27, 2024
 * 
 * If you decide to use this code please cite: 10.1093/mnras/stae1989
 * 
 * For any comments or questions, please contact me at dipto@cmu.edu
 * or open an issue on github at: https://github.com/dipto4/falcon_dm
 *
 * --------------------------------------------------------------------------------------------------------*/

/* -----------------------------------------------------------------------------------------------------
 * Filename: main.cxx
 * 
 * Purpose: This is the driver code.  
 * The driver code calls functions from falcon.h (this sets global parameters)
 *                                      particle.h (this handles the particle system)
 *                                      io.h (this handles the file input/output)
 *                                      diagnostics.h (this prints useful diagnostics of the overall N-body system)
 *                                      integrator_oop.h (this contains the main file/integrator)
 * 
 * -------------------------------------------------------------------------------------------------------
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <cstring>
#include <cstdlib>

#include <src/falcon.hpp>
#include <src/io.hpp>
#include <src/integrator_oop.hpp>
#include <src/diagnostics.hpp>
#include <src/particle.hpp>
using namespace std;

int main()
{
    
    /*global nbodysystem handler, handles global variables like file name, units, softening, etc*/
    auto [input_fname, restart_val, eta, t_end, dt, use_pn, use_precession, use_radiation, c, soft] = falcon::io::read_config(); 
    falcon::Nbodysystem system {eta, soft,0,t_end, dt, input_fname, use_pn, use_precession, use_radiation, c}; // initialize the global system with the parameters
                                                                                                               // read in from config.ini
    
    struct falcon::particles::particle_system mainsys; 
    
    auto [snapnum, t_now] = falcon::io::load_data(&mainsys, restart_val, input_fname);                        // set up the particle system from the input file
                                                                                                              // in config.ini
    std::cout<<"CoM_x \t\t CoM_y \t\t CoM_z \t\t Energy"<<std::endl;                                          // as diagnostics, we print the center of mom. and energy
    std::cout<<"----------------------------------------------------------"<<std::endl;                       // note the CoM will not be exactly conserved because
    falcon::diagnostics::print_system_com_diagnostics(&mainsys);                                              // of the approximated force calculation of DM particles
    falcon::diagnostics::print_system_energy_diagnostics(&mainsys);                                           // Nevertheless, tests show that the conservation is good enough
    
    
    /* Setting up a new integrator of type HOLD_DKD
     * If different integrators are desired, this line must be changed */
    falcon::integrator::GenericHHSIntegrator *integrator = new falcon::integrator::HOLD_DKD(&system,&mainsys);  
    
    /* run until t_end */
    if(t_end > 0)
    {
        while(t_end > t_now)
        {
            
            fflush(stdout);
            falcon::diagnostics::print_system_com_diagnostics(&mainsys);

            falcon::diagnostics::print_system_energy_diagnostics(&mainsys);
            
            falcon::io::write_hdf5_snapshot(snapnum, t_now,mainsys);
            snapnum++;

            double start_time = omp_get_wtime();
         
            if(mainsys.n > 0)
                integrator->integrate();

            double total_time = omp_get_wtime()-start_time;
            printf("time=%21.16f\n",total_time);

            t_now += dt;
        }
    }


    return 0;
}

