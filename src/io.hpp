/* --------------------------------------------------------------------------------------------------
 * Filename: io.h
 * 
 * Purpose: Handles all the important input/output operations, including reading/writing snapshots
 *          and reading config.ini
 *          NOTE: You need the boost library to be able to use the IO operations
 *          Boost is utilized to read ini files
 *
 * Main functions:
 *
 * std::tuple<std::string, int, real_t, real_t, real_t, bool, bool, bool, real_t, real_t> read_config()
 *
 * read_hdf5_snapshot(std::string input_fname, int *snapnum, real_t *t_now,struct falcon::particles::particle_system *s)
 *
 * write_hdf5_snapshot(int snapshot_num, real_t t_now,struct falcon::particles::particle_system s)
 *
 * std::tuple<int, real_t> load_data(struct falcon::particles::particle_system *s, const size_t restart_val, std::string input_fname)
 *
 --------------------------------------------------------------------------------------------------*/
#pragma once

#include <numeric>		
#include <algorithm>		
#include "falcon.hpp"
#include "particle.hpp"
#include "hdf5.h"
#include <filesystem>
#include <tuple>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


namespace falcon::io
{


    /* ---------------------------------------------------------------------------------
     * Function: std::tuple<std::string, int, real_t, real_t, real_t, bool, bool, bool, real_t, real_t> read_config()
     * 
     * Input: None (reads directly from the config.ini file. The config.ini MUST BE in the simulation directory)
     * Output: Tuple of type <std::string, int, real_t, real_t, real_t, bool, bool, bool, real_t, real_t>
     *         containing values of <filename, restart_val, eta, t_end, dt, use_pn, use_precession, use_radiation, c, soft>
     -----------------------------------------------------------------------------------*/


    std::tuple<std::string, int, real_t, real_t, real_t, bool, bool, bool, real_t, real_t> read_config() {


        //check if config.ini exists



        boost::property_tree::ptree pt;

        /* The code tries to read from a file called config.ini. If it is not present, the code CANNOT be run! 
         * Check the documentation to see how to write a config.ini file
         * */
        try {

            boost::property_tree::ini_parser::read_ini("config.ini", pt);
        }
        catch (boost::property_tree::ptree_bad_path& ex) {
            std::cerr<<ex.what()<<std::endl<<"ERROR: config.ini does not exist! Please refer to the documentation to read how to make a proper config file!"
                <<std::endl;
            //std::cout<<"ERROR: config.ini does not exist! Please refer to the documentation to read how to make a proper config file!"<<std::endl;
            //std::cout<<"Falcon exiting!"<<std::endl;
        }
        std::string input_fname { pt.get<std::string>("Globals.filename") };
        int restart_val {std::stoi(pt.get<std::string>("Globals.restart"))};
        real_t eta {std::stod(pt.get<std::string>("Globals.eta"))};
        real_t t_end {std::stod(pt.get<std::string>("Globals.t_end"))};
        real_t dt {std::stod(pt.get<std::string>("Globals.dt"))};
        real_t c {std::stod(pt.get<std::string>("Globals.c"))};
        real_t soft {std::stod(pt.get<std::string>("Globals.soft"))};


        /*for PN terms */
        bool use_pn  = std::stoi(pt.get<std::string>("Globals.use_pn"));
        bool use_precession = std::stoi(pt.get<std::string>("Globals.use_precession"));
        bool use_radiation =  std::stoi(pt.get<std::string>("Globals.use_radiation"));



        return {input_fname, restart_val, eta,t_end,dt,use_pn,use_precession, use_radiation, c, soft};
    }



    /* -----------------------------------------------------------------------------------
     * Function: read_hdf5_snapshot(std::string input_fname, int *snapnum, real_t *t_now,struct falcon::particles::particle_system *s)
     *
     * Inputs: std::string (input filename)
     *         int* (pointer to the snapnumber that we are trying to read in case of a restart)
     *         real_t* (pointer to the current time from the restart file)
     *         falcon::particles::particle_system* (pointer to the particle system to read in from the restart file)
     * Output: None/void (particle dataset is loaded on to the passed array)
     --------------------------------------------------------------------------------------*/

    /* returns the snapnum and t_now*/
    void read_hdf5_snapshot(std::string input_fname, int *snapnum, real_t *t_now,struct falcon::particles::particle_system *s)
    {
        /*needed for some legacy io stuff*/
        static char dummyString[200];
        size_t numBodies {};
        strncpy(dummyString, input_fname.c_str() + 9, input_fname.length() - 5 - 9);
        *snapnum = atoi(dummyString);

        herr_t status;
        hid_t file_id, hdf5_headergrp, hdf5_attribute, dataset;

        file_id = H5Fopen(input_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hdf5_headergrp = H5Gopen(file_id, "/Header");

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPartTotal");
        H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &numBodies);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
        H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, t_now);
        H5Aclose(hdf5_attribute);

        H5Gclose(hdf5_headergrp);

        s->n = numBodies;
        s->part = new struct falcon::particles::particle[numBodies];//(struct falcon::particles::particle *) malloc(numBodies * sizeof(struct falcon::particles::particle));

        s->last = &s->part[numBodies - 1];

        real_t *data = new real_t[numBodies];

        dataset = H5Dopen(file_id, "/Posx/");
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        status = H5Dclose(dataset);
        for(unsigned int b = 0; b < numBodies; b++)
            s->part[b].pos[0] = data[b];

        dataset = H5Dopen(file_id, "/Posy/");
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        status = H5Dclose(dataset);
        for(unsigned int b = 0; b < numBodies; b++)
            s->part[b].pos[1] = data[b];

        dataset = H5Dopen(file_id, "/Posz/");
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        status = H5Dclose(dataset);
        for(unsigned int b = 0; b < numBodies; b++)
            s->part[b].pos[2] = data[b];

        dataset = H5Dopen(file_id, "/Velx/");
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        status = H5Dclose(dataset);
        for(unsigned int b = 0; b < numBodies; b++)
            s->part[b].vel[0] = data[b];

        dataset = H5Dopen(file_id, "/Vely/");
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        status = H5Dclose(dataset);
        for(unsigned int b = 0; b < numBodies; b++)
            s->part[b].vel[1] = data[b];

        dataset = H5Dopen(file_id, "/Velz/");
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        status = H5Dclose(dataset);
        for(unsigned int b = 0; b < numBodies; b++)
            s->part[b].vel[2] = data[b];

        dataset = H5Dopen(file_id, "/Mass/");
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        status = H5Dclose(dataset);
        for(unsigned int b = 0; b < numBodies; b++)
            s->part[b].mass = data[b];



        delete[] data;

        unsigned int *id = new unsigned int[numBodies];
        dataset = H5Dopen(file_id, "/ID/");
        status = H5Dread(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        status = H5Dclose(dataset);
        for(unsigned int b = 0; b < numBodies; b++)
            s->part[b].id = id[b];
        delete[] id;

        status = H5Fclose(file_id);

        if(status < 0)
            printf("Reading snapshot error\n");

        printf("Reading IC done with %d particles \n", numBodies);
        fflush(stdout);
    }

    /*-----------------------------------------------------------------------------------
     * Function: write_hdf5_snapshot(int snapshot_num, real_t t_now,struct falcon::particles::particle_system s)
     *
     * Inputs: int (snapshot number of the HDF5 file)
     *         real_t (current time)
     *         falcon::particles::particle_system (particle system to write to the HDF5 file)
     * Outputs: None
     --------------------------------------------------------------------------------------*/
    void write_hdf5_snapshot(int snapshot_num, real_t t_now,struct falcon::particles::particle_system s)
    {
        size_t numBodies = s.n;
        hsize_t dims[1] = { s.n };
        herr_t status;
        hid_t file_id, space_id, dset_id, memspace, handle = 0;
        char buf[500];
        sprintf(buf, "snapshot_%d.hdf5", snapshot_num);

        double *data = new double[s.n];

        file_id = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // first write the header information
        handle = H5Gcreate(file_id, "/Header", 0);

        hid_t hdf5_dataspace, hdf5_attribute;

        //    printf("Writing %d particles at t=%g \n", numBodies, t_now);
        //    fflush(stdout);

        hdf5_dataspace = H5Screate(H5S_SCALAR);
        hdf5_attribute = H5Acreate(handle, "NumPartTotal", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &numBodies);
        H5Aclose(hdf5_attribute);
        H5Sclose(hdf5_dataspace);

        hdf5_dataspace = H5Screate(H5S_SCALAR);
        hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
        H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &t_now);
        H5Aclose(hdf5_attribute);
        H5Sclose(hdf5_dataspace);


        //Header writing done 

        std::vector < size_t >ids(s.n);
        std::iota(ids.begin(), ids.end(), 0);
        sort(ids.begin(), ids.end(),[&s] ( size_t i1, size_t i2)
                {
                return s.part[i1].id < s.part[i2].id;}
            );

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);
        for(unsigned int b = 0; b < s.n; b++)
            data[b] = s.part[ids[b]].pos[0];
        dset_id = H5Dcreate(file_id, "Posx", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);
        for(unsigned int b = 0; b < s.n; b++)
            data[b] = s.part[ids[b]].pos[1];
        dset_id = H5Dcreate(file_id, "Posy", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);
        for(unsigned int b = 0; b < s.n; b++)
            data[b] = s.part[ids[b]].pos[2];
        dset_id = H5Dcreate(file_id, "Posz", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);
        for(unsigned int b = 0; b < s.n; b++)
            data[b] = s.part[ids[b]].vel[0];
        dset_id = H5Dcreate(file_id, "Velx", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);
        for(unsigned int b = 0; b < s.n; b++)
            data[b] = s.part[ids[b]].vel[1];
        dset_id = H5Dcreate(file_id, "Vely", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);
        for(unsigned int b = 0; b < s.n; b++)
            data[b] = s.part[ids[b]].vel[2];
        dset_id = H5Dcreate(file_id, "Velz", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);
        for(unsigned int b = 0; b < s.n; b++)
            data[b] = s.part[ids[b]].mass;
        dset_id = H5Dcreate(file_id, "Mass", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);
        for(unsigned int b = 0; b < s.n; b++)
            data[b] = s.part[ids[b]].pot;
        dset_id = H5Dcreate(file_id, "Potential", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        delete[]data;

        memspace = H5Screate_simple(1, dims, NULL);
        space_id = H5Screate_simple(1, dims, NULL);

        unsigned int *id = new unsigned int[s.n];

        for(unsigned int b = 0; b < s.n; b++)
            id[b] = s.part[ids[b]].id;
        dset_id = H5Dcreate(file_id, "ID", H5T_NATIVE_UINT, space_id, H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_UINT, memspace, space_id, H5P_DEFAULT, id);
        status = H5Sclose(space_id);
        status = H5Sclose(memspace);
        status = H5Dclose(dset_id);

        delete[]id;

        status = H5Gclose(handle);
        status = H5Fclose(file_id);
        if(status < 0)
            printf("Writing snapshot error\n");
    }
    /*-------------------------------------------------------------------------------------
     * Function: std::tuple<int, real_t> load_data(struct falcon::particles::particle_system *s, const size_t restart_val, std::string input_fname)
     *
     * Inputs: struct falcon::particles::particle_system * (pointer to the particle system)
     *         const size_t (if restarting then the number of the snapshot to restart from)
     *         std::string (filename of the initial conditions)
     * Output: std::tuple<int, real_t> (returns the snapshot number from the restart file and the current time which is required for future IO)
     *         
     --------------------------------------------------------------------------------------*/
    std::tuple<int, real_t> load_data(struct falcon::particles::particle_system *s, const size_t restart_val, std::string input_fname) {
        int snapnum {0};
        real_t t_now {0.0};

        if(restart_val == 0) {
            size_t numBodies {};
            std::vector <real_t> array;

            //size_t number_of_lines = 0;
            // for parsing
            std::string line;

            //load the file, should be in ASCII form with [m x y z vx vy vz] form
            std::ifstream file(input_fname);

            if(file.fail()) {
                std::cerr<<"Error opening "<<input_fname<<" "<<std::endl;
                exit(1);
            }

            numBodies = std::count(std::istreambuf_iterator<char>(file), 
                    std::istreambuf_iterator<char>(), '\n');


            /*while(std::getline(file, line))
              ++number_of_lines;
              numBodies = number_of_lines;
              file.close();*/

            std::ifstream fin(input_fname);

            if(fin.fail()) {
                std::cerr<<"Error opening "<<input_fname<<" "<<std::endl;
                exit(1);
            }


            array.resize(numBodies * 7);
            array.assign(std::istream_iterator < real_t >(fin), std::istream_iterator < real_t >());
            fin.close();

            s->n = numBodies;
            s->part = new struct falcon::particles::particle[numBodies];//(struct falcon::particles::particle *) malloc(numBodies * sizeof(struct falcon::particles::particle));
            s->last = &s->part[numBodies - 1];

            for(size_t b = 0; b < numBodies; b++)
            {
                s->part[b].id = b;
                s->part[b].mass = array[b * 7+0];
                s->part[b].pos[0] = array[b * 7+1];
                s->part[b].pos[1] = array[b * 7+2];
                s->part[b].pos[2] = array[b * 7+3];
                s->part[b].vel[0] = array[b * 7+4];
                s->part[b].vel[1] = array[b * 7+5];
                s->part[b].vel[2] = array[b * 7+6];
                /*  for pn terms, set acc_pn to 0 initially and w0 = v0*/
                s->part[b].acc_pn[0] = 0.0;
                s->part[b].acc_pn[1] = 0.0;
                s->part[b].acc_pn[2] = 0.0;
                s->part[b].w[0] = s->part[b].vel[0];
                s->part[b].w[1] = s->part[b].vel[1];
                s->part[b].w[2] = s->part[b].vel[2];
            }

        } else {
            read_hdf5_snapshot(input_fname,&snapnum,&t_now,s);
            for(size_t b=0; b<s->n;b++) {
                /* for pn terms*/
                s->part[b].acc_pn[0] = 0.0;
                s->part[b].acc_pn[1] = 0.0;
                s->part[b].acc_pn[2] = 0.0;
                s->part[b].w[0] = s->part[b].vel[0];
                s->part[b].w[1] = s->part[b].vel[1];
                s->part[b].w[2] = s->part[b].vel[2];
            }
        }  
        return {snapnum, t_now};
    }
}

