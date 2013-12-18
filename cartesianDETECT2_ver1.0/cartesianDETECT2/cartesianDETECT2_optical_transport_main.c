////////////////////////////////////////////////////////////////
//							      //
//	                cartesianDETECT2                      //	
//  (Monte Carlo code for modeling pixelated scintillators)   //    
//							      //
//			MAIN PROGRAM			      //
//							      //
////////////////////////////////////////////////////////////////
//
// 
//
//  Model a detector (with dimensions x = xdetector um, y = ydetector um, z = H um) with multiple rectangular/ cuboid columns.
//  Top surface of detector being isotropic reflector with absorption fraction between [0,1].
//  Bottom surface is an ideal sensor (detects all the optical photons that hit the sensor plane).
//  
//  Model rectangular columns (height=H, side=2a) made of transparent material with specular reflecting walls, top surface being isotropic reflector,
//  bottom surface being an ideal sensor and material has bulk absorption properies (um^-1).
//
//  Roughness of column surface given by coefficient 'beta' (range [0,0.5]).
//  Photon can get lost if it hits the boundaries of detector or goes out of it in x or y direction (x = 0 or xdetector) or (y = 0 or ydetector). 
//  The average path-length of a photon is calculated as the total distance travelled by all
//  photons in the simulation/total # photons generated for all x-rays.
//
//  RANECU RNG used for generating random numbers.
//
//  Input parameters given by command line arguments: 
//	xdetector, ydetector: 			x and y detector dimensions.
//	side, height: 				column side and height.
//	ICgap:					intercolumnar gap.
//	n_C, n_IC:				refractive index for columns and intercolumnar space.
//	top_absfrac, bulk_abscoeff:		Top, bulk absorption coefficient.
//	beta:					roughness coefficient for column surface walls.
//	lbound_x, lbound_y, ubound_x, ubound_y:	lower and upper x, y dimensions for the point response function (PRF).
//	pixelsize:				pixel size for calculating pixel number of the PRF.
//	seed_input:				RNG input seed.
//	input_scintillation_events_filename:	Scintillation events file name.
//	PRF_filename:				PRF filename.
//
//
//
// ****Disclaimer****
//  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in
//  the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection
//  and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software
//  without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the
//  Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
//  parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality,
//  reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory
//  decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are
//  derived from it, and any modified versions bear some notice that they have been modified. 
//
//
//			File:     	cartesianDETECT2_optical_transport_main.c 			
//			Authors:   	Diksha Sharma (CDRH/OSEL/DIAM, US FDA)
//					Aldo Badano   (CDRH/OSEL/DIAM, US FDA)
//			Emails: 	diksha.sharma@fda.hhs.gov, aldo.badano@fda.hhs.gov
//			Date :    	Oct 29, 2013
// 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// header files

// GNU scientific library RNG - used for Poisson sampling of optical photons from energy deposited
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// optical transport kernel
#include "cartesianDETECT2_optical_transport_kernel.c"

// define constants
#define yield 0.06 			// light yield - 60 per keV for Cesium Iodide (CsI)
#define num_primary 1e0			// number of x-ray primaries to be run scintillation events file
#define max_photon_per_EDE 90000000	// maximum number of optical photons that can generate from an energy deposition event	



////////////////////////////////////////////////////////////////////////////
//				MAIN PROGRAM			          //
////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{    
	// variables declaration	

	// command line arguments
	float xdetector, ydetector, side, height, ICgap, n_C, n_IC, top_absfrac, bulk_abscoeff, beta;
	float lbound_x, lbound_y, ubound_x, ubound_y;
	int pixelsize;
        int seed_input;
	char input_fname[80];
	char input_fname1[80];
	char input_fname2[80];

        float dcos[3]={0}; 		// directional cosines
        float normal[3]={0}; 		// surface normal
        float pos[3] = {0}; 		// current position
	float old_pos[3] = {0};  	// old position
	
	int nbytes = 2000*sizeof(struct start_info);	// memory bytes for storing energy deposition events. Assuming there are maximum 2000 per primary

	int seed[2];			// seed used by RANECU once initialized	
	const gsl_rng_type * Tgsl;	// gsl variables
        gsl_rng * rgsl;
	double mu_gsl;	

	clock_t start, end;		// timing variables
	float num_sec;

	char new_line[250];		// used in reading from input file
	char *new_line_ptr;
	char first_char = 'd';
	int rows_read = 0;


	float rr=0.0f, theta=0.0f;	// misc variables
	float r=0.0f, div_numphotons = 0.0f;			
	float norm=0.0f;
	int binx=0, biny=0;
	long int final_counter=0;
	int index=0, jj=0, ii=0;
	int number_particles = 0;
	int result_algo = 0;
	int rownum = 0, colnum = 0;
	int Rndindex = 0;

	printf("\n\t Input arguments: %d\n\n", argc-1);

	// check number of arguments inputted
	if((argc-1) < 19)
	{
		printf("\t Insufficient arguments!! 18 arguments required !! \n\n");
		printf("**************************************************************************************\n");
		printf("\t X, Y  detector dimensions (um). \n");
		printf("\t Side and thickness of column (um). \n");
		printf("\t Intercolumnar gap (um). \n");
		printf("\t Refractive index of column and intercolumnar space. \n");
		printf("\t Top surface absorption, bulk absorption (um^-1), surface roughness coefficient. \n");
		printf("\t Lower and upper x,y bounds for PRF (um). \n");
		printf("\t Pixel pitch (um). \n");
		printf("\t Seed input for RNG. \n");
		printf("\t Scintillation events file name (input) and PRF file name (output). \n");
		printf("\t Radial response file name (output). \n");
		printf("**************************************************************************************\n\n");
		printf("\t PROGRAM EXITING !!!\n\n");
		exit(0);
	}


	// read command line arguments
	xdetector = atof(argv[1]);	// x dimension of detector (in um). x in (0,xdetector)
	ydetector = atof(argv[2]);	// y dimension of detector (in um). y in (0,ydetector)
	side = atof(argv[3]);		// side of cube = 2*a (in um)
	height = atof(argv[4]);		// height of column and thickness of detector (in um). z in range (-H/2, H/2)
	ICgap = atof(argv[5]);		// intercolumnar distance (in um). Assuming it to be equal in x and y direction.
	n_C = atof(argv[6]);		// refractive index of columns
	n_IC = atof(argv[7]);		// refractive index of intercolumnar material
	top_absfrac = atof(argv[8]);	// column's top surface absorption fraction (0.0, 0.5, 0.98)
	bulk_abscoeff = atof(argv[9]);	// column's bulk absorption coefficient (in um^-1) (0.001, 0.1 cm^-1) 
	beta = atof(argv[10]);		// roughness coefficient of column walls
	lbound_x = atof(argv[11]);	// lower bound x of PRF (um)
	lbound_y = atof(argv[12]);	// lower bound y of PRF (um)
	ubound_x = atof(argv[13]);	// upper bound x of PRF (um)
	ubound_y = atof(argv[14]);	// upper bound y of PRF (um)
	pixelsize = atoi(argv[15]);	// pixel size (um)
	seed_input = atoi(argv[16]);	// seed input
	strcpy(input_fname,argv[17]);	// Scintillation events file name - specifying the starting locations and energy deposited
	strcpy(input_fname1,argv[18]);	// PRF filename
	strcpy(input_fname2,argv[19]);  // Radial response filename

	// print out the arguments
	printf("\n**************************************************************************************\n\n");
	printf("\t X detector dimension:\t\t\t %f um \n", xdetector);
	printf("\t Y detector dimension:\t\t\t %f um \n", ydetector);
	printf("\t Side of a column:\t\t\t %f um \n", side);
	printf("\t Thickness of detector:\t\t\t %f um \n", height);
	printf("\t Intercolumnar gap:\t\t\t %f um \n", ICgap);
	printf("\t Refractive index of column:\t\t %f \n", n_C);
	printf("\t Refractive index of intercolumnar space:%f \n", n_IC);
	printf("\t Top surface absorption:\t\t %f \n", top_absfrac);
	printf("\t Bulk absorption coefficient:\t\t %f um^-1 \n", bulk_abscoeff);
	printf("\t Surface roughness coefficient:\t\t %f \n", beta);
	printf("\t Lower x,y bounds for PRF:\t\t %f, %f um \n", lbound_x, lbound_y);
	printf("\t Upper x,y bounds for PRF:\t\t %f, %f um \n", ubound_x, ubound_y);
	printf("\t Pixel pitch:\t\t\t\t %d um \n", pixelsize);
	printf("\t Seed input for RNG:\t\t\t %d \n", seed_input);
	printf("\t Scintillation events file name (input): %s \n", input_fname);
	printf("\t PRF file name (output):\t\t %s \n", input_fname1);
	printf("\t Radial response file name (output):\t\t %s \n", input_fname2);
	printf("\n**************************************************************************************\n\n");

	fflush(stdout);

	// dynamically allocate memory for all the scintillation events (x,y,z,E)
	struct start_info *structa;
	structa = (struct start_info*) malloc(nbytes);
	if( structa == NULL )
		printf("\n Struct start_info array CANNOT BE ALLOCATED !!");


      	// create a generator chosen by the environment variable GSL_RNG_TYPE 
       	gsl_rng_env_setup();

       	Tgsl = gsl_rng_default;
       	rgsl = gsl_rng_alloc (Tgsl);

	gsl_rng_set(rgsl,gsl_rng_default_seed);


	// open files
	FILE *fpout, *fbin, *fradialbin;
	fpout = fopen(input_fname,"rt");
	if(fpout == NULL)
		printf("\n Cannot open deposition events file from penelope for reading!!!");

	fbin = fopen(input_fname1,"w");
	if(fbin == NULL)
		printf("\n Cannot open PRF file for writing!!!");
		
	fradialbin = fopen(input_fname2,"w");
	if(fradialbin == NULL)
		printf("\n Cannot open Radial Bin file for writing!!!");

	fflush(stdout);

	// get current tim for initializing RNG
	time_t Sseconds;
	Sseconds = time (NULL);
	struct timeval Stv;

	// initialize PRF array
	for(binx = 0; binx < bin_arraysizeX; binx++)
	 for(biny = 0; biny < bin_arraysizeY; biny++)
	  {
		bin_matrix[binx][biny] = 0;
	  }


	start = clock();		// start the clock

	for(index = 0; index < num_primary; index++)		// iterate over all primaries
	{

		// reading from file into start_info struct
		first_char = 'd';
		rows_read = 0;


		do{
			new_line_ptr = fgets(new_line, 250, fpout);
			first_char = new_line[0];
			if ( (strcmp(new_line, "\n") != 0) && (first_char != '#') )		// discard empty lines and comments
			{

				// read a line of data (x,y,z,En)
				sscanf(new_line, "%f %f %f %f",&structa[rows_read].str_x, &structa[rows_read].str_y, &structa[rows_read].str_z, &structa[rows_read].str_E);	
			
				mu_gsl = (double)structa[rows_read].str_E * yield;
				structa[rows_read].str_N = gsl_ran_poisson(rgsl,mu_gsl);  // Poisson Sampling - obtain number of optical photons
			
				if(structa[rows_read].str_N > max_photon_per_EDE)
				{
					printf("\n\n Number of optical photons exceeds max. photons defined. Program is EXITING !! \n\n");
					exit(0);
				}
				rows_read++;	// increment number of rows read from the file
			}
		  }while(first_char != '#');	// read until reach "# past hist" line

				
		// reset the variables and global counters
		number_particles = 0;	// number of photons simulated for a primary
		num_generated = 0;	// number of photons generated
		counter=0;		// number detected at the bottom surface of a column
		num_refl_top=0;		// number of particles reflected from top of a column
		num_abs_top=0;		// number of particles absorbed at the top surface
		num_abs_bulk=0;		// number of particles absorbed in the bulk of a column
		num_lost=0;		// number lost at detector boundaries


		for(ii=0; ii<rows_read; ii++)		// iterate over energy deposition events
		{
			// Initialize the RANECU generator in a position far away from the previous history:
			gettimeofday(&Stv,NULL);
			Rndindex = (int)(Sseconds/3600+Stv.tv_usec);	
			init_PRNG(Rndindex, 80000, seed_input, seed);

			int *num_rebound;		// dynamically allocate memory for number of rebounds undergone by each photon
			num_rebound = (int*) malloc(structa[ii].str_N*sizeof(int));

			for(jj=0; jj<structa[ii].str_N; jj++)
				num_rebound[jj] = 0;

			// reset the vectors
			dcos[0]=0.0f; dcos[1]=0.0f; dcos[2]=0.0f;
			normal[0]=0.0f; normal[1]=0.0f; normal[2]=0.0f;
		
			// set starting location of photon
			pos[0] = structa[ii].str_x; pos[1] = structa[ii].str_y; pos[2] = structa[ii].str_z;	
			old_pos[0] = structa[ii].str_x; old_pos[1] = structa[ii].str_y; old_pos[2] = structa[ii].str_z;
	
			// initializing the direction cosines for the first photon
			r = (ranecu(seed) * 2.0f) - 1.0f; // random number between (-1,1)
			 	
			while(fabs(r) <= 0.01f)	
			 {
			   	r = (ranecu(seed) * 2.0f) - 1.0f;  	
			 }

			dcos[2] = r;		// random number between (-1,1)
			rr = sqrt(1.0f-r*r);
			theta=ranecu(seed)*twopipen;
			dcos[0]=rr*cos(theta);
			dcos[1]=rr*sin(theta);

			norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

			if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))	// normalize directional cosines
			 {
				dcos[0] = dcos[0]/norm;
				dcos[1] = dcos[1]/norm;
				dcos[2] = dcos[2]/norm;
			 }


			local_counter=0;				// number of photons that were simulated
			while(local_counter < structa[ii].str_N)	// until all the photons get transported
			 { 
				// reset the flags
				absorbed = 0;		// if absorbed at top surface or lost (1), else 0
				detect = 0;		// if detected at sensor (1), else 0
				bulk_abs = 0;		// if absorbed in bulk (1), else 0

				// set starting location of photon
				pos[0] = structa[ii].str_x; pos[1] = structa[ii].str_y; pos[2] = structa[ii].str_z;	
				old_pos[0] = structa[ii].str_x; old_pos[1] = structa[ii].str_y; old_pos[2] = structa[ii].str_z;
			
				// calculate the row and column number of the cartesian grid where photon starts
				rownum = floor(structa[ii].str_y/(side+ICgap));
				colnum = floor(structa[ii].str_x/(side+ICgap));

				num_generated++;	// number of photons generated
				result_algo = 0;	// 0 - photon still alive; 1 - either detected/absorbed/lost

				while(result_algo == 0)		// call the kernel
				 {
				  	result_algo = algo(normal, old_pos, pos, dcos, num_rebound, seed, &rownum, &colnum, structa[ii], xdetector, ydetector, side, height, ICgap, n_C, n_IC, top_absfrac, bulk_abscoeff, beta, lbound_x, lbound_y, ubound_x, ubound_y, pixelsize);      
				 }
		
			 }	


			// total particles generated for an x-ray
			number_particles = number_particles + structa[ii].str_N;

			free(num_rebound);	// free resources

		}	// rows_read loop ends


		 printf("X-ray number %d \n",index);
		 final_counter = (long int)(final_counter + number_particles);			// total number of particles generated for all x-rays
		 printf("Generated - %d \n", number_particles);
		 printf("\n\t Statistics for columns: \n\n");
		 printf("\t\t Detected - %d \n", counter);
		 printf("\t\t Reflected from top surface - %d\n", num_refl_top);
		 printf("\t\t Absorbed at top surface - %d\n", num_abs_top);
		 printf("\t\t Absorbed in bulk - %d\n\n", num_abs_bulk);
		 printf("\t Statistics for intercolumnar space: \n\n");
		 printf("\t\t Detected - %d\n", global_num_detect);
		 printf("\t\t Reflected from top surface - %d\n", global_num_refl_top);
		 printf("\t\t Absorbed at top surface - %d\n", global_num_abs_top);
		 printf("\t\t Lost at detector boundary - %d \n\n\n\n", num_lost);



	}	// index loop ends

	// write PRF to file

	div_numphotons = 1.0f/num_primary;

	for(binx = 0; binx < bin_arraysizeX; binx++)
	{
	 for(biny = 0; biny < bin_arraysizeY; biny++)
	  {
		fprintf(fbin,"%f\n", (bin_matrix[binx][biny]*div_numphotons) );		// average number of photons detected at each pixel
	  }
	  fprintf(fbin,"\n");
	}
	
	// write radial bin to file
	for(binx = 0; binx < 21 ; binx++)
	  {
		fprintf(fradialbin,"%d %f %f\n", binx, radialbin[binx], start_pos[2] );		// average number of photons detected at each pixel
	  }
	
	

	end = clock();		// end the clock


	printf("*********************************************************************************\n");
	printf("\n  Average path length of photons = %f um \n", photon_distance/number_particles);

	num_sec = (((float)end - start)/CLOCKS_PER_SEC);
	printf("\n  Total number of photons simulated for %e primaries = %ld \n", num_primary, final_counter);
	printf("\n  Total time %f sec --  Speed = %f photons per sec\n\n", num_sec, final_counter/num_sec );

	printf("\t SIMULATION ENDED !!\n\n");

	free(structa);
	fclose(fpout);
	fclose(fbin);
	fclose(fradialbin);

     
        return 0;
}	// main() ends


