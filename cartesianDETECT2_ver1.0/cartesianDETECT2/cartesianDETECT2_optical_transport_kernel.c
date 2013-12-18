////////////////////////////////////////////////////////////////
//							      //
//	                cartesianDETECT2                      //	
//  (Monte Carlo code for modeling pixelated scintillators)   //    
//							      //
//			KERNEL PROGRAM			      //
//							      //
////////////////////////////////////////////////////////////////
//
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
//			File:        cartesianDETECT2_optical_transport_kernel.c 
//			Included in: cartesianDETECT2_optical_transport_main.c			
//			Authors:     Diksha Sharma (CDRH/OSEL/DIAM, US FDA)
//				     Aldo Badano   (CDRH/OSEL/DIAM, US FDA)
//			Emails:      diksha.sharma@fda.hhs.gov, aldo.badano@fda.hhs.gov
//			Date :       Oct 29, 2013
//
// 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////
//
//      Header libraries
//
/////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>


/////////////////////////////////////////
//
//       Constants
//
/////////////////////////////////////////

#define pi 3.14159265		// pi
#define twopipen 6.283185308	// 2*pi
#define epsilon 8.1929093e-6	// very small number for floating point comparisons
#define bin_arraysizeX 501	// ( (ubound_x - lbound_x)/pixelsize ) + 1
#define bin_arraysizeY 501	// ( (ubound_y - lbound_y)/pixelsize ) + 1


/////////////////////////////////////////
//
//     Structure for storing energy deposition events
//
/////////////////////////////////////////

struct start_info
{
	float str_x;	// x
	float str_y;	// y
	float str_z;	// z
	float str_E;	// energy
	int str_N;	// number of photons to be simulated
};


/////////////////////////////////////////
//
//       Function declarations
//
/////////////////////////////////////////

// moves photon from its generation to detection at the sensor plane or absorption or lost.
int algo(float *normal, float *old_pos, float *pos, float *dcos, int *num_rebound, int* seed, int *rowi, int *colj, struct start_info info, float xdetector, float ydetector, float side, float height, float g, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, float lbound_x, float lbound_y, float ubound_x, float ubound_y, int pixelsize);

// photon within a column. calculate if it gets absorbed or specularly reflects within the column.
int specular_refl(float *old_pos, float *pos, float *dcos, int* seed, int *rowi, int *colj, float bulk_abscoeff, float a, float H, float g, float xdetector, float ydetector);

// photon within a column. find distance to next position and move it.
float dist_to_surface(float *pos, float *dcos, int *rowi, int *colj, float a, float H, float g, float xdetector, float ydetector);

// photon within/between columns. calculate if it gets reflected or transmitted.
int mirror(float *normal, float *old_pos, float *pos, float *dcos, int* seed, int *rowi, int *colj, float xdetector, float ydetector, float a, float H, float g, float n1, float n2, float top_absfrac, int rowi_max, int colj_max, float beta);

// transmit photon to another column. calculates new position.
int transmit(float* pos, float* dcos, float* normal, int* seed, int *rowi, int *colj, float xdetector, float ydetector, float a, float H, float g, float top_absfrac, int rowi_max, int colj_max);

// calculate directional cosines of reflected/refracted vector.
void next_dir_cos(float* dcos, float* normal, float refl_theta, float trans_theta, int flag_ref);

// calculate rough normal depending on value of 'beta'.
void RoughSurface(float* normal, int* seed, float beta);

// determine if photon gets detected at the sensor plane.
int detection(float *pos, int num_rebound, float H, float lbound_x, float lbound_y, float ubound_x, float ubound_y, int pixelsize, struct start_info info);

// calculate dot product of two vectors to give cosine of angle between them.
float dot_product(float *aa, float *b);

// find minimum of the three floats
float min(float a1, float a2);
float min3(float a1, float a2, float a3);

// RANECU random number generator
void init_PRNG(int history_batch, int histories_per_thread, int seed_input, int* seed);
int abMODm(int m_par, int a_par, int s_par);
float ranecu(int* seed);


/////////////////////////////////////////
//
//       Global variables
//
/////////////////////////////////////////

int num_generated=0;		// number of photons generated	
int counter=0;			// number detected at the bottom surface of a column
int local_counter=0;		// number of photons transported - either detected / absorbed at top or bulk / lost
int num_refl_top=0;		// number of particles reflected from top of a column
int num_abs_top=0;		// total number of particles absorbed at the top surface of a column (using 'top_absfrac')
int num_abs_bulk=0;		// total number of particles absorbed in a column (using 'bulk_abscoeff')
int num_lost=0;			// total number of particles lost at the detector boundary in x/y direction
int bin_matrix[bin_arraysizeX][bin_arraysizeY]={{0}};	// PRF: number of photons detected from lbound to ubound_x/y
float radialbin[21]={0.0f};
float start_pos[3] = {0.0f};

int absorbed=0;			// flag for top surface absorption of a column
int detect=0, bulk_abs=0;	// flag for particle detected/absorbed in the material in a column
int global_num_detect=0;	// number of particles detected at the bottom of detector in the intercolumnar space
int global_num_abs_top=0;	// number of particles absorbed at the top of detector in the intercolumnar space
int global_num_refl_top=0;	// number of particles reflected from the top of detector in the intercolumnar space
float photon_distance=0.0f;     // total distance travelled by all the photons


/////////////////////////////////////////
//
//    Functions definition
//
/////////////////////////////////////////

int algo(float *normal, float *old_pos, float *pos, float *dcos, int *num_rebound, int* seed, int *rowi, int *colj, struct start_info info, float xdetector, float ydetector, float side, float H, float g, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, float lbound_x, float lbound_y, float ubound_x, float ubound_y, int pixelsize)
{

        int myresult = 0;		// flag returned by algo to indicate if the photon is still alive (0) else (1).
	float a = (float)(side/2.0f);	// side/2
	float rr=0.0f, theta=0.0f;	// used for calculating initial dir. cosines
	float r=0.0f;	
	float norm = 0.0f;		// used for normalizing

	
	// compute the maximum number of rows and columns in the detector
	int rowi_max = floor((ydetector-2*a)/(2*a+g));	// rowi in (0,rowi_max)
	int colj_max = floor((xdetector-2*a)/(2*a+g));	// colj in (0,colj_max)

	// Start the optical transport

	if(absorbed == 0)        	// if not absorbed at top or lost
	 {
		bulk_abs = specular_refl(old_pos, pos, dcos, seed, rowi, colj, bulk_abscoeff, a, H, g, xdetector, ydetector);

		if(bulk_abs == 0)	// if not absorbed in the bulk
		{
			detect = detection(pos, num_rebound[local_counter], H, lbound_x, lbound_y, ubound_x, ubound_y, pixelsize, info);
		}
	 }

 
        if( (detect == 1) || (absorbed == 1) || (bulk_abs == 1) ) 	// particle terminated
         {

		local_counter++;		// increment number of photons transported

		// generate initial directional cosines for next photon
		r = (ranecu(seed) * 2.0f) - 1.0f; 
 	
		while(fabs(r) <= 0.01)	
		 {
	   		r = (ranecu(seed) * 2.0f) - 1.0f;  	
	 	 }


        	dcos[2] = r;		// random number between (-1,1)
	        rr = sqrt(1.0-r*r);
	        theta=ranecu(seed)*twopipen;
	        dcos[0]=rr*cos(theta);
	        dcos[1]=rr*sin(theta);

		norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

		if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))	// normalize
		 {
			dcos[0] = dcos[0]/norm;
			dcos[1] = dcos[1]/norm;
			dcos[2] = dcos[2]/norm;
		 }

		// set starting location of photon
		pos[0] = info.str_x; pos[1] = info.str_y; pos[2] = info.str_z;
		old_pos[0] = info.str_x; old_pos[1] = info.str_y; old_pos[2] = info.str_z;
		start_pos[0] = info.str_x; start_pos[1] = info.str_y; start_pos[2] = info.str_z;

		normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 0.0f;	// initialize the surface normal
		
		RoughSurface(normal, seed, beta);	// perturb to get rough normal according to 'beta'

		// reset the flags
		absorbed = 0;		// if absorbed at top surface or lost (1), else 0
		detect = 0;		// if detected at sensor (1), else 0
		bulk_abs = 0;		// if absorbed in bulk (1), else 0

		myresult = 1;	

         }
	else if( (detect == 0) && (absorbed == 0) && (bulk_abs == 0) && (fabs(dcos[2] - 0.0f) < epsilon) )  // checking for particle going back and forth with dcos(z)=0 - TRAPPED
	 {
		// we need to kill the particle and generate a new one instead - do not increment the counter

		// re-initialize all the arrays
                r = (ranecu(seed) * 2.0f) - 1.0f;

                while(fabs(r) <= 0.01)	
		 {
	   		r = (ranecu(seed) * 2.0f) - 1.0f;  	
	 	 }


                dcos[2] = r;            // random number between (-1,1)
                rr = sqrt(1.0-r*r);
                theta=ranecu(seed)*twopipen;
                dcos[0]=rr*cos(theta);
                dcos[1]=rr*sin(theta);

		

		norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

		if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
		 {
			dcos[0] = dcos[0]/norm;
			dcos[1] = dcos[1]/norm;
			dcos[2] = dcos[2]/norm;
		 }

		myresult = 0;
	 }
	else	// photon still alive. increment the number of rebounds and call mirror()
	 {
		num_rebound[local_counter]++;
            	absorbed = mirror(normal, old_pos, pos, dcos, seed, rowi, colj, xdetector, ydetector, a, H, g, n1, n2, top_absfrac, rowi_max, colj_max, beta);

		myresult = 0;
	 }


 return myresult;	// if particle alive (0), else (1)
}


// photon within a column. calculate if it gets absorbed or specularly reflects within the column.
int specular_refl(float *old_pos, float *pos, float *dcos, int* seed, int *rowi, int *colj, float bulk_abscoeff, float a, float H, float g, float xdetector, float ydetector)
{

	float dsurf = 999.0f;	// distance to surface within column
	float dabs = 999.0f;	// distance to absorption
	float d = 999.0f;	// final distance
	int flag_bulkabs = 0;	// flag for bulk absorption (1) - photon absorbed

	old_pos[0] = pos[0];
	old_pos[1] = pos[1];
	old_pos[2] = pos[2];

	dsurf = dist_to_surface(pos, dcos, rowi, colj, a, H, g, xdetector, ydetector);		// distance to surface
	dabs = (-1.0f/bulk_abscoeff) * log(ranecu(seed));					// distance to absorption

	if (fabs(dsurf-(-99.0f)) < epsilon)	// photon lost because went out of detector boundary in dist_to_surface()
	{
		flag_bulkabs = 1;
	}
	else if (dsurf < dabs)			// photon reflects
	 {		
		d = dsurf;
		flag_bulkabs = 0;
	 }
	else if (dsurf >= dabs)			// photon absorbed
	 {
		d = dabs;
		flag_bulkabs = 1;

		num_abs_bulk++;
	 }

   return flag_bulkabs;		// flag  - photon absorbed (1), else still alive (0)
}

// photon within a column. find distance to surface and move it.
float dist_to_surface(float *pos, float *dcos, int *rowi, int *colj, float a, float H, float g, float xdetector, float ydetector)
{

	float d=0.0f;
	float dx = 0.0f, dy = 0.0f, dz = 0.0f;		// distance to x, y, z plane
	int f1=0, f2=0, f3=0;
	float temp_pos[3] = {0.0f};
	float vsmall_num = 0.0001f;

	temp_pos[0] = pos[0];
	temp_pos[1] = pos[1];
	temp_pos[2] = pos[2];
	
	// check if particle travel straight in +/- z direction. 
	// particle cannot travel straight in +/- x or y direction (will keep rebounding back and forth). This case been checked for in the algo().
	if((fabs(dcos[0] - 0.0f) < epsilon) && (fabs(dcos[1] - 0.0f) < epsilon) && (fabs(dcos[2] - 1.0f) < epsilon))
	{
		d = H/2.0f - pos[2];
		pos[2] = H/2.0f;

		pos[0] = temp_pos[0] + d*dcos[0];
		pos[1] = temp_pos[1] + d*dcos[1];
	}
	else if((fabs(dcos[0] - 0.0f) < epsilon) && (fabs(dcos[1] - 0.0f) < epsilon) && (fabs(dcos[2] - (-1.0f)) < epsilon))
	{
		d = fabs(-H/2.0f - pos[2]);
		pos[2] = -H/2.0f;
	}
	// else the particle is going to hit any plane at an angle
	else
	{
		// calculating distance to the sides of column

		if (dcos[0] > 0.0f)
		{
			dx = fabs(((2*a*(*colj+1)+((*colj)+1)*g)-pos[0])/dcos[0]);	// +x plane
			f1 = 1;
		}
		else if (dcos[0] < 0.0f)
		{
			dx = fabs(((2*a*(*colj)+((*colj)+1)*g)-pos[0])/dcos[0]);		
			f1 = -1;
		}

		if (dcos[1] > 0.0f)
		{
			dy = fabs(((2*a*(*rowi+1)+((*rowi)+1)*g)-pos[1])/dcos[1]);	// +y plane
			f2 = 1;
		}
		else if (dcos[1] < 0.0f)
		{
			dy = fabs(((2*a*(*rowi)+((*rowi)+1)*g)-pos[1])/dcos[1]);
			f2 = -1;
		}

		if (dcos[2] > 0.0f)
		{
			dz = fabs((H/2.0f-pos[2])/dcos[2]); 			// +z plane
			f3 = 1;
		}
		else if (dcos[2] < 0.0f)
		{
			dz = fabs((-H/2.0f-pos[2])/dcos[2]);
			f3 = -1;
		}
	

		// find the min of these three distances to get where particle hist the column
		if ( (dx > vsmall_num) && (dy > vsmall_num) && (dz > vsmall_num) )
			d = min3(dx,dy,dz);
		else if ( (dx > vsmall_num) && (dy > vsmall_num) && (dz <= vsmall_num) )
			d = min(dx,dy);
		else if ( (dx > vsmall_num) && (dy <= vsmall_num) && (dz > vsmall_num) )
			d = min(dx,dz);
		else if ( (dx <= vsmall_num) && (dy > vsmall_num) && (dz > vsmall_num) )
			d = min(dy,dz);


		if( fabs(d-dx) < epsilon )
		{
			if(f1 == 1)	// hits x=(2*a*(*colj+1)+(*colj+1)*g) plane. 'g' is the inter columnar gap. a = side/2.
				pos[0] = (2*a*(*colj+1)+((*colj)+1)*g);
			else if(f1 == -1)
				pos[0] = (2*a*(*colj)+((*colj)+1)*g);

			pos[1] = temp_pos[1] + d*dcos[1];
			pos[2] = temp_pos[2] + d*dcos[2];
		}
		else if(fabs(d-dy) < epsilon)
		{
			if(f2 == 1)	// hits y=(2*a*(*rowi+1)+(*rowi+1)*g) plane
				pos[1] = (2*a*(*rowi+1)+((*rowi)+1)*g);
			else if(f2 == -1)
				pos[1] = (2*a*(*rowi)+((*rowi)+1)*g);

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[2] = temp_pos[2] + d*dcos[2];


		}
		else if(fabs(d-dz) < epsilon)
		{
			if(f3 == 1)	// hits z=H/2 plane
				pos[2] = H/2.0f;
			else if(f3 == -1)
				pos[2] = -H/2.0f;

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];

		}

	// condition to check that pos is within detector boundaries - if true, particle LOST
	if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) || (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
		{
			d = -99.0f;
			num_lost++;
		}
	else
		photon_distance = photon_distance + d;		// add distance travelled to global variable


		
	} // else ends

 return d;
}


// photon within/between columns. calculate if it gets reflected ot transmitted.
int mirror(float *normal, float *old_pos, float *pos, float *dcos, int* seed, int *rowi, int *colj, float xdetector, float ydetector, float a, float H, float g, float n1, float n2, float top_absfrac, int rowi_max, int colj_max, float beta)
{
    float angle=0.0f;
    float dcos_temp[3] = {0.0f};
    int i=0;
    float Pr = 0.0f, Pt = 0.0f;		// Prob. of reflection and transmission
    float theta1 = 0.0f, theta2 = 0.0f;	// theta1: angle between surface normal and incoming vector. theta2: between normal and refracted vector.
    int trans_flag = 0.0f;		// flag from transmit()
    int flag_abs = 0;		// flag to indicate if particle got absorbed at top surface or exited during the transmission to another boundary
    int flag_call_transmit = 1;	// flag to indicate if the particle is going to move within a column (flag = 0) [call specular_refl()] 
				// or between columns (flag = 1) [call transmit()]
    float temp_norm = 0.0f;
    float old_normal[3] = {0.0f};
    float old_dcos[3] = {0.0f};
    float angle_oldN_R = 0.0f;
    int reperturb_ctr = 0;


        // determine the coordinates of surface normal
	if ( (fabs(pos[2] - (float)(H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	// reached top surface and dir. cosine in z-direction is positive
	{
	
		// top surface absorption - using absorption coefficient 'top_absfrac'
		if ( (top_absfrac > 0.0f) && (ranecu(seed) < top_absfrac) )	// particle gets absorbed
		{
			flag_abs = 1;
			num_abs_top++;
		}
		else
		{
			num_refl_top++;

			// specular reflector		
			normal[0] = 0.0f;
			normal[1] = 0.0f;
			normal[2] = -1.0f;
	
			// -dcos -> inverted the incident vector to get the smaller angle, else would have to do angle = 180-angle
			dcos_temp[0] = -dcos[0];
			dcos_temp[1] = -dcos[1];
			dcos_temp[2] = -dcos[2];

			angle = dot_product(dcos_temp,normal);	// cosine of angle between incident in opposite direction and normal
			for(i=0;i<3;i++)
			 {
			  	dcos[i] = 2.0f*angle*normal[i] + dcos[i];  // specular ray
			 }

	
			flag_abs = 0;
		}

	}	
	else 	// compute the normal and check if gets reflected or transmitted
	{	
	prpt:

		// within the column
		if( (fabs(pos[0] - (2*a*(*colj+1)+((*colj)+1)*g)) < epsilon) && (dcos[0] > 0.0f) ) // hits x=a plane and dir.cos. in x-dir is positive
		{
			normal[0] = -1.0f;
			normal[1] = 0;
			normal[2] = 0;

			flag_call_transmit = 1;		// to indicate that photon is currrently within a column
			flag_abs = 0;
		}
		else if( (fabs(pos[0] - (2*a*(*colj)+((*colj)+1)*g)) < epsilon) && (dcos[0] < 0.0f) ) // hits x=-a plane and dir.cos. in x-dir is negative
		{
			normal[0] = 1.0f;
			normal[1] = 0;
			normal[2] = 0;

			flag_call_transmit = 1;
			flag_abs = 0;
		}
		else if( (fabs(pos[1] - (2*a*(*rowi+1)+(*rowi+1)*g)) < epsilon) && (dcos[1] > 0.0f) ) // hits y=a plane and dir.cos. in y-dir is positive
		{
			normal[0] = 0;
			normal[1] = -1.0f;
			normal[2] = 0;

			flag_call_transmit = 1;
			flag_abs = 0;
		}
		else if( (fabs(pos[1] - (2*a*(*rowi)+((*rowi)+1)*g)) < epsilon) && (dcos[1] < 0.0f) ) // hits y=-a plane and dir.cos. in y-dir is negative
		{
			normal[0] = 0;
			normal[1] = 1.0f;
			normal[2] = 0;

			flag_call_transmit = 1;
			flag_abs = 0;
		}
		// outside the column
                else if( (fabs(pos[0] - (2*a*(*colj+1)+((*colj)+1)*g)) < epsilon) && (dcos[0] < 0.0f) )  // hits x=a plane and dir.cos. in x-dir is negative
                {
                        normal[0] = 1.0f;
                        normal[1] = 0;
                        normal[2] = 0;

			flag_call_transmit = 0;		// to indicate that the photon is currently between columns and has not entered any column yet
                        flag_abs = 0;
                }
                else if( (fabs(pos[0] - (2*a*(*colj)+((*colj)+1)*g)) < epsilon) && (dcos[0] > 0.0f) ) // hits x=-a plane and dir.cos. in x-dir is positive
                {
                        normal[0] = -1.0f;
                        normal[1] = 0;
                        normal[2] = 0;
			
			flag_call_transmit = 0;
                        flag_abs = 0;
                }
                else if( (fabs(pos[1] - (2*a*(*rowi+1)+((*rowi)+1)*g)) < epsilon) && (dcos[1] < 0.0f) ) // hits y=a plane and dir.cos. in y-dir is negative
                {
                        normal[0] = 0;
                        normal[1] = 1.0f;
                        normal[2] = 0;

			flag_call_transmit = 0;
                        flag_abs = 0;
                }
                else if( (fabs(pos[1] - (2*a*(*rowi)+((*rowi)+1)*g)) < epsilon) && (dcos[1] > 0.0f) ) // hits y=-a plane and dir.cos. in y-dir is positive
                {
                        normal[0] = 0;
                        normal[1] = -1.0f;
                        normal[2] = 0;

			flag_call_transmit = 0;
                        flag_abs = 0;
                }

		// Using Snell's law, calculate theta1 (angle between normal and reflected) and theta2 (angle between normal and transmitted)
		// -dcos -> inverted the incident vector to get the smaller angle, else would have to do angle = 180-angle
		dcos_temp[0] = -dcos[0];
		dcos_temp[1] = -dcos[1];
		dcos_temp[2] = -dcos[2];

		old_normal[0] = normal[0];
		old_normal[1] = normal[1];
		old_normal[2] = normal[2];
	
		old_dcos[0] = dcos[0];
		old_dcos[1] = dcos[1];
		old_dcos[2] = dcos[2];

	reperturb:
		normal[0] = old_normal[0];
		normal[1] = old_normal[1];
		normal[2] = old_normal[2];

		dcos[0] = old_dcos[0];
		dcos[1] = old_dcos[1];
		dcos[2] = old_dcos[2];

		dcos_temp[0] = -dcos[0];
		dcos_temp[1] = -dcos[1];
		dcos_temp[2] = -dcos[2];

		RoughSurface(normal, seed, beta);	// new normal for rough surface

	no_perturbation:
		theta1 = dot_product(dcos_temp, normal);// cosine of angle between incident in opposite direction and normal (in radians)

		if ( (theta1 > 1.0f) || (theta1 < 0.0f) )	// if incidence angle > 1.57 radian or < 0 radian, then recalculate normal
			goto prpt;
		else
			theta1 = acosf(theta1);

		// check for conditions where photon can only reflect
		if (flag_call_transmit == 1)		// only valid when photon within the column and can transmit outside the column. asin(n1/n2) -> nan
		{
			if (theta1 > asin(n2/n1))	// critical angle condition for TIR
			{
				Pr = 1.0f;		// TIR occurs
				Pt = 0.0f;
			}
               		else if ( theta1 < epsilon ) // theta1 ~= 0, then always reflect
                	{
                        	theta1 = 0.00042;       // make theta1 a very smal number, to avoid getting nan probabilities
	                        theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians

        	                // Using Fresnel's law, compute probability of reflection and transmission 
                	        Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
                        	Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
	                }
        	        else    // the ray will transmit
                	{
                        	if (flag_call_transmit == 1)
                                	theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians
	                        else if (flag_call_transmit == 0)
        	                        theta2 = asinf((float)(n2/n1)*sin(theta1));

                	        // Using Fresnel's law, compute probability of reflection and transmission 
                        	Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
	                        Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
        	        }


		}
		else if (flag_call_transmit == 0)
		{		
			if ( theta1 < epsilon )	// theta1 ~= 0, then always reflect
			{
				theta1 = 0.00042;	// make theta1 a very smal number, to avoid getting nan probabilities
				theta2 = asinf((float)(n1/n2)*sin(theta1)); 	// refracted/transmitted angle in radians
	
				// Using Fresnel's law, compute probability of reflection and transmission 
				Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
				Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
			}
			else	// the ray will transmit
			{
				if (flag_call_transmit == 1)
					theta2 = asinf((float)(n1/n2)*sin(theta1)); 	// refracted/transmitted angle in radians
				else if (flag_call_transmit == 0)
					theta2 = asinf((float)(n2/n1)*sin(theta1));

				// Using Fresnel's law, compute probability of reflection and transmission 
				Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
				Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
			}
		}

		// normalize Pr and Pt
		temp_norm = Pr + Pt;
		Pr = Pr/temp_norm;
		Pt = Pt/temp_norm;

		// check for Pr+Pt = 1
		if ( (fabs(Pr+Pt)-1.0f) < epsilon)
			;
		else
		 {
			printf("Pr + Pt IS NOT EQUAL TO ONE!!!!!!!!!!!!!! theta1 %f theta2 %f Pr %f, Pt %f dcos %f %f %f normal %f %f %f flag %d\n\n", theta1, theta2, Pr, Pt, dcos[0], dcos[1], dcos[2], normal[0], normal[1], normal[2], flag_call_transmit);
			exit(0);
		}

		if(ranecu(seed) < Pr)	// reflection
		{
			next_dir_cos(dcos, normal, theta1, theta2, 0);

			// condition to check that reflected vector is within 90 degrees from original normal
			angle_oldN_R = dot_product(old_normal, dcos);
			angle_oldN_R = acosf(angle_oldN_R);


			if (angle_oldN_R > 1.57f) // > 90 degrees, reperturb the normal
			{
				reperturb_ctr++;
				
				if(reperturb_ctr < 4)
					goto reperturb;
				else
				{
					normal[0] = old_normal[0];
					normal[1] = old_normal[1];
					normal[2] = old_normal[2];

					dcos[0] = old_dcos[0];
					dcos[1] = old_dcos[1];
					dcos[2] = old_dcos[2];

					dcos_temp[0] = -dcos[0];
					dcos_temp[1] = -dcos[1];
					dcos_temp[2] = -dcos[2];

					reperturb_ctr = 0;
				
					goto no_perturbation;
				}

			}

			if (flag_call_transmit == 0)	// it is reflecting between columns, so need to calculate distance using transmit()
			{
				trans_flag = transmit(pos, dcos, normal, seed, rowi, colj, xdetector, ydetector, a, H, g, top_absfrac, rowi_max, colj_max);

				if (trans_flag == 1)	// photon exited
					flag_abs = 1;
				else if (trans_flag == 0)
					goto prpt;				
			}

		}
		else			// transmission
		{
			next_dir_cos(dcos, normal, theta1, theta2, 1);

			if (flag_call_transmit == 1)	// photon travels between columns
			{
				trans_flag = transmit(pos, dcos, normal, seed, rowi, colj, xdetector, ydetector, a, H, g, top_absfrac, rowi_max, colj_max);

				if (trans_flag == 1)	// particle exited
					flag_abs = 1;
				else if (trans_flag == 0)	// hits a column
					goto prpt;	// check again to see if it gets reflected or transmitted
			}
		}

	} //else ends


   return flag_abs;

}


// transmit photon to another column. calculates new position.
int transmit(float* pos, float* dcos, float* normal, int* seed, int *rowi, int *colj, float xdetector, float ydetector, float a, float H, float g, float top_absfrac, int rowi_max, int colj_max)
{
	int myplane = 0;		// flag to indicate the plane from where photon emits
	float temp_pos[3] = {0.0f};
	float d_transmit = 0.0f;	// distance to next surface	
	int particle_exit = 0;		// flag to indicate if photon enters another column or gets lost/detected/absorbed
	int index = 0;
	float newangle = 0.0f;
	float newdcos_temp[3] = {0.0f};
	float top_pos[3] = {0.0f};
	int i=0;

	if ( fabs( pos[0] - (2*a*((*colj)+1) + ((*colj)+1)*g) ) < epsilon )			// emits from x = a
		myplane = 1;
	else if ( fabs( pos[0] - (2*a*(*colj) + ((*colj)+1)*g) ) < epsilon )			// emits from x = -a
		myplane = -1;
	else if ( fabs( pos[1] - (2*a*((*rowi)+1) + ((*rowi)+1)*g) ) < epsilon )		// emits from y = a
		myplane = 2;
	else if ( fabs( pos[1] - (2*a*(*rowi) + ((*rowi)+1)*g) ) < epsilon )			// emits from y = -a
		myplane = -2;

	temp_pos[0] = pos[0];
	temp_pos[1] = pos[1];
	temp_pos[2] = pos[2];

	// compute the new position of the photon. 
	// Assuming the intercolumnar distance 'g' to be negligible, thus when photon emits from any plane it will enter one of the columns in the opposite col/row.

	if (myplane == 1)
	{
		if(*colj != colj_max)
			pos[0] = temp_pos[0] + g;					// hits pos[0]+g plane
		else
			pos[0] = xdetector;						// hits detector boundary

		d_transmit = (pos[0] - temp_pos[0])/dcos[0];

		pos[1] = temp_pos[1] + d_transmit*dcos[1];
		pos[2] = temp_pos[2] + d_transmit*dcos[2];
	}
	else if (myplane == -1)
	{
		pos[0] = temp_pos[0] - g;
		d_transmit = (pos[0] - temp_pos[0])/dcos[0];

		pos[1] = temp_pos[1] + d_transmit*dcos[1];
		pos[2] = temp_pos[2] + d_transmit*dcos[2];
	}
	else if (myplane == 2)
	{
		if(*rowi != rowi_max)
			pos[1] = temp_pos[1] + g;
		else
			pos[1] = ydetector;

		d_transmit = (pos[1] - temp_pos[1])/dcos[1];

		pos[0] = temp_pos[0] + d_transmit*dcos[0];
		pos[2] = temp_pos[2] + d_transmit*dcos[2];
	}
	else if (myplane == -2)
	{
		pos[1] = temp_pos[1] - g;
		d_transmit = (pos[1] - temp_pos[1])/dcos[1];

		pos[0] = temp_pos[0] + d_transmit*dcos[0];
		pos[2] = temp_pos[2] + d_transmit*dcos[2];
	}

	// condition to check that pos is within detector boundaries - if true, particle LOST
	if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) || (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
		{
			d_transmit = -99.0f;
			num_lost++;
			particle_exit = 1;
			goto exitnow;
		}
	else
		photon_distance = photon_distance + d_transmit;		// add distance travelled to global variable

	
	// check if the photon enters another column or got lost (hit detector side)/ reflected (detector top)/ detected (detector bottom)
	
	// hit side of detector?
	if ( ( fabs(pos[0]-0.0f) < epsilon ) || ( fabs(pos[0]-xdetector) < epsilon ) )		// hits yz plane, gets lost
	{
		num_lost++;
		particle_exit = 1;
		goto exitnow;

	}
	else if ( ( fabs(pos[1]-0.0f) < epsilon ) || ( fabs(pos[1]-ydetector) < epsilon ) )	// hits xz plane, gets lost
	{
		num_lost++;
		particle_exit = 1;
		goto exitnow;
	}

	// hit top?
	if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )				// gets specularly reflected or absorbed
	{
		normal[0] = 0.0f;
		normal[1] = 0.0f;
		normal[2] = -1.0f;

		// top surface absorption - using absorption coefficient 'top_absfrac'
		if ( (top_absfrac > 0.0f) && (ranecu(seed) < top_absfrac) )			// particle gets absorbed
		{
			global_num_abs_top++;
			particle_exit = 1;
			goto exitnow;
		}
		else
		{
			global_num_refl_top++;

			newdcos_temp[0] = -dcos[0];
			newdcos_temp[1] = -dcos[1];
			newdcos_temp[2] = -dcos[2];
			
			
			newangle = dot_product(newdcos_temp,normal);  // cosine of angle between incident in opposite direction and normal
			for(i=0;i<3;i++)
			 {
			  	dcos[i] = 2.0f*newangle*normal[i] + dcos[i];  // specular ray
			 }

			// store current position in another variable
			top_pos[0] = pos[0];
			top_pos[1] = pos[1];
			top_pos[2] = pos[2];			

			// to simplify the computation for determing next surface, we assume it will hit the adjacent row/col based on plane it emits
			// similar to how we calculated position above in 'myplane' loops

			if (myplane == 1)		// emits from x=a plane, so assume will hit x=-a plane
			{
				if(*colj != colj_max)
					pos[0] = temp_pos[0] + g; // temp_pos[0] -> orig. x position before it hits the top surface				
				else
					pos[0] = xdetector;

				d_transmit = (pos[0] - top_pos[0])/dcos[0]; // d_transmit will be distance between top position and new position

				pos[1] = top_pos[1] + d_transmit*dcos[1];
				pos[2] = top_pos[2] + d_transmit*dcos[2];

			}
			else if (myplane == -1)
			{
				pos[0] = temp_pos[0] - g;
				d_transmit = (pos[0] - top_pos[0])/dcos[0];

				pos[1] = top_pos[1] + d_transmit*dcos[1];
				pos[2] = top_pos[2] + d_transmit*dcos[2];

			}
			else if (myplane == 2)
			{
				if(*rowi != rowi_max)
					pos[1] = temp_pos[1] + g;
				else
					pos[1] = ydetector;

				d_transmit = (pos[1] - top_pos[1])/dcos[1];

				pos[0] = top_pos[0] + d_transmit*dcos[0];
				pos[2] = top_pos[2] + d_transmit*dcos[2];

			}
			else if (myplane == -2)
			{
				pos[1] = temp_pos[1] - g;
				d_transmit = (pos[1] - top_pos[1])/dcos[1];

				pos[0] = top_pos[0] + d_transmit*dcos[0];
				pos[2] = top_pos[2] + d_transmit*dcos[2];

			}

			particle_exit = 0;
		}
	}	// hit top ends
	

	// hit bottom? z of detector can be in the range (-H/2, H/2).
	if ( fabs(pos[2] - (-H/2.0f)) < epsilon )	// gets detected
	{
		global_num_detect++;
		particle_exit = 1;
		goto exitnow;
	}

	
	// if it enters a column, find the column row and col index
	if (myplane == 1)	// emits plane x = a
	{
		if ( (pos[1] >= (2*a*(*rowi)+((*rowi)+1)*g)) && (pos[1] <= (2*a*((*rowi)+1)+((*rowi)+1)*g)) )	// enters adjacent column on right
		{
			(*colj) = (*colj)+1;	
			particle_exit = 0;
		}	
		else if ( pos[1] > (2*a*((*rowi)+1)+((*rowi)+1)*g) )	// moves upward +y direction
		{
			if (*rowi == rowi_max)	// goes outside last row - taken as getting to the detector boundary and lost
			{
				num_lost++;
				particle_exit = 1;
				goto exitnow;
			}
			for(index=(*rowi)+1; index<=rowi_max; index++)	// it could enter any column from rowi+1 to rowi_max
			{		
				if ( (pos[1] >= (2*a*index+(index+1)*g)) && (pos[1] <= (2*a*(index+1)+(index+1)*g)) )	// enters column (index,colj+1)
				{				
					(*rowi) = index;
					(*colj) = (*colj)+1;
					particle_exit = 0;
				}
				else if ( (pos[1] <= (2*a*index+(index+1)*g)) && (pos[1] >= (2*a*((index-1)+1)+(index)*g)) )	
				// enters intercolumnar space, move this particle to the below column (index-1,colj+1) -> approximation assuming intercolumnar space is negligible
				{
					(*rowi) = index-1;
					(*colj) = (*colj)+1;
					particle_exit = 0;
				}
			}

			
		}
		else if ( pos[1] < (2*a*(*rowi)+((*rowi)+1)*g) )		// moves downward in -y direction
		{
			if (*rowi == 0)	// goes outside last row - taken as getting to the detector boundary and lost
			{
				num_lost++;
				particle_exit = 1;
				goto exitnow;
			}
			for(index=(*rowi)-1; index>=0; index--)	// it could enter any column from rowi+1 to rowi_max
			{		
				if ( (pos[1] >= (2*a*index+(index+1)*g)) && (pos[1] <= (2*a*(index+1)+(index+1)*g)) )	// enters column (index,colj+1)
				{				
					(*rowi) = index;
					(*colj) = (*colj)+1;
					particle_exit = 0;
				}
				else if ( (pos[1] >= (2*a*(index+1)+(index+1)*g)) && (pos[1] <= (2*a*(index+1)+(index+2)*g))  )	
				// enters intercolumnar space, move this particle to the above column (index+1,colj+1)
				{
					(*rowi) = index+1;
					(*colj) = (*colj)+1;
					particle_exit = 0;
				}
			}

			
		}
	}	// myplane=1 loop ends
	else if (myplane == -1)	// emits plane x = -a
	{
		if ( (pos[1] >= (2*a*(*rowi)+((*rowi)+1)*g)) && (pos[1] <= (2*a*((*rowi)+1)+((*rowi)+1)*g)) )	// enters adjacent column on left
		{
			(*colj) = (*colj)-1;	
			particle_exit = 0;
		}	
		else if ( pos[1] > (2*a*((*rowi)+1)+((*rowi)+1)*g) )	// moves upward +y direction
		{
			if (*rowi == rowi_max)	// goes outside last row - taken as getting to the detector boundary and lost
			{
				num_lost++;
				particle_exit = 1;
				goto exitnow;
			}
			for(index=(*rowi)+1; index<=rowi_max; index++)	// it could enter any column from rowi+1 to rowi_max
			{	
					
				if ( (pos[1] >= (2*a*index+(index+1)*g)) && (pos[1] <= (2*a*(index+1)+(index+1)*g)) )	// enters column (index,colj-1)
				{				
					(*rowi) = index;
					(*colj) = (*colj)-1;
					particle_exit = 0;
				}
				else if ( (pos[1] <= (2*a*index+(index+1)*g)) && (pos[1] >= (2*a*((index-1)+1)+(index)*g)) )	
				// enters intercolumnar space, move this particle to the below column (index-1,colj-1) -> approximation assuming intercolumnar space is negligible
				{
					(*rowi) = index-1;
					(*colj) = (*colj)-1;
					particle_exit = 0;
				}
			}

			
		}
		else if ( pos[1] < (2*a*(*rowi)+((*rowi)+1)*g) )		// moves downward in -y direction
		{
			if (*rowi == 0)	// goes outside last row - taken as getting to the detector boundary and lost
			{
				num_lost++;
				particle_exit = 1;
				goto exitnow;
			}
			for(index=(*rowi)-1; index>=0; index--)	// it could enter any column from 0 to rowi-1
			{		
				if ( (pos[1] >= (2*a*index+(index+1)*g)) && (pos[1] <= (2*a*(index+1)+(index+1)*g)) )	// enters column (index,colj-1)
				{				
					(*rowi) = index;
					(*colj) = (*colj)-1;
					particle_exit = 0;
				}
				else if ( (pos[1] >= (2*a*(index+1)+(index+1)*g)) && (pos[1] <= (2*a*(index+1)+(index+2)*g))  )	
				// enters intercolumnar space, move this particle to the above column (index+1,colj-1)
				{
					(*rowi) = index+1;
					(*colj) = (*colj)-1;
					particle_exit = 0;
				}
			}

			
		}
	}	// myplane=-1 loop ends
	else if (myplane == 2)	// emits plane y = a
	{

		if ( (pos[0] >= (2*a*(*colj)+((*colj)+1)*g)) && (pos[0] <= (2*a*((*colj)+1)+((*colj)+1)*g)) )	// enters adjacent column above

		{
			(*rowi) = (*rowi)+1;	
			particle_exit = 0;
		}	
		else if ( pos[0] > (2*a*((*colj)+1)+((*colj)+1)*g) )	// moves forward +x direction
		{
			if ((*colj) == colj_max) // goes outside last ending column - taken as getting to the detector boundary and lost
			{
				num_lost++;
				particle_exit = 1;
				goto exitnow;
			}

			for(index=(*colj)+1; index<=colj_max; index++)	// it could enter any column from colj+1 to colj_max

			{	
				if ( (pos[0] >= (2*a*index+(index+1)*g)) && (pos[0] <= (2*a*(index+1)+(index+1)*g)) )	// enters column (rowi+1,index)
				{				
					(*rowi) = (*rowi)+1;
					(*colj) = index;
					particle_exit = 0;
				}
				else if ( (pos[0] <= (2*a*index+(index+1)*g)) && (pos[0] >= (2*a*((index-1)+1)+(index)*g)) )	
				// enters intercolumnar space, move this particle to the column on its left (rowi+1, index-1) -> approximation assuming intercolumnar space is negligible
				{
					(*rowi) = (*rowi)+1;
					(*colj) = index-1;
					particle_exit = 0;
				}

			}


		}
		else if ( pos[0] < (2*a*(*colj)+((*colj)+1)*g) )		// moves backward in -x direction
		{
			if (*colj == 0)
			{
				num_lost++;
				particle_exit = 1;
				goto exitnow;
			}

			for(index=(*colj)-1; index>=0; index--)	// it could enter any column from 0 to colj-1
			{	
					
				if ( (pos[0] >= (2*a*index+(index+1)*g)) && (pos[0] <= (2*a*(index+1)+(index+1)*g)) )	// enters column (rowi+1,index)
				{				
					(*rowi) = (*rowi)+1;
					(*colj) = index;
					particle_exit = 0;
				}
				else if ( (pos[0] >= (2*a*(index+1)+(index+1)*g)) && (pos[0] <= (2*a*(index+1)+(index+2)*g))  )	
				// enters intercolumnar space, move this particle to the column on its right (rowi+1,index+1)
				{
					(*rowi) = (*rowi)+1;
					(*colj) = index+1;
					particle_exit = 0;
				}

			}


		}

	}	// myplane=2 loop ends
	else if (myplane == -2)	// emits plane y = -a
	{
		if ( (pos[0] >= (2*a*(*colj)+((*colj)+1)*g)) && (pos[0] <= (2*a*((*colj)+1)+((*colj)+1)*g)) )	// enters adjacent column below
		{
			(*rowi) = (*rowi)-1;	
			particle_exit = 0;
		}	
		else if ( pos[0] > (2*a*((*colj)+1)+((*colj)+1)*g) )	// moves forward +x direction
		{
			if ((*colj) == colj_max) // goes outside last ending column - taken as getting to the detector boundary and lost
			{
				num_lost++;
				particle_exit = 1;
				goto exitnow;
			}
			for(index=(*colj)+1; index<=colj_max; index++)	// it could enter any column from colj+1 to colj_max
			{	
				if ( (pos[0] >= (2*a*index+(index+1)*g)) && (pos[0] <= (2*a*(index+1)+(index+1)*g)) )	// enters column (rowi-1,index)
				{				
					(*rowi) = (*rowi)-1;
					(*colj) = index;
					particle_exit = 0;
				}
				else if ( (pos[0] <= (2*a*index+(index+1)*g)) && (pos[0] >= (2*a*((index-1)+1)+(index)*g)) )	
				// enters intercolumnar space, move this particle to the column on its left (rowi-1, index-1) -> approximation assuming intercolumnar space is negligible
				{
					(*rowi) = (*rowi)-1;
					(*colj) = index-1;
					particle_exit = 0;
				}

			}

			
		}
		else if ( pos[0] < (2*a*(*colj)+((*colj)+1)*g) )		// moves backward in -x direction
		{
			if ((*colj) == 0) // goes outside last ending column - taken as getting to the detector boundary and lost
			{
				num_lost++;
				particle_exit = 1;
				goto exitnow;
			}
			for(index=(*colj)-1; index>=0; index--)	// it could enter any column from 0 to colj-1
			{	
				if ( (pos[0] >= (2*a*index+(index+1)*g)) && (pos[0] <= (2*a*(index+1)+(index+1)*g)) )	// enters column (rowi-1,index)
				{				
					(*rowi) = (*rowi)-1;
					(*colj) = index;
					particle_exit = 0;
				}
				else if ( (pos[0] >= (2*a*(index+1)+(index+1)*g)) && (pos[0] <= (2*a*(index+1)+(index+2)*g))  )	
				// enters intercolumnar space, move this particle to the column on its right (rowi-1,index+1)
				{
					(*rowi) = (*rowi)-1;
					(*colj) = index+1;
					particle_exit = 0;
				}

			}


		}
	}	// myplane=-2 loop ends


exitnow:
 return particle_exit;	
}	// transmit function ends


// calculate directional cosines of reflected/refracted vector.
void next_dir_cos(float* dcos, float* normal, float refl_theta, float trans_theta, int flag_ref)
{
	float cos_angle = 0.0f;
	float norm = 0.0f;
	float dcos_temp[3] = {0.0f};
	int i=0;

	dcos_temp[0] = -dcos[0];
	dcos_temp[1] = -dcos[1];
	dcos_temp[2] = -dcos[2];
	
	cos_angle = dot_product(dcos_temp,normal);	// cosine of angle between incident in opposite direction and normal

	if (flag_ref == 0)		// reflection
	{
		for(i=0;i<3;i++)
		{
			dcos[i] = 2.0f*cos_angle*normal[i] + dcos[i];  // specular ray
		}
	}
	else if (flag_ref == 1)		// transmission	
	{
		 dcos[0]= -normal[0]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[0]+(cos_angle*normal[0]));
		 dcos[1]= -normal[1]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[1]+(cos_angle*normal[1]));
		 dcos[2]= -normal[2]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[2]+(cos_angle*normal[2]));
	}

	norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

	if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
	 {
		dcos[0] = dcos[0]/norm;
		dcos[1] = dcos[1]/norm;
		dcos[2] = dcos[2]/norm;
	 } 

return;	
}

// calculate rough normal depending on value of 'beta'.
void RoughSurface(float* normal, int* seed, float beta)
{

	float theta = 0.0f;
	float status = 0.0f;
	float rr = 0.0f;
	float normalpert1 = 0.0f;
	float normalpert2 = 0.0f;
	float normalpert3 = 0.0f;
	float rough_normal_x = 0.0f;
	float rough_normal_y = 0.0f;
	float rough_normal_z = 0.0f;
	float normalize_base = 0.0f;
	float normalize_pert=0.0f;


	status = ranecu(seed);
	normalpert3 = 2.0f*status - 1.0f;
	rr = sqrt(1.0f - status*status);
	status = ranecu(seed);
	theta = status * 2.0f * pi;
	normalpert1 = rr * cos(theta);
	normalpert2 = rr * sin(theta);
	
	normalize_pert = sqrt(pow(normalpert1,2) + pow(normalpert2,2 )+pow(normalpert3,2));

	normalpert1=normalpert1/normalize_pert;
	normalpert2=normalpert2/normalize_pert;
	normalpert3=normalpert3/normalize_pert;

	rough_normal_x = beta*normalpert1 + normal[0];	
	rough_normal_y = beta*normalpert2 + normal[1];
	rough_normal_z = beta*normalpert3 + normal[2];

	normalize_base = sqrt(pow(rough_normal_x,2) + pow(rough_normal_y,2 )+pow(rough_normal_z,2));

	normal[0] = rough_normal_x / normalize_base; 
	normal[1] = rough_normal_y / normalize_base;
	normal[2] = rough_normal_z / normalize_base;


return;
}

// determine if the photon gets detected at the sensor
int detection(float *pos, int num_rebound, float H, float lbound_x, float lbound_y, float ubound_x, float ubound_y, int pixelsize, struct start_info info)
{
        int result = 0;
	int pixelx = 0;
	int pixely = 0;
	int bin = 0;
	float distance = 0;
	float radius = 10.0f; // 10 um radius bins
	float area_diff = 0;

        //equation of sensor is z = -H/2
        // if a point satisfies above equation, it is detected	

        if (fabs(pos[2] - (float)(-H/2.0f)) < epsilon) 
	 {
		
                result = 1;
		counter++; 
		
		if( (pos[0] >= lbound_x) && (pos[1] >= lbound_y) && (pos[0] <= ubound_x) && (pos[1] <= ubound_y) )
		{
			pixelx = floor((pos[0]-lbound_x)/pixelsize);
			pixely = floor((pos[1]-lbound_y)/pixelsize);
			bin_matrix[pixelx][pixely]++;
		}
		
		distance = sqrt( ( (pos[0]-start_pos[0])*(pos[0]-start_pos[0]) ) + ( (pos[1]-start_pos[1])*(pos[1]-start_pos[1]) ) );
		bin = floor(distance/radius);
		
		if ((bin >= 0 ) && (bin <= 20)) // only consider first 21 bins
		{
			area_diff = pi*( ( ((bin + 1)*radius)*((bin + 1)*radius) ) - ( (bin*radius)*(bin*radius) ) );
			radialbin[bin] = radialbin[bin] + (1.0f/area_diff);	// radial bin obtained
		}

	 }
        else
            	result = 0;
  
 return result;
}


// find dot product of two vectors
float dot_product(float *aa, float *b)
{
        float result;

        result = aa[0]*b[0] + aa[1]*b[1] + aa[2]*b[2];

  return result;
}


// find minimum of three floats
float min(float a1, float a2)
{
	return (a1 < a2 ? a1 : a2);
}

float min3(float a1, float a2, float a3)
{
	return min(a1, min(a2,a3));
}


////////////////////////////////////////////////////////////////////////////////
//! Initialize the pseudo-random number generator (PRNG) RANECU to a position
//! far away from the previous history (leap frog technique).
//!
//! Each calculated seed initiates a consecutive and disjoint sequence of
//! pseudo-random numbers with length LEAP_DISTANCE, that can be used to
//! in a parallel simulation (Sequence Splitting parallelization method).
//! The basic equation behind the algorithm is:
//!    S(i+j) = (a**j * S(i)) MOD m = [(a**j MOD m)*S(i)] MOD m  ,
//! which is described in:
//!   P L'Ecuyer, Commun. ACM 31 (1988) p.742
//!
//! This function has been adapted from "seedsMLCG.f", see:
//!   A Badal and J Sempau, Computer Physics Communications 175 (2006) p. 440-450
//!
//!       @param[in] history   Particle bach number.
//!       @param[in] seed_input   Initial PRNG seed input (used to initiate both MLCGs in RANECU).
//!       @param[out] seed   Initial PRNG seeds for the present history.
//!
////////////////////////////////////////////////////////////////////////////////
// -- Upper limit of the number of random values sampled in a single track:
#define  LEAP_DISTANCE    1000
// -- Multipliers and moduli for the two MLCG in RANECU:
#define  a1_RANECU       40014
#define  m1_RANECU  2147483563
#define  a2_RANECU       40692
#define  m2_RANECU  2147483399

void init_PRNG(int history_batch, int histories_per_thread, int seed_input, int* seed)
{
  // -- Move the RANECU generator to a unique position for the current batch of histories:
  //    I have to use an "unsigned long long int" value to represent all the simulated histories in all previous batches
  //    The maximum unsigned long long int value is ~1.8e19: if history >1.8e16 and LEAP_DISTANCE==1000, 'leap' will overflow.
  // **** 1st MLCG:
  unsigned long long int leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
  int y = 1;
  int z = a1_RANECU;
  // -- Calculate the modulo power '(a^leap)MOD(m)' using a divide-and-conquer algorithm adapted to modulo arithmetic
  for(;;)
  {
      // printf(" leap, leap>>1, leap&1: %d, %d, %d\n",leap, leap>>1, leap&1);   

    // (A2) Halve n, and store the integer part and the residue
    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  
    {
      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2);  
      y = abMODm(m1_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
      if (0==leap) break;         // (A4) leap==0? ==> finish
    }
    else           // (leap is even)
    {
      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);   
    }
    z = abMODm(m1_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
  }
  // AjMODm1 = y;                 // Exponentiation finished:  AjMODm = expMOD = y = a^j

  // -- Compute and display the seeds S(i+j), from the present seed S(i), using the previously calculated value of (a^j)MOD(m):
  //         S(i+j) = [(a**j MOD m)*S(i)] MOD m
  //         S_i = abMODm(m,S_i,AjMODm)
  seed[0] = abMODm(m1_RANECU, seed_input, y);     // Using the input seed as the starting seed

  // **** 2nd MLCG (repeating the previous calculation for the 2nd MLCG parameters):
  leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
  y = 1;
  z = a2_RANECU;
  for(;;)
  {
    // (A2) Halve n, and store the integer part and the residue
    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number) 
    {
      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2);
      y = abMODm(m2_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
      if (0==leap) break;         // (A4) leap==0? ==> finish
    }
    else           // (leap is even)
    {
      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
    }
    z = abMODm(m2_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
  }
  // AjMODm2 = y;
  seed[1] = abMODm(m2_RANECU, seed_input, y);     // Using the input seed as the starting seed

}


/////////////////////////////////////////////////////////////////////
//!  Calculate "(a1*a2) MOD m" with 32-bit integers and avoiding   **
//!  the possible overflow, using the Russian Peasant approach     **
//!  modulo m and the approximate factoring method, as described   **
//!  in:  L'Ecuyer and Cote, ACM Trans. Math. Soft. 17 (1991)      **
//!                                                                **
//!  This function has been adapted from "seedsMLCG.f", see:       **
//!  Badal and Sempau, Computer Physics Communications 175 (2006)  **
//!                                                                **
//!    Input:          0 < a1 < m                                  **
//!                    0 < a2 < m                                  **
//!                                                                **
//!    Return value:  (a1*a2) MOD m                                **
//!                                                                **
/////////////////////////////////////////////////////////////////////
int abMODm(int m_par, int a_par, int s_par)
{
  // CAUTION: the input parameters are modified in the function but should not be returned to the calling function! (pass by value!)  
  int mval,aval,sval;
  mval=m_par; aval=a_par; sval=s_par;
  
  int qval, kval;
  int pval = -mval;            // p is always negative to avoid overflow when adding

  // ** Apply the Russian peasant method until "a =< 32768":
  while (aval>32768)        // We assume '32' bit integers (4 bytes): 2^(('32'-2)/2) = 32768
  {
    if (0!=(aval&1))        // Store 's' when 'a' is odd    
    {
      pval += sval;
      if (pval>0) pval -= mval;
    }
    aval >>= 1;             // Half a (move bits 1 position right)        
    sval = (sval-mval) + sval;       // float s (MOD m)
    if (sval<0) sval += mval;     // (s is always positive)
  }

  // ** Employ the approximate factoring method (a is small enough to avoid overflow):
  qval = (int) mval / aval;
  kval = (int) sval / qval;
  sval = aval*(sval-kval*qval)-kval*(mval-qval*aval);
  while (sval<0)
    sval += mval;

  // ** Compute the final result:
  pval += sval;
  if (pval<0) pval += mval;

  return pval;
}

////////////////////////////////////////////////////////////////////////////////
//! Pseudo-random number generator (PRNG) RANECU returning a float value
//! (single precision version).
//!
//!       @param[in,out] seed   PRNG seed (seed kept in the calling function and updated here).
//!       @return   PRN float value in the open interval (0,1)
//!
////////////////////////////////////////////////////////////////////////////////
float ranecu(int* seed)
{
  int i1 = (int)(seed[0]/53668);
  seed[0] = 40014*(seed[0]-i1*53668)-i1*12211;

  int i2 = (int)(seed[1]/52774);
  seed[1] = 40692*(seed[1]-i2*52774)-i2*3791;

  if (seed[0] < 0) seed[0] += 2147483563;
  if (seed[1] < 0) seed[1] += 2147483399;

  i2 = seed[0]-seed[1];
  if (i2 < 1) i2 += 2147483562;

  const float USCALE = 1.0/2147483563.0;       
  return ((float)(i2*USCALE));

}


