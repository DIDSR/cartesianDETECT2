##########################################################################################################################################################
#
# ****Disclaimer****
#  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.
#
#	@file    README_example.txt
#       @author  Diksha Sharma (Diksha.Sharma@fda.hhs.gov)
#                Aldo Badano (Aldo.Badano@fda.hhs.gov)
#       @date    Oct 29, 2013
#
##########################################################################################################################################################

cartesianDETECT2 requires an input file giving the information regarding the location where the energy is deposited and further optical photons created, and generates two output files.  The input file lists x, y, z, energy in the first line and a pound sign in the next line to inform the program that the file has reached its end.  Every energy deposition event is stored as a separate input file so that parallel runs can be started on multiple processor cores.

In this example, we have included sample input files at the top, middle and bottom of the scintillator.  All the values are in microns.  

'-490' means at -490 microns and is closest to the sensor plane (bottom), '490' farthest (top) and '0' in the middle.
If the detector thickness if 'H' microns, then the z location goes from -H/2 to +H/2 microns, 0 microns always being in the center of the scintillator.

The code generates two output files - radialbin and prf.  Apart from this one may choose to save the information printed on the console during the simulation run to a file.  We have named this as 'console'.

radialbin*.out - the radial response obtained by aggregating the detected optical photons at the sensor plane in the bin size of 10 microns.  The radialbin images included under the example folder shows one solid blue curve and another dashed red curve.  The solid blue curve is the actual radial response obtained from the prf (details in the Manual_cartesianDETECT2.pdf) and the dashed red curve is the Gaussian fit on the blue curve.

prf*.out - point response function image data from 4500 to 5500 microns in x and y directions with pixel pitch of 2 microns (501 x 501 pixels).  From the PRF images included in the folder, we can see that as we move closer to the sensor plane (-490 microns) teh sharper the PRF peak becomes, as expected.

console*.out - lists the input arguments provided, optical transport statistics, total time and speed of the simulation at that particular depth.
