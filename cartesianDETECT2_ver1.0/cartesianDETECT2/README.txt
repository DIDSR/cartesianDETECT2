##########################################################################################################################################################
#
# ****Disclaimer****
#  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.
#
#	@file    README.txt
#       @author  Diksha Sharma (Diksha.Sharma@fda.hhs.gov)
#                Aldo Badano (Aldo.Badano@fda.hhs.gov)
#       @date    Oct 30, 2013
#
##########################################################################################################################################################

*************
INTRODUCTION
*************

cartesianDETECT2 is a dedicated Monte Carlo optical transport code for modeling pixelated scintillator structures.  It is based on DETECT2 (the optical transport code for MANTIS).


*******************************
CODE COMPILATION AND EXECUTION
*******************************

cartesianDETECT2 has been tested only on the Linux operating system. GNU gcc compiler and GNU scientific library needs to be pre installed before running cartesianDETECT2. A bash script is included to compile the code. A pre-compiled executable is also attached.

cartesianDETECT2 takes in 18 arguments which are as follows:

X, Y  detector dimensions (um)
Side and thickness of column (um)
Intercolumnar gap (um)
Refractive index of column and intercolumnar space
Top surface absorption
Bulk absorption (um^-1)
Surface roughness coefficient
Lower x, y PRF bounds (um)
Upper x, y PRF bounds (um)
Pixel pitch (um)
Seed input for RNG
Scintillation events file name (input)
PRF file name (output)
Radial response file name (output)


File 'sample_run_script.in' is a script containing a list of example run commands which will run in parallel on multiple processor nodes.
	
Sample inputs and outputs are included under the 'example' folder.


**********************************
cartesianDETECT2 PACKAGE CONTENTS
**********************************

1. cartesianDETECT2_optical_transport_main.c - main program in C
2. cartesianDETECT2_optical_transport_kernel.c - kernel in C
3. cartesianDETECT2_script.sh - compilation script
4. cartesianDETECT2.x - executable
5. sample_run_script.in - contains execute commands
6. /example/ - folder containing sample input and output files and figures
7. Manual_cartesianDETECT2.pdf - reference manual
8. README.txt - this file
