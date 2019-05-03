# This file is here to help explain the input file for inham.dat 
# so it can be used for the main program and the fourier transform program

#!*************************************************************************
#!
#!	Main Program		Author: John Kattirtzi
#!
#!	To compile the main program use the compile script in the mainprog 
#!	directory
#! 	The clean script may help clean some files
#!	The compile script will produce the CCSProgO3 program
#!	This executable can be run using the inham.dat and 
#! 	the basisset.dat file.
#!	The inham.dat file must:
#!	specify the system - as shown in the example inham.dat provided
#!	give the initial time and max time
#!	give the names of the output files that are required e.g
#!	outp name.outp 
#!	dn is an integer that determines how often the output is written
#!	if the name of the output is not specified then the file will not
#!	be produced
#!	outp- gives general information about the system
#!	ACF- gives autocorrelation function
#!	Pop- gives the populations for the Hubbard model
#!	Traj- gives the trajectories of the basis functions
#!	Consv -gives the conservation of the norm and the classical energy	
#!	The order in which these are placed does not affect the program
#!*************************************************************************	
#	FT.out - Fourier Transform Program
#
#	This program multiplies the ACF by the winodw function-if required
#	and calculates the Autocorrelation Function
#	It is found in the FT directory
#	The ACF file that is read in is specified in the inham.dat file
#	by ACF name.acf
#	For the window function to be used then
#	'Winop YES' 
#	must be specified
#	Decout option gives the output for the ACF*window function
#	FTout  option gives the output of the Fourier Transform
#	The max value for E for the FT(E) is given by FTEmax
#	The initial value of E is given by E0
#	The time step used in the calculation is given by dt (earlier)
#	The FTcut gives the t_cut for the window function, beyond this point
#	the ACF is multiplied by 0
#	MAXJ is the number of lines that need to be read in from the .acf file
#	This is given at the end of the ACF file from the main simulation.  
#!*************************************************************************	
