#
#************************************************************
#	ReadMe.txt 		author: John Kattirtzi
#		
#	Welcome to the Coupled Coherent States (CCS) CD
#	This file gives an overview of what is  included i
#	on the CD. 
#
#	The CCS software is split into three programs,
#	these are found in the directory CCScode3009.
#	The software is found in the three directories:
#	1. basis/
#	2. FT/
#	3. mainprog/
#       4. run/
#	Comments can be written (#) 
#	in the input.dat and
#	inham.dat files which are ignored when reading these 
#	files
#************************************************************
#	To set up a basis set for the simulation:
#	go into the directory basis/
#	the source code is found in the directory
#	srcbs. To compile this program the compile script
#	can be run. This will generate the executable file
#	Mkprog in this directory and in the run directory. 
#	The Mkprog is used to generate the basis set. 
#	Instructions on how to run the program 
#	can be found in the run directory 
#       details on the individual modules and the
#	subroutines can be found in the source files.
#***********************************************************
#	To run	a simulation
#	The basisset.dat file and the inham.dat file are 
#	required. The same inham.dat as for the basis set 
#	program should be used. 
#	Run the executable CCSProg
#	Example inham.dat and basisset.dat files are given
#	in the run directory
#	Further notes on explaining the inham.dat file are
#	also provided there
#***********************************************************
#	To calculate the Fourier Transform:
#	Run the FT.out program in the FT directory.
#	The input file that is required is inham.dat
#	There is a FTparameters section in the inham.dat
#	file, provided, which must be used to calculate
#	the Fourier Transform. 
#	Additional notes explaining how to use this program
#	are given in the text file in the run directory. 
#***********************************************************
#	The run directory is recommended for running 
#	the CCS simulation and for calculating the 
#	Fourier Transform
#***********************************************************
#	Enjoy using the software!!!!!!!!!!!!!!!!!!!!!!!!!!
#	If you have any difficulties contact:
#	cmjak@leeds.ac.uk or jkattirtzi@yahoo.co.uk
#***********************************************************
