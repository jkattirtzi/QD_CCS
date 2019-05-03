#!	To compile the Fourier Transform program
#!	with gfortran type:
#!	gfortran FT2.f90 Fprog.f90
#!	The following information may help advise
#!	what is needed in the inham.dat file
#!	More information is given in the text file
#!	in the run directory. 
#!	
#!	#decides whether a window function is needed
#!	Winop Wdwf
#!	#gives the maximum energy for the FT
#!	FTEmax ETmax
#!	# decides the name of the Decay output file
#!	Decay Decay
#!	 decides the name of the FT output file
#!	FT FTname
#!	 decides the name of the ACF input file
#!	ACF ACFname
#!	decides the initial E value for FT
#!	FTE0 E0
#!	decides the max time for the FT before window function kills
#!	Tdec Tmax
#!	gives the number of lines need to be read in from the ACF file
#!	MAXJ jmax
#!
#!	The program should be run from the directory run/
