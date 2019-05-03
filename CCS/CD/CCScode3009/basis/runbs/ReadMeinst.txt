#
#************************************************************
#	ReadMeinst.txt		Author:John Kattirtzi
#	This files provides information on how to set up a 
#	basis set.
#	The program that configures the basis set is 
#	called Mkprog and this uses the input.dat 
#	and inham.dat input files.
#	The input.dat file has information regarding the
#	basis set whilst the inham.dat file has information
#	regarding the Hamiltonian of the system. This
#	is required when energy/population restrictions on the 
#	basis functions are required. 
#	The structure of the input files discussed here and in
#	the examples provided are a guide and can be changed. 
#*************************************************************
#	The input.dat file:
#	ndim is an integer for the number of dimensions
#	hbar affects the units used, the standard is 1.0d0
#	nbf is an integer for the number of basis functions
#	SYM is the option for a symmetrical basis set i.e one
#	that follows the symmetrical arguments:
#	Z1(0)=Z2*(0)
#	B1(0)=B2*(0)
#	For a wave function |Psi>=B|Z>
#	gamma,mass and wfreq (omega) are required for
#	the basis functions where gamma=mass*wfreq/hbar
#	these should be set to 1.0d0 but they can be changed
#	SEED is an integer that is used to generate numbers 
#	from a Gaussian function. Setting this to 0 generates
#	a SEED value. To run the program many times and 
#	produce a different SEED value each time requires a 5 
#	second wait.
#	outp gives the file unit that seed is written to (130) is pref.
#	mean p gives p_0 for the Gaussian distribution in each dim
#	mean q gives q_0 for the Gaussian distribution in each dim
#	sigp gives the standard deviation in the Gaussian for selecting p
#	sigq gives the standard deviation in the Gaussian for selecting q
#
#	Syntax when a paramter for each degree of freedom is read in:
#	The gamma for each degree of freedom
#	syntax is :
#	word(gamma) (degree number) gamma(value)
#	if you have more than one degree of freedom
#	put it directly below the first
#
#	Further details on the syntax and structure of the input file can be 
#	found in the input.dat file. 
#*****************************************************************************
#	inham.dat file
#	This file gives information on the system that is being used.
#	The systems programmed so far are the Free Particle, 
#	Harmonic Oscillator, Morse Oscillator and Bose-Hubbard Model. 
#	Examples of how to input these is shown in the comments in
#	the inham.dat file. 
#	An Ebfmin and and Ebfmax value is required, which give the min and max
#	energy of the basis function. 
#	If population constraint are required then POPChk should be YES
#	and Popmin and Popmax must be specified
#	Ntries gives the number of tries that will be used to find
#	a basis function that meets these constrains before giving up
#
#*****************************************************************************
#	output basisset.dat file
#	This provides the basis set file, which can be used multiple times
#	for different simulations. 
#	ndof gives the number of degrees of freedom for the basis function
#	nbasisfns gives the number of basis functions
#	basis N, where N is the integer is the basis function Number
#	C is the coefficient for the basis function
#	the complex z values are given below this
#	with each degree of freedom given on the line below
#	For a symmetrical basis set with N basis functions 
#	Basis x= conjg(basis(N+x)) where x is an integer
#	i.e for 10 basis functions 
#	basis 1 = conjugate of basis 6
#******************************************************************************
#	output General.out
#	This files gives general information such as the SEED (if outp =130)
#	If tight constraints are used it will indicate if choosing the
#	basis functions	with these constraints was possible for ntries. 
