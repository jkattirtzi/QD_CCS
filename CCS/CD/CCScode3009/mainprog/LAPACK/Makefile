# Include platform dependent options.
include make.options

lapack:	dep.f driv.f 
	$(F90) $(F90FLAGS) -c dep.f driv.f 
	ar rcs liblapackJK.a dep.o driv.o 
#	rm -f dep.o driv.o 

lcleana:	
	rm liblapackJK.a
lcleano:
	rm -f dep.o driv.o          
