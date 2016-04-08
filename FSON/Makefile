FSON= fson_string_m.o fson_value_m.o fson_path_m.o   fson.o 

libfson: $(FSON)
	$(AR) rcs libfson.a $(FSON)
%.o:%.f90
	$(FC) $(OPTS) -c $*.f90 -o $*.o 
clean :
	/bin/rm -f *.o *.mod

