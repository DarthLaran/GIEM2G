


UTILITIES=utilities.o
FGMRES_O=fgmres_interfaces.o gmres.o  seq_fgmres.o
ALL_O=$(UTILITIES) $(FGMRES_O)  


fgmres_lib: $(ALL_O)
	$(AR) rcs libgfgmres.a $(ALL_O)

%.o:%.F90
	$(FC) $(FOPTS) -c $*.F90 -o $*.o 
clean:
	rm $(ALL_O)   *.mod  
