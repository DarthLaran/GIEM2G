
ZGMRES = zPackgmres.o zPackfgmres.o  
ZGMRES_F = zPackgmres.f zPackfgmres.f 
all:$(GMRES_F)
	$(F77) $(FOPTS) -c $(ZGMRES_F)
	ar rcs libzfgmres.a $(ZGMRES)
	rm $(ZGMRES)
clean :
	/bin/rm -f *.o *.mod
