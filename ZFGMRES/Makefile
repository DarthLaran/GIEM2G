
ZGMRES = zPackgmres.o zPackfgmres.o  
ZGMRES_F = zPackgmres.f zPackfgmres.f 
all:$(ZGMRES_F)
	$(F77) $(FOPTS) -c $(ZGMRES_F)
	$(AR) rcs libzfgmres.a $(ZGMRES)
	rm $(ZGMRES)
clean :
	/bin/rm -f *.o *.mod
