ZGMRES = zPackgmres.o zPackfgmres.o  
ZGMRES_F = zPackgmres.f zPackfgmres.f
PATCH=gmres_patch fgmres_patch
all:  $(PATCH)
	cp $(FGMRES_PATH)/zPackgmres.f ./
	cp $(FGMRES_PATH)/zPackfgmres.f ./
	patch zPackgmres.f <gmres_patch
	patch zPackfgmres.f <fgmres_patch
	$(F77) $(FOPTS) -c zPackgmres.f 
	$(F77) $(FOPTS) -c zPackfgmres.f 
	$(AR) rcs libzfgmres.a $(ZGMRES)
	rm $(ZGMRES)
clean :
	/bin/rm -f *.o *.mod

