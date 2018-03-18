TAR=tar

all : serial parallel

.PHONY : serial parallel clean

serial :
	$(MAKE) -C $@

parallel :
	$(MAKE) -C $@

delivery : clean
	$(TAR) cvzf inf3380_maalgaba_oblig1.tar.gz *

clean :
	$(MAKE) -C parallel clean
	$(MAKE) -C serial clean
	$(MAKE) -C simple-jpeg clean
	$(RM) *~
