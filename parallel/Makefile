CC = mpicc
CFLAGS = -O2
LDFLAGS = -L../simple-jpeg -lsimplejpeg

PROJ = parallel_main
OBJS = parallel_main.o

all : simple-jpeg $(PROJ)

parallel_main : $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

simple-jpeg :
	$(MAKE) -C ../simple-jpeg

obj-clean:
	$(RM) *.o *.out

exec-clean:
	$(RM) $(PROJ)

autosave-clean:
	$(RM) *~

clean:
	$(MAKE) obj-clean
	$(MAKE) exec-clean
	$(MAKE) autosave-clean
