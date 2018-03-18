CC = gcc
CFLAGS = -O2
LDFLAGS = -L../simple-jpeg/ -lsimplejpeg

PROJ = serial_main
OBJS = serial_main.o

all : simple-jpeg $(PROJ)

serial_main : $(OBJS)
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
