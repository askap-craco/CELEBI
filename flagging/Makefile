DESTDIR = $(PWD)/
CC=gcc
CLIBS= -lgsl -lgslcblas -lm
CFLAGS=-I. -I./include -g -fopenmp
c = $(CC) $(CFLAGS) 
o = -o$(DESTDIR)$@

# default ways to compile
.c.o: ; $c -c  $< 
.f.o: ; $f -c $< 
.c: ; $c $o $< $(CLIBS)
.f: ; $f $o $< 

make: ankflag

ANKOBJ=ankflg.c exc_fns.c flg_fns.c basic_fns.c aux_fns.c init_pars.c
ankflag:$(ANKOBJ)
	$c $(ANKOBJ) $o $(CLIBS)

clean:
	rm ankflag
