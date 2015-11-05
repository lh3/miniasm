CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
CPPFLAGS=
INCLUDES=	-I.
OBJS=		sys.o sdict.o paf.o
PROG=		miniasm
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

miniasm:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE
