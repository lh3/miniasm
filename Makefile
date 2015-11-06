CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
CPPFLAGS=
INCLUDES=	-I.
OBJS=		sys.o sdict.o paf.o asg.o common.o hit.o asm.o
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

asg.o: asg.h kvec.h ksort.h
asm.o: miniasm.h sdict.h asg.h kvec.h kseq.h
common.o: miniasm.h sdict.h asg.h
hit.o: sdict.h paf.h kvec.h sys.h miniasm.h asg.h ksort.h
main.o: kvec.h sys.h paf.h sdict.h miniasm.h asg.h
paf.o: paf.h kseq.h
sdict.o: sdict.h khash.h
