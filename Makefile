# En Makefile som den burde være lavet fra starten
# Lavet til luna og venner.
# Lavet af Kim Brugger Sommeren 2000

PROG   = kmap
WARN   = #-Wshadow  -Wmissing-prototypes
FLAGS  = -pipe -Wall -ggdb -I$(ILIB) -D_POSIX_SOURCE $(WARN) -O6 -pg
CC     = gcc
LIBS   = -lm #-lefence #-ldmalloc  #-lefence #-lmcheck#-lmalloc
RM     = rm -f
MAKE   = gmake
OBJS   = kmap.o ktree.o

all :  $(PROG)  

kmap    :  $(OBJS) 
	$(CC) $(FLAGS)  -o  $@    $(OBJS) $(LIBS)

clean: 
	$(RM) $(PROG) *.o

veryclean : 
	$(RM)  *.o core $(PROG) *~ $(ILIB)*~ *.bak

depend:
	makedepend -- $(FLAGS) -- $(SRCS)
# DO NOT DELETE
