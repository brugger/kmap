# En Makefile som den burde v�re lavet fra starten
# Lavet til luna og venner.
# Lavet af Kim Brugger Sommeren 2000

PROG   = kmap
WARN   = #-Wshadow  -Wmissing-prototypes
FLAGS  = -pipe -Wall -ggdb -I$(ILIB) -D_POSIX_SOURCE $(WARN) -O6 #-pg
CC     = gcc
LIBS   = -lm #-lefence #-ldmalloc  #-lefence #-lmcheck#-lmalloc
RM     = rm -f
MAKE   = gmake


all :  $(PROG)  

kmap    :  kmap.o $(OBJS) 
	$(CC) $(FLAGS)  -o  $@  $<  $(OBJS) $(LIBS)

clean: 
	$(RM) $(PROG)

veryclean : cleandoc
	$(RM)  *.o core $(PROG) *~ $(ILIB)*~ *.bak

depend:
	makedepend -- $(FLAGS) -- $(SRCS)
# DO NOT DELETE