############################################################################
# Makefile for oplist command line library sample
#
#   $Id: Makefile,v 1.1.1.2 2007/09/04 04:45:42 Exp $
#   $Log: Makefile,v $
#   Revision 1.1.1.2  2007/09/04 04:45:42 
#   Added FreeOptList.
#
#   Revision 1.1.1.1  2007/08/07 05:01:48 
#   Initial Release
#
#
############################################################################
CC = gcc
LD = gcc
CFLAGS = -I. -O3 -Wall -pedantic -ansi -c
LDFLAGS = -O3 -o

# libraries
LIBS = -L. -loptlist
EXE =
DEL = rm

all:	liboptlist.a

liboptlist.a:	optlist.o
	ar crv liboptlist.a optlist.o
	ranlib liboptlist.a

optlist.o:	optlist.c optlist.h
	$(CC) $(CFLAGS) $<

clean:
	$(DEL) *.o
	$(DEL) *.a
