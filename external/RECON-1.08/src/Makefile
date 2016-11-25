# Makefile for ReCon
#
#
VERSION = 1.08

## where you want things installed
BINDIR  = ../bin
MANDIR  = $(HOME)/man

## sourcefiles
#SRCDIR = .

## your compiler
#CC = cc    # for SGI Origin200 compiler#
CC = gcc

## any special compiler flags you want
#  -pedantic clashes with -DMEMDEBUG??
#
#CFLAGS = -g -Wall -pedantic
CFLAGS = -O
#CFLAGS = -g
#CFLAGS = -TARG:exc_max -64 -Ofast=ip27   #specific SGI Origin200  compiler#

## machine specific definitions
#  None should be needed
#
# Development:
# -DASSERT        Diagnostics level 1. Assertions. 
# -DVERBOSE       Diagnostics level 2. Verbosity.  (unused at the moment)
# -DMEMDEBUG      Diagnostics level 3. Memory checking. (maybe slow!)
#                                      (must link w/ -ldbmalloc) 
# -DLINTING       Suppress overzealous warnings from gcc -Wall and lint   
#
#MDEFS = -DLINTING -DASSERT -DMEMDEBUG
#MDEFS = -DASSERT 
 
## how to install the man pages 
## cp is fine
INSTMAN   = cp
MANSUFFIX = 1

## Your compression program
COMPRESS = gzip               # GNU compression (better, smaller)
#COMPRESS = compress          # generic UNIX compress (portable)

#######
## should not need to modify below this line
#######
SHELL  = /bin/sh
LIBS   = -lsquid -lm
#LIBS = -ldbmalloc -lsquid -lm 

READMES = README 
#INSTALL COPYING DEMO GNULICENSE FILES

PROGS = imagespread eledef eleredef edgeredef famdef

MANSRC= 

DOCS =  

HDRS =  bolts.h seqlist.h msps.h ele.h

MAIN =  

SRC =   imagespread.c eledef.c eleredef.c edgeredef.c famdef.c


OBJ =   

DISTFILES = $(MAIN) $(SRC) $(HDRS) $(MANSRC) $(READMES) $(DOCS)

TESTSUITE = 

DEMOS = $(HOME)/Demos

all: 	$(PROGS)


install: $(PROGS)
	cp $(PROGS) $(BINDIR)/

clean:
	-rm -f imagespread eledef eleredef edgeredef famdef *.o *~ *.Addrs *.Counts *.pixie Makefile.bak core $(PROGS) TAGS

tags:
	etags -t $(DISTFILES)

dist:
	@if test -d recon-$(VERSION);        then rm -rf recon-$(VERSION);    fi
	@if test -f recon-$(VERSION).tar.gz; then rm recon-$(VERSION).tar.gz; fi
	@if test -f recon-$(VERSION).tar.Z;  then rm recon-$(VERSION).tar.Z;  fi
	mkdir recon-$(VERSION)/
	cp $(DISTFILES) recon-$(VERSION)/
	cp Makefile.dist recon-$(VERSION)/Makefile
	mkdir recon-$(VERSION)/Testsuite
	for testfile in $(TESTSUITE); do \
	  cp Testsuite/$$testfile recon-$(VERSION)/Testsuite/;\
	done
	mkdir recon-$(VERSION)/Demos
	for demofile in $(DEMOS); do\
	  cp Demos/$$demofile recon-$(VERSION)/Demos/;\
	done
	tar cvf recon-$(VERSION).tar recon-$(VERSION)
	$(COMPRESS) recon-$(VERSION).tar

test:   
	(cd Testsuite; testsuite.sh)

lint:
	lint -h -I$(SQUIDDIR) $(RFLAGS) $(MAIN) $(SRC) $(SQUIDDIR)/llib-lsquid.ln /usr/lib/lint/llib-lm.ln

.c.o:
	$(CC) $(CFLAGS) $(MDEFS) -I$(SQUIDDIR) -c $<		

.m.o:
	$(CC) $(CFLAGS) $(MDEFS) -I$(SQUIDDIR) -c $<





