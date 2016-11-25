#####################################################################
#
# RepeatScout Project Build Script
#
#####################################################################
# $Log: Makefile,v $
# Revision 1.4  2008/08/07 21:59:12  rhubley
#   - An error reported by Eric Ganko in the filter-stage-1.prl script.
#
#
####################################################################

# Set the version here
VERSION = 1.0.5

# Installation Directory
INSTDIR = /usr/local/RepeatScout-$(VERSION)

CFLAGS = -O3 -Wall
LIBS = -lm
OBJ= cmd_line_opts.o version.o

HDR= cmd_line_opts.h 

all: RepeatScout build_lmer_table 

RepeatScout: build_repeat_families.o build_repeat_families.h $(OBJ) $(HDR)
	$(CC) build_repeat_families.o $(OBJ) -o $@ $(LIBS)

build_lmer_table: build_lmer_table.o build_lmer_table.h $(OBJ) $(HDR)
	$(CC) build_lmer_table.o $(OBJ) -o $@ $(LIBS)

version.c: Makefile
	echo "char const* Version = \"$(VERSION)\";" > version.c

clean:
	@rm *.o build_lmer_table RepeatScout

.c.o:
	$(CC) $(CFLAGS) -c $< -o $*.o $(CCINCLUDES)

install: all
	@mkdir $(INSTDIR)
	cp RepeatScout $(INSTDIR)
	cp README $(INSTDIR)
	cp build_lmer_table $(INSTDIR)
	cp filter-stage-1.prl $(INSTDIR)
	cp filter-stage-2.prl $(INSTDIR)
	cp merge-lmer-tables.prl $(INSTDIR)
	cp compare-out-to-gff.prl $(INSTDIR)

distribution:
	rm *~
	(cd ../; tar zcvf RepeatScout-$(VERSION).tar.gz RepeatScout-1 --exclude RepeatScout-1/orig  --exclude RepeatScout-1/tests  --exclude RepeatScout-1/CVS  --exclude RepeatScout-1/rc-change-w-debug)
