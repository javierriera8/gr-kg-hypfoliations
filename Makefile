#! /usr/sbin/smake
#

NAME	= potential4_lam100_mu200                                                                                                                              # name of the executable that is produced
NAMP	= main

OBJD    =./obj
EXED    =./exe

FC	= gfortran
F90	= gfortran
CC	= gcc
LD	= $(F90)
MOD	= mod
SWP	=
RM	= /bin/rm -f
MP	=
ABI	=
ISA	=
ARCH	= $(MP) $(ABI) $(ISA)
OLEVEL	= -O3
FOPTS	=
PROFI   = #-pg
STATIC	= #-static -static-libgfortran
F90OPTS = -fimplicit-none -ffree-line-length-none -ffree-form -J ./obj $(OBJD) $(PROFI) #-ggdb -r16 # -Wall -pedantic # these two show warnings
COPTS	=
F90FLAGS= $(ARCH) $(OLEVEL) $(F90OPTS)
FFLAGS	= $(ARCH) $(OLEVEL) $(FOPTS)
CFLAGS	= $(ARCH) $(OLEVEL) $(COPTS)
LIBPATH =
LIBS	=
LDFLAGS	= $(STATIC) $(ARCH) $(OLEVEL) $(PROFI)
LIBFLAGS = $(LIBPATH) $(LIBS)
PROF	=

VPATH = $(OBJD):$(EXED):molsrc

PROG =	./exe/$(NAME)

SRCS =	precision.f90 numservice.f90 command_line.f90 grid.f90 mol_subroutines.f90 mol.f90 molparams.f90 output.f90 \
	parameters.f90 $(NAMP).f90 \
	userparams.f90

OBJS =	precision.o  numservice.o molparams.o userparams.o \
	parameters.o grid.o mol_subroutines.o mol.o output.o \
	command_line.o $(NAMP).o

OBJS2 = $(OBJD)/precision.o \
	$(OBJD)/numservice.o $(OBJD)/molparams.o $(OBJD)/userparams.o \
        $(OBJD)/parameters.o $(OBJD)/grid.o $(OBJD)/mol_subroutines.o $(OBJD)/mol.o $(OBJD)/output.o \
        $(OBJD)/command_line.o $(OBJD)/$(NAMP).o

CNTRL = *.in makemake Make* Scripts

all: $(OBJD) $(EXED) $(PROG)

$(OBJD):
	mkdir $(OBJD)

$(EXED):
	mkdir $(EXED)


$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS2) $(LIBFLAGS)

clean:
	rm -rf  $(OBJD)  $(OBJS2) *.$(MOD) last_data.txt $(EXED)/$(NAME)

cleandat:
	rm *.dat last_data.txt

tar:
	tar cf `basename $(PWD)`.tar $(SRCS) *.in make* Make*
	gzip -f `basename $(PWD)`.tar

zip:
	zip -r `basename $(PWD)`.zip $(SRCS) *.in Make*

gedit:
	gedit main.f90 source.f90 molsrc/mol_subroutines.f90 molsrc/mol.f90 molsrc/output.f90 molsrc/numservice.f90 run/parameter_template.txt Makefile &

wrang:
	open -a /Applications/TextWrangler.app/Contents/MacOS/TextWrangler main.f90 source.f90 molsrc/mol_subroutines.f90 molsrc/mol.f90 molsrc/output.f90 molsrc/numservice.f90 run/parameter_template.txt Makefile &

atom:
	atom main.f90 source.f90 molsrc/mol_subroutines.f90 molsrc/mol.f90 molsrc/output.f90 molsrc/numservice.f90 run/parameter_template.txt Makefile &

.SUFFIXES:
.SUFFIXES: $(SUFFIXES) .f90 .f .o

.f90.o:
	$(F90) $(FREE) $(F90FLAGS) -o $(OBJD)/$@ -c $<

.f.o:
	$(FC) $(FIXED) $(FFLAGS) -c $<
