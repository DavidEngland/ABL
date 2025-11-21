# program name
PROGRAM       =   radiationmodel
DEST          =   /n1/data/linux/bin 

# source files
         
SRCS          = RadParams_flcode.f aerosol_init_flcode.f chou_routines_flcode.f trmm_wnflt_0599_flcode.f flrad0100_flcode.f oned_rad_driv_flcode.f iaform3_flcode.f misc_subs_flcode.f tau_wp_oct96_flcode.f water_hu2_flcode.f dkfrwn_0599_flcode.f aerosol_data_flcode.f model_Rad.f 
          
# object files
 
OBJS          =  RadParams_flcode.o aerosol_init_flcode.o chou_routines_flcode.o trmm_wnflt_0599_flcode.o flrad0100_flcode.o oned_rad_driv_flcode.o iaform3_flcode.o misc_subs_flcode.o tau_wp_oct96_flcode.o water_hu2_flcode.o dkfrwn_0599_flcode.o aerosol_data_flcode.o model_Rad.o

# libraries and include directories
LIBS          = -lm 

USRLIBS         =  

INCLUDES        =  

# C, Fortran, load and flags
CFLAGS	      = -O   -w $(INCLUDES)
FFLAGS	      =  -Mvect=cachesize:524288 -Munroll -Mnoframe -O2 -pc 64 -Mfree -Mx,119,0x200000 -byteswapio
LDFLAGS	      = -O 

# other definitions (must run under sh shell)
CC 	      = cc
F77	      = pgf90
F90	      = pgf90
LD	      = $(F90)
INSTALL	      = /etc/install
MAKEFILE      = Makefile
PRINT	      = lp
SHELL	      = /bin/sh

.f.o:
	$(F90) -c $(FFLAGS) -o $@ $<

all:		$(PROGRAM)

# create the exectuable
$(PROGRAM):     $(OBJS) $(USRLIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(USRLIBS) -o $(PROGRAM) $(LIBS)
		@echo "done"

# remove object code and core files but leave executables
clean:;		@rm -f $(OBJS) core

# remove object code, executables and core
clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

# update dependencies of objects and headers in makefile
depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

# install executable in destination directory
install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(mv $(PROGRAM) $(DEST)); fi

# print source code and headers
print:;		@$(PRINT) $(HDRS) $(SRCS)

# create C tags table
tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

# update program
update:		$(DEST)/$(PROGRAM)
###

