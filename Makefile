###############################################
# global Makefile with automatic dependencies #
###############################################


HOME_BOOST_54 := $(shell test -d $(HOME)/boost_1_54_0 && echo yes)
HOME_BOOST_49 := $(shell test -d $(HOME)/boost_1_49_0 && echo yes)
HOME_BOOST_47 := $(shell test -d $(HOME)/boost_1_47_0 && echo yes)
HOME_BOOST_46.1 := $(shell test -d $(HOME)/boost_1_46_1 && echo yes)
HOME_BOOST_46 := $(shell test -d $(HOME)/boost_1_46_0 && echo yes)
HOME_BOOST_45 := $(shell test -d $(HOME)/boost_1_45_0 && echo yes)
HOME_BOOST_43 := $(shell test -d $(HOME)/boost_1_43_0 && echo yes)
HOME_BOOST_41 := $(shell test -d $(HOME)/boost_1_41_0 && echo yes)

ifeq ($(HOME_BOOST_46),yes)
	BOOST_PATH = $(HOME)/boost_1_46_0
	BOOST_VERSION = 46
else ifeq ($(HOME_BOOST_54),yes)
	BOOST_PATH = $(HOME)/boost_1_54_0
	BOOST_VERSION = 54
else ifeq ($(HOME_BOOST_46.1),yes)
	BOOST_PATH = $(HOME)/boost_1_46_1
	BOOST_VERSION = 46.1
else ifeq ($(HOME_BOOST_47),yes)
	BOOST_PATH = $(HOME)/boost_1_47
	BOOST_VERSION = 45
else ifeq ($(HOME_BOOST_45),yes)
	BOOST_PATH = $(HOME)/boost_1_45_0
	BOOST_VERSION = 43
else
	BOOST_PATH =
endif

ifeq ($(BOOST_PATH),)
BOOST_H = /usr/include
else
BOOST_H = $(BOOST_PATH)
endif

ifndef OSTYPE
  OSTYPE = $(shell uname -s|awk '{print tolower($$0)}')-$(shell uname -m)
endif

ifeq ($(OSTYPE),darwin-i386)
	LDIR=/usr/local/lib/
else
	LDIR=/usr/lib/
endif


GCC_VERSION := $(shell gcc -dumpversion)
GCC_SHORT_VERSION=43

##################
# compiler flags #
##################

INCLUDE_DIRS = -I. -I$(BOOST_H)

THREADING = # -lboost_thread -L$(BOOST_PATH)/stage/lib

# FLAGS        = -O3 -Wall -g -march=$(CPU) -pipe $(INCLUDE_DIRS) -DWno_deprecated # -O2 # -mfpma# th=sse -msse2 -O2
# FLAGS        = -Wall -g -O2 -march=$(CPU) -pipe $(INCLUDE_DIRS)# -O2 # -mfpmath=sse -msse2 -O2
FLAGS        = -Wall -g  $(INCLUDE_DIRS)  -O2 -DNDEBUG  -fpermissive
# FLAGS        = -Wall -O0 -ggdb -march=$(CPU) -pipe $(INCLUDE_DIRS) # -O2 # -mfpmath=sse -msse2 -O2

CXXFLAGS     = $(FLAGS)


################
# object files #
################

DIRS     := .
OBJ := 	main.o \
        BLAS.o \
	Util.o \
	Automaton.o \
	WordTree.o \
	RandomValueMap.o \
	Timer.o \
	mpreal.o \
	dlmalloc.o

####################
# build everything #
####################
all : 	equiv



#####################################
# static library: PASS object files #
#####################################

libequiv.a : $(OBJ)
	@ar -rus libequiv.a $(OBJ)

####################
# build executable #
####################
equiv:	libequiv.a
	echo $(OSTYPE)
	$(CXX) $(FLAGS) $(THREADING) -o equiv -lequiv -L. -lmpfr -lgmp
	@echo Build is complete.


fusion:
	$(CXX) $(FLAGS) -o fusion
	@echo Fusion built.

###########
# cleanup #
###########
clean:
	# @if [ -e Makefile.dep ]; then rm Makefile.dep; fi
	@for obj in $(OBJ); do \
		(if [ -e $$obj ]; then rm $$obj; fi) \
		done
	rm libequiv.a

######################
# build dependencies #
######################
Makefile.dep:
	@if [ ! -e Makefile.dep ]; then echo "# automatic dependencies" > Makefile.dep; fi
	@makedepend -w -a -Y -fMakefile.dep -- $(FLAGS) $(INCLUDE_DIRS) -- main.cpp *.hpp &> /dev/null
	@for dir in $(DIRS); do \
		(makedepend -w -a -Y -fMakefile.dep -- $(FLAGS) $(INCLUDE_DIRS) -- $$dir/*.hpp $$dir/*.h $$dir/*.cpp) \
	done &> /dev/null

depend:
	@if [ -e Makefile.dep ]; then rm Makefile.dep; fi
	make Makefile.dep

########################################
# automatically generated dependencies #
########################################
ifeq "$(wildcard Makefile.dep)" ""
else
include Makefile.dep
endif
