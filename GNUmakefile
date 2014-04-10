#+
#+  $Id: GNUmakefile, 2014-02-20 14:36:15 liudongzy $
#+  Author(s):
#+    Liu Dong (dliu13@mail.ustc.edu.cn) 20/02/2014
#+


name := mirror
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CPPFLAGS += $(ROOTCFLAGS)
EXTRALIBS += $(ROOTLIBS) $(ROOTGLIBS)
#CPPFLAGS += -O -Wall -fPIC -I$(ROOTSYS)/include -g
#OOTLIBS      = -L$(ROOTSYS)/lib -lNew -lBase -lCint -lClib -lCont -lFunc -lGraf -lGraf3d -lHist -lHtml -lMatrix -lMeta -lMinuit -lNet -lPostscript -lProof -lTree -lUnix -lZip
#ROOTGLIBS     = -lGpad -lGui -lGX11 -lX3d
#EXTRALIBS += $(ROOTLIBS) $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11 -lm -ldl

#EXTRALIBS += -L$(ROOTSYS)/lib


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
