# Master Makefile for Sequence Utilities for Genome analysis
# Type make to build

O = o
SRC_DIR = ./



CFLAGS =  -I./              		                     # compiler switches to be applied to every module
OPTIM_SPEED = -O3             	  	                     # switches that give speed priority over size
OPTIM_SIZE = -O1              	  	                     # switches that give size priority over speed
CC = g++
OPTIONS = $(CFLAGS) $(INCLUDE)


all: map_orthology 



MAP_ORTHOLOGY_OBJS = map_orthology.$(O)  \
			gen_dna_funcs.$(O) 


map_orthology: $(MAP_ORTHOLOGY_OBJS) 
	$(CC) $(LINUX_BUILD)  $(LIBRARY_DIR) \
	 -o map_orthology $(OPTIONS) $(MAP_ORTHOLOGY_OBJS) 


%.o: %.cpp
	$(CC)  $(OPTIONS) $(OPTIM_SPEED) -c $<


%.o: %.c
	$(CC)  $(OPTIONS) $(OPTIM_SPEED) -c $<





