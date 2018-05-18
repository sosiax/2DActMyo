#       Makefile
#
#       Copyright 2009 Unknown <fons@arnold.cbm.uam.es>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
# 


####################
# Progam files

#OBJ= mm_data.o mm_pbsa_statistics.o
OBJ= 2DActMyo.o
PROG= 2DActMyo

####################
# Compiler flags

# Intel - Production
#~ CC=icpc
#CFLAGS= -Wall -static -O3 -xSSE4.1 -ipo -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
#~ CFLAGS= -O3 -mp1 -static

# GNU - Debuggin
CC=gcc
#CFLAGS= -Wall -std=c99 -pedantic -static -g -pg -fbounds-check -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
CFLAGS= -Wall -g -fbounds-check -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -DCHECK
LDFLAGS= -lm -lbsd 

# GNU - Production
#CFLAGS= -Wall -std=c99 -pedantic -static -O3 -march=nocona -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
#~ CFLAGS= -Wall  -static -O3 -march=nocona -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
#~ LDFLAGS= -lm -lbsd

all: $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o $(PROG)

%.o: %.c
	$(CC) $(LDFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -fr $(OBJ) $(PROG)

