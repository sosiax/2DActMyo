/*
 * kk.c
 * 
 * Copyright 2018 Alfonso Núñez Salgado <anunez@cbm.csic.es>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <bsd/stdlib.h>
#include <unistd.h>

int randomInt(int Upperbound);
float ramdomFloat();

int main(int argc, char **argv)
{
	int i=0;
	time_t seed=time(NULL);
	printf("RANDOM SEED %ld\n",seed); 
	srandom(1519076756);
	for (i=0;i<10;i++){
	  printf("RAND %d, %f ,%ld\n",randomInt(3),ramdomFloat(),time(NULL));
	  sleep(1);
	  }
	return 0;
}


int randomInt(int Upperbound){
   //~ return arc4random_uniform(Upperbound);
   return (int) random() % (Upperbound + 1);
}

float ramdomFloat(){         
   //~ return (float)arc4random() / UINT32_MAX;
   return (float) random() / RAND_MAX;
}
