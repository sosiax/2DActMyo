/*
 * main.c
 * 
 * Copyright 2018 Pau Casanova <pau.casanova25@gmail.com>
 *           Alfonso Nunez Salgado <anunez@cbm.csic.es>
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

#include <stdio.h>
#include "2DActMyo.h"

/** 
 * name: main
 * 
 *   Execute membrane dynamics
 * 
 * @param 
 * @return O : OK
 *         1 : ERR
 */
int main(int argc, char **argv)
{
  MEMBRANE M;
  long int t;
  
  //Read config
  // TODO : create a function to read config from file
  ReadConfig(&M);
  
  // Check MEMBRANE struct : in case a specific membrane state can be loaded
  //~ if (CheckAllMol(M))
    //~ return 1;
  
  for (t=0; t<M.Teq; t++){
	//
    if (Dynamics(&M)){
      PrintOutP(M);
      fprintf(stderr,"Itearations: %ld / %d - %2.1f \n",t,M.Teq,(double)t*100/M.Teq);
      DEB_PRINT("ERROR: On dynamics!!!!\n");
      return 1;
      }
      
    // Print output each M.Cheese interation
    if (t%M.Cheese== 0){ 
      system("clear");
      PrintOutP(M);
      fprintf(stderr,"Itearations: %ld / %d - %2.1f \n",t,M.Teq,(double)t*100/M.Teq);
      }
    }
  
  //Print out last membrane state
  system("clear");
  PrintOutP(M);
  fprintf(stderr,"Itearations: %ld / %d - %2.1f \n",t,M.Teq,(double)t*100/M.Teq);
  return 0;
}

