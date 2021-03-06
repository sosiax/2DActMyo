/*
 * 2DActMyo.h
 * 
 * Copyright 2018 Pau Casanova <pau.casanova25@gmail.com>
 *                Alfonso Nunez Salgado <anunez@cbm.csic.es>
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
#include <stdlib.h>
#include <bsd/stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define ClearScreen() printf("\033[H\033[J")

#define DEB_PRINT(fmt, args...) fprintf(stderr, "%s:%d:%s(): " fmt, \
    __FILE__, __LINE__, __func__, ##args)

#define Perpendicular(M1,M2) ((M1.x*M2.x+M1.y*M2.y)==0)

//~ Molecule Type
#define NTYPE 4  //Nmber of types of molecules
#define NONE 0  //None
#define ACT  1  //Actine
#define CRL  2  //Crosslinker
#define MYO  3  //Myosine


typedef struct __molecule__ {
  int type;      //Molecule type
  int x;         //Direction X
  int y;         //Direction Y
  int nl[NTYPE]; //number of links per type
  int atp;       //presence or not of atp
  char ldir;     //link direction
  double E;      //energy
  struct __molecule__ * front;
  struct __molecule__ * rear;
  struct __molecule__ * left;
  struct __molecule__ * right;
}MOL;


typedef struct __membrane__ {
  struct __molecule__ ** P;      //Molecule matrix
  int Lx;                        //Longitud Sistema
  int Ly;                        //Anchura Sistema
  double D;                       //Relacion Reservorio/Sistema
  int Term;                    //Pasos de termalizacion
  int Teq;               //Pasos de simulacion
  double VR;                  //Numero de pasos internos de montecarlo (VR*Lx*Ly)
  int Cheese;             //Periodicidad toma de datos

  //Parametros de accion
  double pHid;                 //Probabilidad Hidrolisis
  double pAbs;                 //Probabilidad Absorcion
  double fHid;                    //Fraccion que intenta Hidrolisis
  double fAbs;                    //Fraccion que intenta Absorcion
  double TBreak;                 //Tension de rotura
  double fTDis;                 //Fraccion de la tension que se reparte entre los vecinos

  //Numero de proteinas total
  int NMAX[NTYPE];            
  int N[NTYPE];            

  //Parametros del systema
  double e[NTYPE];                   //Energia de enlace para cada molecula -> Cte
  double muE[NTYPE];                 //Barrera de potencial para cada molecula -> Cte
  double mu[NTYPE];                   //Energia total del sistema
  double Et;                         //Energia total del sistema

  time_t seed;                  //Semilla para random        
}MEMBRANE;

int InitMolecule(MOL * M);
/****
 * P *** int : Membrane 
 * P[x][y][0]: Molecule Type
 * P[x][y][1]: Molecule direction
 * P[x][y][2]: ATP presence [0,1]
 * P[x][y][3]: Link presence :  0
 * 				1
 * 				2
 * 				3
 * P[x][y][4]: Actine list tag [1,NFA]
 **/
//~ P Type
#define MTYPE 0
#define DIR   1
#define ATP   2
#define LINK  3
#define TAG   4
 
 
/**
 * FA ** int : Actine filamentent data
 * 
 * *FA[0] : General data
 * FA [0][0] : Number of filaments (NFA)
 * FA [0][1] : total sum length of filaments
 * FA [0][2] : total Crosslinkers connections
 * FA [0][3] : total Myosine connections
 * 
 * *FA[1..NFA] : filament data
 * FA[1..NFA][0] : Initial x 
 * FA[1..NFA][1] : Initial y 
 * FA[1..NFA][2] : Filement direction 
 * FA[1..NFA][3] : Filement lenght 
 **/
//~ FA Globals
#define GLOBAL 0
//~ FA types
#define NFACT   0
#define SLEN  1
#define NCRL  2
#define NMYO  3

//FA types
#define X  0
#define Y  1
#define FDIR  2
#define FLEN  3

 /**
  * FM ** int : Myosine  data
  * 
  **/
  
#define CRT 3
#define LLEFT  0
#define LRIGHT 1

/*
 * 
 * name: desconocido
 * @param
 * @return
 * 
 */
int ReadConfig(MEMBRANE * M);
int InitMolecule(MOL * M);
int Dynamics(MEMBRANE * M);
int PrintOutP(MEMBRANE M);
