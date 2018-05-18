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

/****
 * P *** int : Membrane 
 * P[x][y][0]: Molecule Type
 * P[x][y][1]: Molecule direction
 * P[x][y][2]: ATP presence [0,1]
 * P[x][y][3]: Link presence [0,1]
 * P[x][y][4]: Actine list tag [1,NFA]
 **/
//~ P Type
#define MTYPE 0
#define DIR   1
#define ATP   2
#define LINK  3
#define TAG   4
 
//~ Molecule Type
#define NONE 0  //None
#define ACT  1  //Actine
#define CRL  2  //Crosslinker
#define MYO  3  //Myosine
 
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
  * *FA[0] : General data
  * FM[0][0] : Total Number of connected filaments
  * FM[0][1] : NA
  * FM[0][2] : NA
  * FM[0][3] : NA
  * 
  * *FA[1..NFM] : Myosine data
  * FA[1..NFA][0] : Number of connected filaments [0,1]
  * FA[1..NFA][1] : Number of connected filaments [0,1]
  * FA[1..NFA][2] : ATP
  * FA[1..NFA][3] : Chronometer
  **/
#define CRT 3

/*
 * 
 * name: desconocido
 * @param
 * @return
 * 
 */

int En(int Lx,int Ly,double eA,double eL,double eM,int x,int y,int*** P,double** E,int** FM,int C,int R,double* DESis);
int EnAA(int Lx,int Ly,double eA,int x,int y,int*** P,double** E,int C,int R,double* DESis);
double EnAL(int Lx,int Ly,double eL,int x,int y,int*** P,double** E,int C);
double EnAM(int Lx,int Ly,double eM,int x,int y,int*** P,double** E,int** FM,int C);
double EnLA(int Lx,int Ly,double eL,int x,int y,int*** P,double** E,int C,int R,double* DESis);
double EnMA(int Lx,int Ly,double eM,int x,int y,int*** P,double** E,int** FM,int C);void HATP(int Lx,int Ly,double eA,double eL,double eM,int*** P,double** E,int** FM,int x,int y,long* HA);
void ResV(int Lx,int Ly,double eA,double eL,double eM,double mu,int*** P,double** E,int** FM,int x,int y,long* RI,long* RO,int* N,int* NX,int tpro,int OldT,int OldS);
inline double RegMu(double muE,int N,int NX,int D);
void FilNew(int Lx,int Ly,int*** P,int x,int y,int** FA,int** FL,int** FM,double* TA,int t,double TBreak);
void FilOld(int Lx,int Ly,int OldT,int OldS,int*** P,int x,int y,int** FA,int** FL,int** FM,double* TA);
inline void Neighbor(int Lx, int Ly,int*** P,int** FA,int** FL,int** FM,int** NH);
void FilBreak(int Lx,int Ly,double eA,double eL,double eM,double fTDis,int nz,int*** P,double** E,int** FA,double* TA,int** FL,int** FM,int** NH);
void ControlBundle(int* DP,int* CDP,int Lx, int Ly,int*** P,int** FA,int** FL,int** FM,int** NH,int** BA,double* TA,int NFA,int Cheese,int TBreak);

int randomInt(int Upperbound);
float ramdomFloat();
int DireccionXY(int Dir, int * x, int * y);
