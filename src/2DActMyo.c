/*
 * 2DActMyo.c
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

#include "2DActMyo.h"

/** 
 * name: Inbound
 * 
 *   Checks if the position x,y in inside the membrane limits
 * 
 * @param MOLECULE * M : molecule struct to be inicialized
 * @return 0 : OK
 *      1 : ERR
 */
static inline int Inbound(MEMBRANE M, int x, int y){
  return ((x >= M.Lx) || (y >= M.Ly) || (x < 0) || (y < 0))? 0 : 1;
}

/** 
 * name: UpdateMu
 * 
 *  Modify links asociated to the molecule x,y when removed.
 * 
 * @param MEMBRANE M
 *     int x, int y : Curren position of the molecule
 * 
 * @return 0 : OK
 *      1 : ERR
 */
static inline double UpdateMu(MEMBRANE M,int Mtype){
  return (M.N[Mtype]!=0 )? \
          M.muE[Mtype]-log((double)(M.D+1)*M.N[Mtype]/M.NMAX[Mtype]) : \
          M.muE[Mtype]-log((double)(M.D+1)*0.001/(M.NMAX[Mtype]+0.001));
}

/** 
 * name: InitMolecule
 * 
 *   Inicialize molecue struct
 * 
 * @param MOLECULE * M : molecule struct to be inicialized
 * @return 0 : OK
 *      1 : ERR
 */
int InitMolecule(MOL * M){
  int i;
  if (M==NULL)
   return 1;
  M->type=0;    //Molecule type
  M->x=0;      //Direction X
  M->y=0;      //Direction Y
  for (i=0;i<NTYPE;i++)M->nl[i]=0;   //number of links  
  M->atp=0;  //presence or not of atp
  M->E=0.0;  //energy
  M->front=NULL;
  M->rear=NULL;
  M->left=NULL;
  M->right=NULL;
  return 0;
}

/** 
 * name: PrintOutP
 * 
 *   print out the current P status
 * 
 * @param MEMBRANE M
 * @return O : OK
 *         1 : ERR
 */
int PrintOutP(MEMBRANE M){               
  int x,y;
  int s,i;
  for(y=0;y<M.Ly;y++){
    for(x=0;x<M.Lx;x++){
      switch(M.P[x][y].type){
        case ACT:
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y== 0))  fprintf(stderr,"\u2b05");
          if ((M.P[x][y].x== 0)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u2b06");
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y== 0))  fprintf(stderr,"\u27a1");
          if ((M.P[x][y].x== 0)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u2b07");
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u2b09");
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u2b08");
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u2b0a");
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u2b0b");
          break;                                                 
        case MYO:                                                
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y== 0))  fprintf(stderr,"\u21e6");         
          if ((M.P[x][y].x== 0)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u21e7");         
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y== 0))  fprintf(stderr,"\u21e8");         
          if ((M.P[x][y].x== 0)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u21e9");         
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u2b01");         
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u2b00");         
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u2b02");         
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u2b03");         
          break;                                                     
        case CRL :                                                   
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y== 0))  fprintf(stderr,"\u2190"); 
          if ((M.P[x][y].x== 0)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u2191"); 
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y== 0))  fprintf(stderr,"\u2192"); 
          if ((M.P[x][y].x== 0)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u2193"); 
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u2196"); 
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y== 1))  fprintf(stderr,"\u2197"); 
          if ((M.P[x][y].x== 1)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u2198"); 
          if ((M.P[x][y].x==-1)&&(M.P[x][y].y==-1))  fprintf(stderr,"\u2199"); 
          break;
        default :
          fprintf(stderr," ");
        }
      }
    //~ fprintf(stderr,"\t%4d\n",y);
    fprintf(stderr,"\t%4d\t",y);
    
    for(x=0;x<M.Lx;x++,s=0){
      for(i=1;i<NTYPE;i++) s+=M.P[x][y].nl[i];
      (s==0)? fprintf(stderr," "): fprintf(stderr,"%d",s);
      }
    
    fprintf(stderr,"\n");
    }
  
  //~ fprintf(stderr,"=============================================================\n");
  //~ for(y=0;y<M.Ly;y++){
    //~ for(x=0;x<M.Lx;x++,s=0){
      //~ for(i=1;i<NTYPE;i++) s+=M.P[x][y].nl[i];
      //~ fprintf(stderr,"%d",s);
      //~ }
     //~ fprintf(stderr,"\t%4d\n",y);
    //~ }
  
  //~ for(x=0;x<M.Lx;x++)
    //~ fprintf(stderr,"%4d",x);
  //~ fprintf(stderr,"\n");

  fprintf(stderr,"=============================================================\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"ACT : %d/%d - CLR : %d/%d - MYO : %d/%d\n",M.N[ACT],M.NMAX[ACT],M.N[CRL],M.NMAX[CRL],M.N[MYO],M.NMAX[MYO]);
  fprintf(stderr,"ENERGY : %f\n",M.Et);
  return 0;
}

/** 
 * name: randomInt
 * 
 *   returns a random int in [0;Upperbound]
 * 
 * @param int Upperbound 
 * @return returns a random int in [0;Upperbound]
 */
int randomInt(int Upperbound){
  //~ return arc4random_uniform(Upperbound);
  return (int) random() % (Upperbound);
  //~ return (int) random() % (Upperbound + 1); --> deberia llevar el +1!
}

/** 
 * name: randomFloat
 * 
 *   returns a random double in [0;1]
 * 
 * @param none
 * @return returns a random double in [0;1]
 */
double randomFloat(){      
  //~ return (double)arc4random() / UINT32_MAX;
  return (double) random() / RAND_MAX;
}

/** 
 * name: randomType
 * 
 *   returns a random int molecule type in [0;NTYPE] in concordance 
 * whith the molecules in the reservoir.
 * 
 * @param MEMBRANE M
 * @return random int molecule type in [0;NTYPE]
 */
int randomType(MEMBRANE M){
  int nt=0;
  double p=randomFloat();
  double pm=0;
  int i=0;
  
  for (i=1;i<NTYPE; i++)
    nt+=M.NMAX[i]-M.N[i];
  
  for (i=1;i<NTYPE; i++){
    pm+=(M.NMAX[i]-M.N[i])/(double)nt;
    if (p < pm)
      return i;
    }
    
  return 0;
}

/** 
 * name: ReadConfig
 * 
 *   Allocate memory and inicialize simulation parameters
 * 
 * TODO:
 *   Read parameters from file
 * 
 * @param MEMBRANE * M : membrane struct to be inicialized
 * @return 0 : OK
 *      1 : ERR
 */
int ReadConfig(MEMBRANE * M){
  int x,y;
  
  if (M==NULL)
   return 1;

  M->Lx=90;              //Longitud Sistema
  M->Ly=70;              //Anchura Sistema
  M->D=1;                //Relacion Reservorio/Sistema
  M->Term=0;              //Pasos de termalizacion
  M->Teq=1000000000;          //Pasos de simulacion
  M->VR=0.5;              //Numero de pasos internos de montecarlo (VR*Lx*Ly)
  M->Cheese=1000000;          //Periodicidad toma de datos

  //Parametros de accion
  M->pHid=0.05;            //Probabilidad Hidrolisis
  M->pAbs=0.05;            //Probabilidad Absorcion
  M->fHid=2;              //Fraccion que intenta Hidrolisis
  M->fAbs=2;              //Fraccion que intenta Absorcion
  M->TBreak=50;            //Tension de rotura
  M->fTDis=0.6;            //Fraccion de la tension que se reparte entre los vecinos
  
  //Numero de proteinas total
  M->NMAX[NONE]=0;
  M->NMAX[ACT]=3000;             //Actinas
  M->NMAX[CRL]=000;             //Cross-Linkers
  M->NMAX[MYO]=800;             //Miosinas
  
  //Numero de proteinas en membrana
  M->N[NONE]=0;
  M->N[ACT]=0;             //Actinas
  M->N[CRL]=0;             //Cross-Linkers
  M->N[MYO]=0;             //Miosinas

  //Parametros del systema
  M->e[NONE]=0;
  M->e[ACT]=(double) 8;               //Energia de enlace Actina
  M->e[CRL]=(double) 2.5;               //Energia de enlace CrossLinker
  M->e[MYO]=(double) 2.5;               //Energia de enlace Myosina
  M->muE[NONE]=0;
  M->muE[ACT]=(double) -11;            //Barrera de potencial Actina
  M->muE[CRL]=(double) -6;             //Barrera de potencial CrossLinker
  M->muE[MYO]=(double) -6;             //Barrera de potencial Myosina
  
  
  M->P=(struct __molecule__ **)(malloc(M->Lx*sizeof (struct __molecule__ *)));
  if(M->P==NULL){
   DEB_PRINT("Error. Allocation was unsuccessful.\n");
   return 1;
   }

  for(x=0;x<M->Lx;x++){
    M->P[x]=(struct __molecule__ *)malloc(M->Ly*sizeof(struct __molecule__));
    if(M->P[x]==NULL){
     DEB_PRINT("Error. Allocation was unsuccessful.\n");
     return 1;
     }
    }

  for(x=0;x<M->Lx;x++)
    for(y=0;y<M->Ly;y++)
      InitMolecule(&(M->P[x][y]));

  //Potencial quimico inicial
  M->mu[NONE]=0;
  M->mu[ACT]=UpdateMu(*M,ACT);        //Actinas
  M->mu[CRL]=UpdateMu(*M,CRL);        //Cross-Linkers
  M->mu[MYO]=UpdateMu(*M,MYO);        //Miosinas
  M->Et=0;                   //Energia total del sistema
  
  //random seed -> set to 0 for diff random sequences
  //~ M->seed=time(NULL);
  M->seed=151907839;
  srandom(M->seed);
  DEB_PRINT("RANDOM SEED %ld\n",M->seed); 
  
  return 0;
}

/** 
 * name: CheckMol
 * 
 *   Function made for debug purposes that checks molecule connections.
 * 
 * @param MEMBRANE * M : membrane struct 
 * 		  int       x,y: molecule coordinates
 * 
 * @return 0 : OK
 *      1 : ERR
 */
int CheckMol(MEMBRANE M,int x, int y){
  int dx=M.P[x][y].x;
  int dy=M.P[x][y].y;
  int N[NTYPE];
  int i;
  int nx,ny;
  int ret=0;
  
  for(i=0;i<NTYPE; i++) N[i]=0;
  
  switch (M.P[x][y].type){
    case ACT:
      //check rear
      nx=x-dx; ny=y-dy;
      if (Inbound(M,nx,ny)){
        if ((M.P[nx][ny].type == ACT) &&     //Actine (-dx,-dy)
            (M.P[nx][ny].x == M.P[x][y].x) &&  //Same dir
            (M.P[nx][ny].y == M.P[x][y].y)){
          N[ACT]++;
          if (M.P[x][y].rear != &(M.P[nx][ny])){
            //~ DEB_PRINT("ERROR: Check ACT[%d][%d] link REAR not updated\n",x,y);
            ret=1;
            }
          }
	      else
          if (M.P[x][y].rear != NULL){
            DEB_PRINT("ERROR: Check ACT link REAR Shuld be NULL\n");
            ret=2;
            }    
          
        }
      
      //check left
      nx=x-dy; ny=y+dx;
      if (Inbound(M,nx,ny)){ 
        if ((M.P[nx][ny].type >ACT ) &&      //MYO (-dy,dx) 
            ((M.P[nx][ny].x*M.P[x][ny].x+M.P[nx][ny].y*M.P[x][y].y) == 0)){ //Perpendicular
	        N[M.P[nx][ny].type]++;
	        if (M.P[x][y].left != &(M.P[nx][ny])){
	          //~ DEB_PRINT("ERROR: Check ACT link LEFT to %d not updated\n",M.P[nx][ny].type);
	          ret=1;
	          }
	        }
        else
          if (M.P[x][y].left != NULL){
            DEB_PRINT("ERROR: Check ACT link LEFT Shuld be NULL\n");
            ret=2;
            }
          
        }

      //check right
      nx=x+dy; ny=y-dx;
      if (Inbound(M,nx,ny)){ 
        if ((M.P[nx][ny].type >ACT ) &&      //MYO (-dy,dx) 
             Perpendicular(M.P[nx][ny],M.P[x][y])){ //Perpendicular
	        N[M.P[nx][ny].type]++;
	        if (M.P[x][y].right != &(M.P[nx][ny])){
	          //~ DEB_PRINT("ERROR: Check ACT link RIGHT to %d not updated\n",M.P[nx][ny].type);
	          ret=1;
	          }
	        }
        else
          if (M.P[x][y].right != NULL){
            DEB_PRINT("ERROR: Check ACT link RIGHT Shuld be NULL\n");
            ret=2;
            }
        }
      break;
    
    case MYO:
    case CRL:
      //check front
      nx=x+dx; ny=y+dy;
      if (Inbound(M,nx,ny)){
        if ((M.P[nx][ny].type == ACT) &&     //Actine (-dx,-dy)
            ((M.P[nx][ny].x*M.P[x][y].x + M.P[nx][ny].y*M.P[x][y].y) == 0)){ //Perpendicular
          N[M.P[x][y].type]++;
          if (M.P[x][y].front != &(M.P[nx][ny])){
            //~ DEB_PRINT("ERROR: Check CRL link FRONT not updated\n");
            ret=1;
            }
	        }
	      else
          if (M.P[x][y].front != NULL){
            DEB_PRINT("ERROR: Check CRL link FRONT Shuld be NULL\n");
            ret=2;
            }
        }
      
      //check rear
      nx=x-dx; ny=y-dy;
      if (Inbound(M,nx,ny)){
        if ((M.P[nx][ny].type == ACT) &&     //Actine (-dx,-dy)
            ((M.P[nx][ny].x*M.P[x][y].x + M.P[nx][ny].y*M.P[x][y].y) == 0)){ //Perpendicular
          N[M.P[x][y].type]++;
          if (M.P[x][y].rear != &(M.P[nx][ny])){
            //~ DEB_PRINT("ERROR: Check CRL link REAR not updated\n");
            ret=1;
            }
          }
	      else 
          if (M.P[x][y].rear != NULL){
            DEB_PRINT("ERROR: Check CRL link REAR Shuld be NULL\n");
            ret=2;
            }
        }
    
      break;
    } //end switch
      
  for (i=0;i<NTYPE;i++){
    if (M.P[x][y].nl[i] > 2){
      DEB_PRINT("FATAL ERROR!!! : Number of links is bad nl[%d]=%d!!!!\n",i,M.P[x][y].nl[i]);
      ret=2;
      }
    if (M.P[x][y].nl[i] != N[i]){
      //~ if (N[i] < M.P[x][y].nl[i]){
        //~ DEB_PRINT("ERROR : Number of links seams bad nl[%i]=%d > %d ...\n",i,M.P[x][y].nl[i],N[i]);
        //~ ret=2;
        //~ }
      //~ else
        ret=1;
      }
    }
  
  if (ret>1)
    DEB_PRINT("MOL [%d][%d].type=%d nl[ACT]=%d nl[CRL]=%d nl[MYO]=%d\n", \
                      x,y,M.P[x][y].type,M.P[x][y].nl[ACT],M.P[x][y].nl[CRL],M.P[x][y].nl[MYO]);
  return ret;  
}

/** 
 * name: CheckAllMol
 * 
 *   Check all molecules in membrane
 * 
 * @param MEMBRANE M
 * @return 0: OK
 *         1: ERR
 */
int CheckAllMol(MEMBRANE M){
  int x,y;
  for (x=0;x<M.Lx;x++)
    for (y=0;y<M.Ly;y++)
      if (CheckMol(M,x,y)>1)
	      return 1;
  return 0;
}

/** 
 * name: IsThereLinkerLinks
 * 
 *   Returns how many links a crosslink have
 * 
 * @param MEMBRANE M
 *        x,y current crosslink coordinates
 * @return Number of links
 */
int IsThereLinkerLinks(MEMBRANE M,int x, int y){
  int nlinks=0;
  int dx=0;
  int dy=0;
  int nx,ny;
  
   
  dx=M.P[x][y].x;
  dy=M.P[x][y].y;
  
  if (M.P[x][y].type < ACT+1){ // **** addoc
    //DEB_PRINT("ERROR: This function should be called whith another molecule type!! (%d)\n",M.P[x][y].type);
    return 0;
    }
  
  nx=x+dx;
  ny=y+dy;
  //Checking front
  if ( (Inbound(M,nx,ny))&&
       (M.P[nx][ny].type == ACT) &&
      ((M.P[nx][ny].x*M.P[x][y].x + M.P[nx][ny].y*M.P[x][y].y) == 0)){//Perpendicular
      nlinks+=1;
      }
  
  nx=x-dx;
  ny=y-dy;
  if ((Inbound(M,nx,ny))&&
      (M.P[nx][ny].type == ACT) &&
      ((M.P[nx][ny].x*M.P[x][y].x + M.P[nx][ny].y*M.P[x][y].y) == 0)){//Perpendicular
      if (nlinks > 0){
        if ((M.P[nx][ny].x == -M.P[x+dx][y+dy].x) &&  // Diff sense
            (M.P[nx][ny].y == -M.P[x+dx][y+dy].y))
          nlinks+=1;
        else
          nlinks-=1;
        }
      else
        nlinks+=1;
      }
  return nlinks;
}


/** 
 * name: ThereIsFront
 * 
 *   Depending on the molecule it return modify lx and ly if there is a 
 * possible bond at the front of the molecule.
 * 
 * @param MEMBRANE M
 *     int x, int y : Curren position of the molecule
 *     int lx, int ly : position of the rear molecule if bond is possible
 * 
 * @return O : False
 *         1 : True
 */
int ThereIsCentral(MEMBRANE M,int x, int y, int * lx, int * ly, char D){
  int nl=0;
  int nx=(D=='F')?x+M.P[x][y].x:x-M.P[x][y].x;
  int ny=(D=='F')?y+M.P[x][y].y:y-M.P[x][y].y;
  
  *lx=nx;
  *ly=ny;
  
  if (Inbound(M,nx,ny)==0)
    return 0;
  
  switch(M.P[x][y].type){
    case ACT:
      if ((M.P[nx][ny].type == ACT) && 
	  (M.P[nx][ny].x == M.P[x][y].x) &&
          (M.P[nx][ny].y == M.P[x][y].y))
	  return 1;
      break;

    case CRL:
    case MYO:
	if ((M.P[nx][ny].type == ACT) &&
	    Perpendicular(M.P[nx][ny],M.P[x][y]) &&
	    (IsThereLinkerLinks(M,x,y) > 0))
	    return 1;
        break;

    default :
      DEB_PRINT("ERROR: no contempled molecule type (%d)!!",M.P[x][y].type);
    }
 
  return 0;//False
}

/** 
 * name: ThereIsLateral
 * 
 *   Depending on the molecule it return modify lx and ly if there is a 
 * possible bond in the latera of the molecule depending on parameter 'D'
 * 
 * @param MEMBRANE M
 *     int x, int y : Curren position of the molecule
 *     int lx, int ly : position of the rear molecule if bond is possible
 *     char D : 'L' or 'R' depending on the lateral to examine
 * 
 * @return O : False
 *         1 : True
 */  
int ThereIsLateral(MEMBRANE M,int x, int y, int * lx, int * ly,char D){
   int nx,ny=0;
   int nl;
   
   if (M.P[x][y].type != ACT){
      fprintf(stderr,"ERROR No se contempla lateral para otra molecula distinta de ACT: P[%d][%d][MTYPE]=%d\n",x,y,M.P[x][y].type);
      return 0;
      }
      
  //If left
  nx=(D == 'L')?x-M.P[x][y].y :x+M.P[x][y].y;
  ny=(D == 'L')?y+M.P[x][y].x :y-M.P[x][y].x;

  if (Inbound(M,nx,ny) == 0)
    return 0;
 
  *lx=nx;
  *ly=ny;
  
  if ((M.P[nx][ny].type > ACT) &&
      Perpendicular(M.P[nx][ny],M.P[x][y]) &&
      (IsThereLinkerLinks(M,nx,ny) > 0))
      return 1;

  return 0; 
}

/** 
 * name: DeltaE
 * 
 *  returns the energy asociated to the molecule defined in x,y position
 * 
 * @param MEMBRANE M
 *     int x, int y : Curren position of the molecule
 * 
 * @return double : energy asociated to the molecule
 */  

double DeltaE(MEMBRANE M, int x, int y){
  double delta=0;
  int i,lx=0,ly=0;
  
  // Energy depends on number of links
  for (i=ACT;i<NTYPE;i++)
    delta+=M.P[x][y].nl[i]*2*M.e[i];
  
  // if no links
  if (delta==0){
    switch(M.P[x][y].type){
      case ACT:
       if (ThereIsCentral(M,x,y,&lx,&ly,'R'))
         delta-=2*M.e[ACT];
       
       if (ThereIsLateral(M,x,y,&lx,&ly,'L'))
         delta-=2*M.e[M.P[lx][ly].type];
       
       if (ThereIsLateral(M,x,y,&lx,&ly,'R'))
         delta-=2*M.e[M.P[lx][ly].type];
       break;
      
      case MYO:
      case CRL:
        if (ThereIsCentral(M,x,y,&lx,&ly,'F'))
          delta-=2*M.e[M.P[x][y].type];
        if (ThereIsCentral(M,x,y,&lx,&ly,'R'))
          delta-=2*M.e[M.P[x][y].type];
        break; 
      
      default :
        DEB_PRINT("ERROR: no contempled molecule type (%d)!!\n",M.P[x][y].type);
      }
    }
  return delta;
}

/** 
 * name: InsertMol
 * 
 *  Modify links asociated to the molecule x,y when inserted
 * 
 * @param MEMBRANE M
 *     int x, int y : Curren position of the molecule
 * 
 * @return 0 : OK
 *      1 : ERR
 */
double InsertMol(MEMBRANE * M, int x, int y){
  int lx=0,ly=0;
  int dx=M->P[x][y].x,dy=M->P[x][y].y;

#ifdef CHECK
  if (Inbound(*M,x,y)==0){
    DEB_PRINT("Extreme ERROR!!! x,y, not in bound!!!\n");
    }
#endif

  switch(M->P[x][y].type){
    case ACT:
      if (ThereIsCentral(*M,x,y,&lx,&ly,'R')){
        M->P[x][y].rear=&(M->P[lx][ly]);
        M->P[x][y].nl[ACT]++;
        M->P[lx][ly].front=&(M->P[x][y]);
        M->P[lx][ly].nl[ACT]++;
        }
      
      if (ThereIsLateral(*M,x,y,&lx,&ly,'L')){
        M->P[x][y].left=&(M->P[lx][ly]);
        M->P[x][y].nl[M->P[lx][ly].type]++;
        M->P[lx][ly].nl[M->P[lx][ly].type]++;
        if ((lx+M->P[lx][ly].x==x)&&(ly+M->P[lx][ly].y==y)) // I'm it's front?
          M->P[lx][ly].front=&(M->P[x][y]);
        else
          M->P[lx][ly].rear=&(M->P[x][y]);
        }

      if (ThereIsLateral(*M,x,y,&lx,&ly,'R')){
        M->P[x][y].right=&(M->P[lx][ly]);
        M->P[x][y].nl[M->P[lx][ly].type]++;
        M->P[lx][ly].nl[M->P[lx][ly].type]++;
        if ((lx+M->P[lx][ly].x==x)&&(ly+M->P[lx][ly].y==y)) // I'm it's front?
          M->P[lx][ly].front=&(M->P[x][y]);
        else
          M->P[lx][ly].rear=&(M->P[x][y]);
        }
      break;
    
    case MYO:
    case CRL:
      if (ThereIsCentral(*M,x,y,&lx,&ly,'F')){
        M->P[x][y].nl[M->P[x][y].type]++;
        M->P[lx][ly].nl[M->P[x][y].type]++;
        M->P[x][y].front=&(M->P[lx][ly]);
        if ((lx-M->P[lx][ly].y == x) &&  //CRL is on the left of the actine
            (ly+M->P[lx][ly].x == y))
          M->P[lx][ly].left=&(M->P[x][y]);
        else
          M->P[lx][ly].right=&(M->P[x][y]);
        }
      if (ThereIsCentral(*M,x,y,&lx,&ly,'R')){
        M->P[x][y].nl[M->P[x][y].type]++;
        M->P[lx][ly].nl[M->P[x][y].type]++;
        M->P[x][y].rear=&(M->P[lx][ly]);
        if ((lx-M->P[lx][ly].y == x) &&  //CRL is on the left of the actine
            (ly+M->P[lx][ly].x == y))
          M->P[lx][ly].left=&(M->P[x][y]);
        else
          M->P[lx][ly].right=&(M->P[x][y]);
        }
      break;
    default :
      DEB_PRINT("ERROR: no contempled molecule type (%d)!!",M->P[x][y].type);
    }
  M->N[M->P[x][y].type]++; //Molecule count

#ifdef CHECK
  if (CheckMol(*M,x,y)>1)
    return 1;
#endif
  return 0;
}

/** 
 * name: RemoveMol
 * 
 *  Modify links asociated to the molecule x,y when removed.
 * 
 * @param MEMBRANE M
 *     int x, int y : Curren position of the molecule
 * 
 * @return 0 : OK
 *      1 : ERR
 */
double RemoveMol(MEMBRANE * M, int x, int y){
  int dx=M->P[x][y].x,dy=M->P[x][y].y;
  switch(M->P[x][y].type){
    case ACT:
    
      if (M->P[x][y].rear != NULL){
        //~ DEB_PRINT("Cannot remove rear ... ");
        //~ return 0;
        //~ }
        M->P[x][y].nl[ACT]--;
        M->P[x][y].rear->nl[ACT]--;
        M->P[x][y].rear->front=NULL;
        M->P[x][y].rear=NULL;
        }

      if (M->P[x][y].front != NULL){
        M->P[x][y].nl[ACT]--;
        M->P[x][y].front->nl[ACT]--;
        M->P[x][y].front->rear=NULL;
        M->P[x][y].front=NULL;
        }
        
      if (M->P[x][y].left != NULL){
        M->P[x][y].left->nl[M->P[x][y].left->type]--;
        M->P[x][y].nl[M->P[x][y].left->type]--;
        if (M->P[x][y].left->front==&(M->P[x][y])) // I'm it's front?
          M->P[x][y].left->front=NULL;
        else if (M->P[x][y].left->rear==&(M->P[x][y]))
          M->P[x][y].left->rear=NULL;
        else DEB_PRINT("ERROR : No correct pointer found on ACT LEFT\n");
        M->P[x][y].left=NULL;
        }
      if (M->P[x][y].right != NULL){
        M->P[x][y].right->nl[M->P[x][y].right->type]--;
        M->P[x][y].nl[M->P[x][y].right->type]--;
        if (M->P[x][y].right->front==&(M->P[x][y])) // I'm it's front?
          M->P[x][y].right->front=NULL;
        else if (M->P[x][y].right->rear==&(M->P[x][y]))
          M->P[x][y].right->rear=NULL;
        else DEB_PRINT("ERROR : No correct pointer found on ACT RIGHT\n");
        M->P[x][y].right=NULL;
        }
      break;
    
    case MYO:
    case CRL:
      if (M->P[x][y].front != NULL){
        M->P[x][y].nl[M->P[x][y].type]--;
        M->P[x][y].front->nl[M->P[x][y].type]--;
        if (M->P[x][y].front->left==&(M->P[x][y])) // I'm it's left?
          M->P[x][y].front->left=NULL;
        else if (M->P[x][y].front->right==&(M->P[x][y]))
          M->P[x][y].front->right=NULL;
        else DEB_PRINT("ERROR : No correct pointer found on CRL FRONT\n");
        M->P[x][y].front=NULL;
        }
      if (M->P[x][y].rear != NULL){
        M->P[x][y].nl[M->P[x][y].type]--;
        M->P[x][y].rear->nl[M->P[x][y].type]--;
        if (M->P[x][y].rear->left==&(M->P[x][y])) // I'm it's left?
          M->P[x][y].rear->left=NULL;
        else if (M->P[x][y].rear->right==&(M->P[x][y]))
          M->P[x][y].rear->right=NULL;
        else DEB_PRINT("ERROR : No correct pointer found on CRL REAR\n");
        M->P[x][y].rear=NULL;
        }
      break;
    
    default :
      DEB_PRINT("ERROR: no contempled molecule type (%d)!!",M->P[x][y].type);
    }
  M->N[M->P[x][y].type]--; //Molecule count
#ifdef CHECK
  if (M->P[x][y].nl[ACT]>1){
      DEB_PRINT("ERROR: too much ACT links\n");
      return 1;
     }
#endif
  return 0;
}

/** 
 * name: DirectionXY
 * 
 *  Transform integer code to x,y direction
 * 
 * @param MEMBRANE M
 *     int * x, int * y : variables where result will be stored
 * 
 * @return 0 : OK
 *      1 : ERR
 */
int DirectionXY(int Dir, int * x, int * y){
  switch(Dir) {
    case 1 :  *x= 0; *y= 1; break;
    case 2 :  *x= 1; *y= 0; break;
    case 3 :  *x= 0; *y=-1; break;
    case 4 :  *x=-1; *y= 0; break;
    //Not implemeted **** addoc
    case 5 :  *x=-1; *y=-1; DEB_PRINT("Oooh\n");break;
    case 6 :  *x=-1; *y= 1; break;
    case 7 :  *x=1 ; *y=-1; break;
    case 8 :  *x=1 ; *y= 1; break;
    
    default : 
       fprintf(stderr,"Error en el codigo de direccion: %d\n",Dir);
       return 1;
    }
  return 0;
}
/** 
 * name: Dynamics
 * 
 *   Execute one step of the membrane dynamics
 * 
 * @param MEMBRANE M
 * @return O : OK
 *         1 : ERR
 */
int Dynamics(MEMBRANE * M){
  int x=randomInt(M->Lx),y=randomInt(M->Ly);
  double p=randomFloat();
  double delta=0.0;
  double Exp=0.0;
  double Dt=0.0;
  
  if (M->P[x][y].type==NONE){
    //random Type and direction
    M->P[x][y].type=randomType(*M);

    if (M->P[x][y].type == 0){
      //~ DEB_PRINT("Warning: randomType must return a type, there is no enough molecules in the reservoir\n");
      return 0;
      }      

    DirectionXY(randomInt(4)+1,&(M->P[x][y].x),&(M->P[x][y].y));
    
    delta=DeltaE(*M,x,y);  // signo - cuando se crean enlaces
    Exp=exp(-delta+M->mu[M->P[x][y].type]);
    if(p<Exp){
      M->Et+=delta;
      if (InsertMol(M,x,y)){
        DEB_PRINT("ERROR: Inserting molecule\n");
	      return 1;
        }
        
#ifdef CHECK
      Dt=DeltaE(*M,x,y);
      if (delta != Dt){
        DEB_PRINT("ERROR: delta E not consistent!!!\n");
	      return 1;
        }
      if (CheckAllMol(*M)){
        DEB_PRINT("ERROR: Checking ALL molecules\n");
        return 1;
        }
#endif
      }
    else{
      InitMolecule(&(M->P[x][y]));
#ifdef CHECK
      if (CheckAllMol(*M)){
        DEB_PRINT("ERROR: Checking ALL molecules\n");
	      return 1;
        }
#endif
      }
    }
  else {
    delta=DeltaE(*M,x,y);   // !!! cambio de signo!
    Exp=exp(-delta-M->mu[M->P[x][y].type]);
    if( p < Exp ){
      M->Et+=delta;
      RemoveMol(M,x,y);
#ifdef CHECK
      Dt=DeltaE(*M,x,y);
      if (delta != Dt){
        DEB_PRINT("ERROR: delta E not consistent!!!\n");
	      return 1;
        }
#endif
      InitMolecule(&(M->P[x][y]));
#ifdef CHECK
      if (CheckAllMol(*M)){
        DEB_PRINT("ERROR: Checking molecule\n");
	      return 1;
        }
#endif
      }
    }
  M->mu[M->P[x][y].type]=UpdateMu(*M,M->P[x][y].type);
  return 0;
}


