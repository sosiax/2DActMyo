/*
 * 2DActMyo.c
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
#include "2DActMyo.h"
#include <stdio.h>
#include <stdlib.h>
#include <bsd/stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>


int randomInt(int Upperbound){
   //~ return arc4random_uniform(Upperbound);
   //~ return (int) random() % (Upperbound + 1); --> deberia llevar el +1!
   return (int) random() % (Upperbound);
}

float ramdomFloat(){         
   //~ return (float)arc4random() / UINT32_MAX;
   return (float) random() / RAND_MAX;
}

int main()
{
   FILE* docC;             //Archivo historial de acciones
   FILE* docN;             //Num de prot en el sistema
   FILE* docNF;            //Num de fil en el sistema
   FILE* docmu;            //Archivo potencial quimico del sistema

   int n=0,t=0,x=0,y=0;   //Variables de bucles
   int a=0,nz=0;
   double p=0,pro=0;      //Variable random
                                                      //Datos del sistema:
   int Lx,Ly,D,Teq,Term,Cheese;                       //Size(Lx*Ly),Size ResV,It Eq,It Term,Freq FOTO
   double pHid,pAbs;                                  //Probabilidad Hidrolisis y Absorcion
   double fHid,fAbs;                                  //Fraccion Hidrolisis y Absorcion
   double TBreak,fTDis;                               //Tension rotura y Fraccion de tension distribuida
   double VR;
   int NA=0,NL=0,NM=0,NAX,NLX,NMX;                    //Numero de proteinas sistema/reservorio
   double eA,eL,eM;                                   //Energias de enlace
   double muA,muL,muM,muEA,muEL,muEM;                 //Potenciales quimicos total/intrinsico

   int*** P;               //Coordenadas de las red
   double** E;             //Energia particulas
   int** FA;               //Control de filamentos
   double* TA;
   int** FL;
   int** FM;
   int** NH;
   int** BA;
   int NFA,NFL,NFM;        //Memoria listas

   double ATM=0;           //Tiempo absorcion ATP
   int BOOM=0;             //Flag rotura

   int CDP=0,DP=0,alert=1,palert=0;            //Control de impresion
   int XL=0;
   int AeF=0,L=0;
   float* DL;

   int OldT=0;             //Memoria proteina previa
   int OldS=0;             //Memoria estado previo

   long RAI=0;             //Contadores de acciones Actina
   long RAO=0;
   long RLI=0;             //Contadores de acciones C-Linker
   long RLO=0;
   long RMI=0;             //Contadores de acciones Myosina
   long RMO=0;

   long AA=0;              //Contador de acciones del ATP
   long HA=0;
   long RM=0;

   int CNFA=0,CNFL=0,CNFM=0;     //Control Memoria Listas

   /* Random Seed */
   // Seg Fault 1519078387
   //~ time_t seed=time(NULL);
   time_t seed=1519078387;
   srandom(1519078387);
   printf("RANDOM SEED %ld\n",seed); 
   
//Parametros de la simulacion
   Lx=51;                     //Longitud Sistema
   Ly=51;                     //Anchura Sistema
   D=1;                       //Relacion Reservorio/Sistema
   Term=0;                    //Pasos de termalizacion
   Teq=1000000;               //Pasos de simulacion
   VR=0.5;                    //Numero de pasos internos de montecarlo (VR*Lx*Ly)
   Cheese=100000;               //Periodicidad toma de datos
   Cheese=100000;               //Periodicidad toma de datos

//Parametros de accion
   pHid=0.05;                 //Probabilidad Hidrolisis
   pAbs=0.05;                 //Probabilidad Absorcion
   fHid=2;                    //Fraccion que intenta Hidrolisis
   fAbs=2;                    //Fraccion que intenta Absorcion
   TBreak=50;                 //Tension de rotura
   fTDis=0.6;                 //Fraccion de la tension que se reparte entre los vecinos

//Numero de proteinas total
   NAX=600;                   //Actinas
   NLX=300;                   //Cross-Linkers
   NMX=100;                   //Miosinas

//Parametros del systema
   eA=8;                      //Energia de enlace Actina
   eL=3;                      //Energia de enlace CrossLinker
   eM=3;                      //Energia de enlace Myosina
   muEA=-11;                  //Barrera de potencial Actina
   muEL=-6;                   //Barrera de potencial CrossLinker
   muEM=-6;                   //Barrera de potencial Myosina

//Memoria reservada para las listas
   NFA=NAX/3+1;
   NFL=NLX+1;
   NFM=NMX+1;

//Reserva de memoria:
   E=(double**) (malloc(Lx*sizeof (double*)));
   P=(int***) (malloc(Lx*sizeof (int**)));
   for(x=0;x<Lx;x++){
      E[x]=(double*)malloc(Ly*sizeof(double));
      P[x]=(int**)malloc(Ly*sizeof(int*));
      for(y=0;y<Ly;y++){
          P[x][y]=(int*)malloc(5*sizeof(int));
      }
   }
   FA=(int**) (malloc(NFA*sizeof (int*)));
   for(x=0;x<NFA;x++)
   {
      FA[x]=(int*)malloc(4*sizeof(int));
   }
   TA=(double*) (malloc(NFA*sizeof (double)));

   FL=(int**) (malloc(NFL*sizeof (int*)));
   for(x=0;x<NFL;x++)
   {
      FL[x]=(int*)malloc(2*sizeof(int));
   }
   FM=(int**) (malloc(NFM*sizeof (int*)));
   for(x=0;x<NFM;x++)
   {
      FM[x]=(int*)malloc(4*sizeof(int));
   }
   NH=(int**) (malloc(NFA*sizeof (int*)));
   for(x=0;x<NFA;x++)
   {
      NH[x]=(int*)malloc(NFA*sizeof(int));
   }
   BA=(int**) (malloc((NFA/3)*sizeof (int*)));
   for(x=0;x<(NFA/3);x++)
   {
      BA[x]=(int*)malloc(NFA*sizeof(int));
   }

   if(Lx>Ly) XL=Lx;
   else XL=Ly;
   DL=(float*)malloc(XL*sizeof(float));

//Inicializacion potencial con el sistema vacio
   muA=RegMu(muEA,NA,NAX,D);
   muL=RegMu(muEL,NL,NLX,D);
   muM=RegMu(muEM,NM,NMX,D);

//Archivo de informacion sobre la simulacion
   docC=fopen("AA0Control.dat","w");
//Archivo de numero de proteinas
   docN=fopen("DataN.dat","w");
//Archivo de los potenciales quimicos
   docmu=fopen("DataMU.dat","w");
//Archivo de numero de actina en filamento
   docNF=fopen("DataNF.dat","w");

//Print de los parametros
   fprintf(docC,"Datos introducidos:\n\n");
   fprintf(docC,"Dimension de la red: %ix%i \n",Lx,Ly);
   fprintf(docC,"Dimension del reservorio: %ix%ix%i \n",D,Lx,Ly);
   fprintf(docC,"Numero de actinas: %i \n",NAX);
   fprintf(docC,"Numero de C-Linkers: %i \n",NLX);
   fprintf(docC,"Numero de myosinas: %i \n",NMX);
   fprintf(docC,"Energia de enlace de actinas: %lf \n",eA);
   fprintf(docC,"Energia de enlace actina-C-Linker: %lf \n",eL);
   fprintf(docC,"Energia de enlace actina-Myosina: %lf \n",eM);
   fprintf(docC,"Barrera energetica de las actinas: %lf \n",muEA);
   fprintf(docC,"Barrera energetica de las C-Linkers: %lf \n",muEL);
   fprintf(docC,"Barrera energetica de las myosinas: %lf \n",muEM);
   fprintf(docC,"Termalizacion  : %i \n\n",Term);
   fprintf(docC,"Iteraciones  : %i \n\n",Teq);
   fprintf(docC,"Intervalo entre imagenes: %i \n",Cheese);

   fprintf(docC,"Posibilidades de accion:\n");
   fprintf(docC,"Probabilidad de intentar que la actina hidrolize ATP: %lf \n\n",pHid);
   fprintf(docC,"Fraccion de la actina que intenta hidrolizar ATP: %lf \n",fHid);
   fprintf(docC,"Probabilidad de intentar que la myosina absorba ATP: %lf \n",pAbs);
   fprintf(docC,"Fraccion de la myosina que intenta absorber ATP: %lf \n",fAbs);
   fprintf(docC,"Tension de Rotura: %lf \n",TBreak);
   fprintf(docC,"Fraccion de tension propagada: %lf \n",fTDis);

//Inicializacion de el puntero de las posiciones del espacio y del de las energias
   for(x=0;x<Lx;x++)
   {
      for(y=0;y<Ly;y++)
      {
         E[x][y]=0;
         for(t=0;t<5;t++) P[x][y][t]=0;
      }
   }

//Inicializacion de los punteros de control de filamentos
   for(x=0;x<NFA;x++)
   {
      for(y=0;y<4;y++) FA[x][y]=0;
      TA[x]=0;
   }
   for(x=0;x<NFL;x++) for(y=0;y<2;y++) FL[x][y]=0;

   for(x=0;x<NFM;x++) for(y=0;y<4;y++) FM[x][y]=0;

   for(x=0;x<XL;x++) DL[x]=0;

//Inicializacion contadores
   DP=-1;
   CDP=Cheese-1;

   alert=Term/5;
   palert=alert;

//Termalizacion
   for(t=0;t<Term;t++)
   {
      if(t==alert)
      {
         printf("Termalizacion: %i/%i \n",alert/palert,Term/palert);
         alert+=palert;
      }

      for(n=0;n<(Lx*Ly);n++) //En cada paso de tiempo se realizan Lx*Ly acciones (slots del sistema)
      {
         x=randomInt(Lx);
         y=randomInt(Ly);

         OldT=P[x][y][0];
         OldS=P[x][y][3];

         p=ramdomFloat();

         if(P[x][y][0]==0)
         {
            pro=randomInt(3)+1;

            if(pro==1)
            {
               ResV(Lx,Ly,eA,eL,eM,muA,P,E,FM,x,y,&RAI,&RAO,&NA,&NAX,pro,OldT,OldS);
               muA=RegMu(muEA,NA,NAX,D);
               if((P[x][y][0]!=0)&&(P[x][y][3]!=0)) FilNew(Lx,Ly,P,x,y,FA,FL,FM,TA,t-Term,TBreak);
            }
            else if(pro==2)
            {
               ResV(Lx,Ly,eA,eL,eM,muL,P,E,FM,x,y,&RLI,&RLO,&NL,&NLX,pro,OldT,OldS);
               muL=RegMu(muEL,NL,NLX,D);
               if((P[x][y][0]!=0)&&(P[x][y][3]!=0)) FilNew(Lx,Ly,P,x,y,FA,FL,FM,TA,t-Term,TBreak);
            }
            else if(pro==3)
            {
               ResV(Lx,Ly,eA,eL,eM,muM,P,E,FM,x,y,&RMI,&RMO,&NM,&NMX,pro,OldT,OldS);
               muM=RegMu(muEM,NM,NMX,D);
               if((P[x][y][0]!=0)&&(P[x][y][3]!=0)) FilNew(Lx,Ly,P,x,y,FA,FL,FM,TA,t-Term,TBreak);
            }
         }
         else if(P[x][y][0]!=0)
         {
            if(P[x][y][0]==1)
            {
               ResV(Lx,Ly,eA,eL,eM,muA,P,E,FM,x,y,&RAI,&RAO,&NA,&NAX,P[x][y][0],OldT,OldS);
               muA=RegMu(muEA,NA,NAX,D);
               if((P[x][y][0]==0)&&(OldS!=0)&&(P[x][y][4]!=0)) FilOld(Lx,Ly,OldT,OldS,P,x,y,FA,FL,FM,TA);
            }
            else if(P[x][y][0]==2)
            {
               ResV(Lx,Ly,eA,eL,eM,muL,P,E,FM,x,y,&RLI,&RLO,&NL,&NLX,P[x][y][0],OldT,OldS);
               muL=RegMu(muEL,NL,NLX,D);
               if((P[x][y][0]==0)&&(OldS!=0)&&(P[x][y][4]!=0)) FilOld(Lx,Ly,OldT,OldS,P,x,y,FA,FL,FM,TA);
            }
            else if(P[x][y][0]==3)
            {
               ResV(Lx,Ly,eA,eL,eM,muM,P,E,FM,x,y,&RMI,&RMO,&NM,&NMX,P[x][y][0],OldT,OldS);
               muM=RegMu(muEM,NM,NMX,D);
               if((P[x][y][0]==0)&&(OldS!=0)&&(P[x][y][4]!=0)) FilOld(Lx,Ly,OldT,OldS,P,x,y,FA,FL,FM,TA);
            }
         }
         //Control/Revision memoria listas
         if(CNFA<FA[0][0]) CNFA=FA[0][0];
         if(CNFL<FL[0][0]) CNFL=FL[0][0];
         if(CNFM<FM[0][0]) CNFM=FM[0][0];
      }

      for(a=0;a<=(FA[0][0]*fHid);a++)
      {
         x=randomInt(Lx);
         y=randomInt(Ly);

         p=ramdomFloat();

         if((p<pHid)&&(P[x][y][0]==1)&&(P[x][y][3]==1))
         {
            OldT=P[x][y][0];
            OldS=P[x][y][3];
            HATP(Lx,Ly,eA,eL,eM,P,E,FM,x,y,&HA);
            if((P[x][y][3]==0)&&(P[x][y][4]!=0)) FilOld(Lx,Ly,OldT,OldS,P,x,y,FA,FL,FM,TA);
         }
      }

      for(a=0;a<=(FA[0][3]*fAbs);a++)
      {
         x=randomInt(Lx);
         y=randomInt(Ly);

         p=ramdomFloat();

         if((p<pAbs)&&(P[x][y][0]==3)&&(P[x][y][3]==1))
         {
            P[x][y][2]++;
            nz=P[x][y][4];
            FM[nz][2]=P[x][y][2];
            TA[FM[nz][0]]++;
            if(TA[FM[nz][0]]>TBreak) BOOM=1;
            TA[FM[nz][1]]++;
            if(TA[FM[nz][1]]>TBreak) BOOM=1;
            ATM+=(t-Term)-FM[nz][3];
            FM[nz][3]=t-Term;
            AA++;
         }
      }
      while(BOOM==1)
      {
         for(x=0;x<NFA;x++) for(y=0;y<NFA;y++) NH[x][y]=0;
         Neighbor(Lx,Ly,P,FA,FL,FM,NH);
         for(nz=1;nz<=TA[0];nz++) if(TA[nz]>TBreak) break;
         FilBreak(Lx,Ly,eA,eL,eM,fTDis,nz,P,E,FA,TA,FL,FM,NH);
         RM++;
         BOOM=0;
         for(nz=1;nz<=TA[0];nz++)
         {
            if(TA[nz]>TBreak)
            {
               BOOM=1;
               break;
            }
         }
      }
   }

   if (Term>0){
      printf("Termalizacion: %i/%i \n",Term/palert,Term/palert);
      printf("\nControl Memoria:\n\tFAX:%d/%d\tFLX:%d/%d\tFMX:%d/%d\n\n",CNFA,NFA-1,CNFL,NFL-1,CNFM,NFM-1);
      }

//Reinicializacion contadores
   alert=Teq/5;
   palert=alert;

   RAI=0;
   RAO=0;
   RLI=0;
   RLO=0;
   RMI=0;
   RMO=0;
   HA=0;
   RM=0;

//Simulacion
   for(t=0;t<Teq;t++)
   {
      if(t==alert)
      {
         printf("Simulacion: %i/%i \n",alert/palert,Teq/palert);
         alert+=palert;
      }

      AeF=0;
      L=0;
      for(x=0;x<XL;x++) DL[x]=0;

      for(a=1;a<=FA[0][0];a++)
      {
         if(FA[a][3]!=0)
         {
            DL[FA[a][3]]+=FA[a][3];
            AeF+=FA[a][3];
            L++;
         }
      }

      fprintf(docNF,"%i   %i   %i\n",t,L,AeF);

      fprintf(docN,"%i   %i   %i   %i   %i   %i   %i\n",t,NA,NL,NM,NAX,NLX,NMX);

      fprintf(docmu,"%i  %lf  %lf  %lf\n",t,muA,muL,muM);

      CDP++;
      ControlBundle(&DP,&CDP,Lx,Ly,P,FA,FL,FM,NH,BA,TA,NFA,Cheese,TBreak);

      for(n=0;n<(Lx*Ly)*VR;n++) //En cada paso de tiempo se realizan Lx*Ly*VR acciones (slots del sistema)
      {
         x=randomInt(Lx);
         y=randomInt(Ly);

         OldT=P[x][y][0];
         OldS=P[x][y][3];

         p=ramdomFloat();

         if(P[x][y][0]==0)
         {
            pro=randomInt(3)+1;

            if(pro==1)
            {
               ResV(Lx,Ly,eA,eL,eM,muA,P,E,FM,x,y,&RAI,&RAO,&NA,&NAX,pro,OldT,OldS);
               muA=RegMu(muEA,NA,NAX,D);
               if((P[x][y][0]!=0)&&(P[x][y][3]!=0)) FilNew(Lx,Ly,P,x,y,FA,FL,FM,TA,t,TBreak);
            }
            else if(pro==2)
            {
               ResV(Lx,Ly,eA,eL,eM,muL,P,E,FM,x,y,&RLI,&RLO,&NL,&NLX,pro,OldT,OldS);
               muL=RegMu(muEL,NL,NLX,D);
               if((P[x][y][0]!=0)&&(P[x][y][3]!=0)) FilNew(Lx,Ly,P,x,y,FA,FL,FM,TA,t,TBreak);
            }
            else if(pro==3)
            {
               ResV(Lx,Ly,eA,eL,eM,muM,P,E,FM,x,y,&RMI,&RMO,&NM,&NMX,pro,OldT,OldS);
               muM=RegMu(muEM,NM,NMX,D);
               if((P[x][y][0]!=0)&&(P[x][y][3]!=0)) FilNew(Lx,Ly,P,x,y,FA,FL,FM,TA,t,TBreak);
            }
         }
         else if(P[x][y][0]!=0)
         {
            if(P[x][y][0]==1)
            {
               ResV(Lx,Ly,eA,eL,eM,muA,P,E,FM,x,y,&RAI,&RAO,&NA,&NAX,P[x][y][0],OldT,OldS);
               muA=RegMu(muEA,NA,NAX,D);
               if((P[x][y][0]==0)&&(OldS!=0)&&(P[x][y][4]!=0)) FilOld(Lx,Ly,OldT,OldS,P,x,y,FA,FL,FM,TA);
            }
            else if(P[x][y][0]==2)
            {
               ResV(Lx,Ly,eA,eL,eM,muL,P,E,FM,x,y,&RLI,&RLO,&NL,&NLX,P[x][y][0],OldT,OldS);
               muL=RegMu(muEL,NL,NLX,D);
               if((P[x][y][0]==0)&&(OldS!=0)&&(P[x][y][4]!=0)) FilOld(Lx,Ly,OldT,OldS,P,x,y,FA,FL,FM,TA);
            }
            else if(P[x][y][0]==3)
            {
               ResV(Lx,Ly,eA,eL,eM,muM,P,E,FM,x,y,&RMI,&RMO,&NM,&NMX,P[x][y][0],OldT,OldS);
               muM=RegMu(muEM,NM,NMX,D);
               if((P[x][y][0]==0)&&(OldS!=0)&&(P[x][y][4]!=0)) FilOld(Lx,Ly,OldT,OldS,P,x,y,FA,FL,FM,TA);
            }
         }
         //Control/Revision memoria listas
         if(CNFA<FA[0][0]) CNFA=FA[0][0];
         if(CNFL<FL[0][0]) CNFL=FL[0][0];
         if(CNFM<FM[0][0]) CNFM=FM[0][0];
      }

      for(a=0;a<=(FA[0][0]*fHid);a++)
      {
         x=randomInt(Lx);
         y=randomInt(Ly);

         p=ramdomFloat();

         if((p<pHid)&&(P[x][y][0]==1)&&(P[x][y][3]==1))
         {
            OldT=P[x][y][0];
            OldS=P[x][y][3];
            HATP(Lx,Ly,eA,eL,eM,P,E,FM,x,y,&HA);
            if((P[x][y][3]==0)&&(P[x][y][4]!=0)) FilOld(Lx,Ly,OldT,OldS,P,x,y,FA,FL,FM,TA);
         }
      }

      for(a=0;a<=(FA[0][3]*fAbs);a++)
      {
         x=randomInt(Lx);
         y=randomInt(Ly);

         p=ramdomFloat();

         if((p<=pAbs)&&(P[x][y][0]==3)&&(P[x][y][3]==1))
         {
            P[x][y][2]++;
            nz=P[x][y][4];
            FM[nz][2]=P[x][y][2];
            TA[FM[nz][0]]++;
            if(TA[FM[nz][0]]>TBreak) BOOM=1;
            TA[FM[nz][1]]++;
            if(TA[FM[nz][1]]>TBreak) BOOM=1;
            ATM+=t-FM[nz][3];
            FM[nz][3]=t;
            AA++;
         }
      }
      while(BOOM==1)
      {
         for(x=0;x<NFA;x++) for(y=0;y<NFA;y++) NH[x][y]=0;
         Neighbor(Lx,Ly,P,FA,FL,FM,NH);
         for(nz=1;nz<=TA[0];nz++) if(TA[nz]>TBreak) break;
         FilBreak(Lx,Ly,eA,eL,eM,fTDis,nz,P,E,FA,TA,FL,FM,NH);
         RM++;
         BOOM=0;
         for(nz=1;nz<=TA[0];nz++)
         {
            if(TA[nz]>TBreak)
            {
               BOOM=1;
               break;
            }
         }
      }
   }


   printf("Simulacion: %i/%i \n",Teq/palert,Teq/palert);
   printf("\nControl Memoria:\n\tFAX:%d/%d\tFLX:%d/%d\tFMX:%d/%d\n\n",CNFA,NFA-1,CNFL,NFL-1,CNFM,NFM-1);

   AeF=0;
   L=0;
   for(a=1;a<=FA[0][0];a++)
   {
      if(FA[a][3]!=0)
      {
         DL[FA[a][3]]+=FA[a][3];
         AeF+=FA[a][3];
         L++;
      }
   }

   fprintf(docNF,"%i   %i   %i\n",t,L,AeF);

   fprintf(docN,"%i   %i   %i   %i   %i   %i   %i\n",t,NA,NL,NM,NAX,NLX,NMX);

   fprintf(docmu,"%i  %lf  %lf  %lf\n",t,muA,muL,muM);

   CDP++;
   ControlBundle(&DP,&CDP,Lx,Ly,P,FA,FL,FM,NH,BA,TA,NFA,Cheese,TBreak);

   ATM/=AA;

   fprintf(docC,"Acciones realizadas durante el equilibrio:\n\n");

   fprintf(docC,"Introducciones de actina aceptadas: %li \n",RAI);
   fprintf(docC,"Salidas de actina aceptadas: %li \n",RAO);
   fprintf(docC,"Introducciones de C-Linker aceptadas: %li \n",RLI);
   fprintf(docC,"Salidas de C-Linker aceptadas: %li \n",RLO);
   fprintf(docC,"Introducciones de myo aceptadas: %li \n",RMI);
   fprintf(docC,"Salidas de myo aceptadas: %li \n",RMO);
   fprintf(docC,"Absorciones de ATP (myosina): %li \n",AA);
   fprintf(docC,"Tiempo medio de absorcion de ATP por la myosina: %lf \n",ATM);
   fprintf(docC,"Hidrolisis de ATP (actina): %li \n",HA);
   fprintf(docC,"Roturas de filamento: %li \n",RM);

   fclose(docC);
   fclose(docN);
   fclose(docmu);
   fclose(docNF);

   for(x=0;x<Lx;x++)
   {
      free(E[x]);

      for(y=0;y<Ly;y++)
         free(P[x][y]);
   }

   free(E);

   for(x=0;x<Lx;x++)
      free(P[x]);

   free(P);

   for(x=0;x<NFA;x++)
   {
      free(FA[x]);
   }
   free(FA);

   free(TA);

   for(x=0;x<NFL;x++)
   {
      free(FL[x]);
   }
   free(FL);

   for(x=0;x<NFM;x++)
   {
      free(FM[x]);
   }
   free(FM);

   for(x=0;x<NFA;x++)
   {
      free(NH[x]);
   }
   free(NH);

   for(x=0;x<NFA/3;x++)
   {
      free(BA[x]);
   }
   free(BA);

   free(DL);

   return 0;
}

/*
   Explicacion de la informacion contenida en el marcador P:
   Limites:
      P[x][y][t] (x,y e [0,L) & t e [0,4])
   Contenido:
      P[x][y][0]---Tipo de particula en la posicion (x,y): 0=Vacia,1=Actina,2=C-Linker,3=Myosin
      P[x][y][1]---Direccion de la particula en la posicion (x,y):
         Actina: 1=Y+,2=X+,3=Y- y 4=X-
         C-Linker o Myosin: 1=X,2=Y
      P[x][y][2]---Presencia o no de APT:
         C-Linker: 0=No
         Actina: 0=No,1=Si
         Myosin: Numero de ATP absorbidos
      P[x][y][3]---Presencia o no de enlace:
         Actina,Myosin: 0=No y 1=Si
         C-Linker 0=No,1=X+,2=Y+,3=X-,4=Y-,5=X y 6=Y
      P[x][y][4]---Marcadores:
         Actina: N Filamento al que pertenece.
         C-Linker: Identificacion en la lista FL (solo si [3]=5|6)
         Myosin: Identificacion en la lista FM (solo si [3]=1)

   Explicacion de la informacion contenida en el marcador FA:
   Limites:
      FA[nf][a] (nf e [0,NFA) & a e [0,3])
   Contenido:
      FA[0][0]: Cota superior al num de filamentos en el sistema (NFA)
           [1]: Sumatorio de la longitud de los filamentos del sistema
           [2]: Sumatorio de C-Linkers que unen filamentos
           [3]: Sumatorio de Myosinas que unen filamentos
      FA[nf>0][0]: Coordenada x de la actina inicial
              [1]: Coordenada y de la actina inicial
              [2]: Direccion del filamento
              [3]: Longitud del filamento

   Explicacion de la informacion contenida en el marcador TA:
   Limites:
      TA[nf] (nf e [0,NFA))
   Contenido:
      TA[0]: Cota superior al num de filamentos en el sistema (NFA)
      TA[nf>0]: Tension en el filamento

   Explicacion de la informacion contenida en el marcador FL:
   Limites:
      FL[nl][b] (nl e [0,NFL) & b e [0,1])
   Contenido:
      FL[0][0]: Cota superior al num de C-Linkers que enlazan filamentos (NFL)
      FL[nl>0][b]: nf de los filamentos enlazados

   Explicacion de la informacion contenida en el marcador FM:
   Limites:
      FM[nm][d] (nm e [0,NFM) & d e [0,3])
   Contenido:
      FM[0][0]: Cota superior al num de C-Linkers que enlazan filamentos (NFM)
      FM[n>0][0&1]: nf de los filamentos enlazados
      FM[n>0][2]: ATP de la Myosina (Ya no es util con el ultimo cambio)
      FM[n>0][3]: Cronometro

   Explicacion de la informacion contenida en el marcador NH:
   Limites:
      NH[nh][c] (nh e [0,NFA] & c e [0,NFA))
   Contenido:
      NH[0][0]: Cota superior al num de filamentos en el sistema (NFA)
      NH[nb>0][0]: N filamentos adyacentes
      NH[nb>0][c>0]: Lista de los filamentos adyacentes

   Explicacion de la informacion contenida en el marcador BA:
   Limites:
      BA[a][b] (a e [0,NFA/2] & c e [0,NFA))
   Contenido:
      BA[0][0]: Cota superior al num de Bundles en el sistema (NFB)
      BA[a>0][0]: N filamentos en el bundle
      BA[a>0][b>0]: Lista de los filamentos adyacentes
*/
//Creacion de los enlaces entre moleculas
int En(int Lx,int Ly,double eA,double eL,double eM,int x,int y,int*** P,double** E,int** FM,int C,int R,double* DESis)
// C(=0 Romper enlaces,=1 Crear enlaces),R(Memoria: =1 AA(Direccion negativa))
{
   *DESis=0;

   if(C==0)
   {
      if(P[x][y][0]==1)
      {
         *DESis=*DESis+EnAM(Lx,Ly,eM,x,y,P,E,FM,C);
         *DESis=*DESis+EnAL(Lx,Ly,eL,x,y,P,E,C);
         R=EnAA(Lx,Ly,eA,x,y,P,E,C,0,DESis);
      }
      else if(P[x][y][0]==2) R=EnLA(Lx,Ly,eL,x,y,P,E,C,0,DESis);
      else if(P[x][y][0]==3) *DESis=EnMA(Lx,Ly,eM,x,y,P,E,FM,C);
   }
   else if(C==1)
   {
      if(P[x][y][0]==1)
      {
         EnAA(Lx,Ly,eA,x,y,P,E,C,R,DESis);
         *DESis=*DESis+EnAL(Lx,Ly,eL,x,y,P,E,C);
         *DESis=*DESis+EnAM(Lx,Ly,eM,x,y,P,E,FM,C);
      }
      else if(P[x][y][0]==2) EnLA(Lx,Ly,eL,x,y,P,E,C,R,DESis);
      else if(P[x][y][0]==3) *DESis=EnMA(Lx,Ly,eM,x,y,P,E,FM,C);
   }
   return R;
}

//Enlaces Actina-Actina
int EnAA(int Lx,int Ly,double eA,int x,int y,int*** P,double** E,int C,int R,double* DESis)
{
   int EP=0,a=0; //Dummies de energia y direccion

   if(C==1) EP=-eA;
   else if(C==0) EP=+eA;

   if((P[x][y][1]==1)||(P[x][y][1]==2)) a=1;          //Para acortar codigo asocio los mismos sentidos
   else if((P[x][y][1]==3)||(P[x][y][1]==4)) a=-1;

   if((P[x][y][1]==1)||(P[x][y][1]==3))
   {
      if((y-a<0)||(y-a>Ly-1)){}                                                  //Comprobacion de fronteras
      else if((P[x][y-a][0]==1)&&(P[x][y-a][1]==P[x][y][1])&&(P[x][y][2]==1)&&(  //Enlace en direccion de crecimiento
      (((R==1)||(R==3)||(R==4))&&(C==1)&&(((P[x][y-a][3]==0)||(P[x][y-a][3]==2))&&((P[x][y][3]==0)||(P[x][y][3]==1))))||
                                ((C==0)&&(((P[x][y-a][3]==1)||(P[x][y-a][3]==3))&&((P[x][y][3]==2)||(P[x][y][3]==3))))))
      {
         E[x][y-a]+=EP;
         E[x][y]+=EP;
         if(C==1)
         {
            if(P[x][y-a][3]==0) P[x][y-a][3]=1;
            else if(P[x][y-a][3]==2) P[x][y-a][3]=3;
            if(P[x][y][3]==0) P[x][y][3]=2;
            else if(P[x][y][3]==1) P[x][y][3]=3;
         }
         if(C==0)
         {
            if(P[x][y-a][3]==1) P[x][y-a][3]=0;
            else if(P[x][y-a][3]==3) P[x][y-a][3]=2;
            if(P[x][y][3]==2) P[x][y][3]=0;
            else if(P[x][y][3]==3) P[x][y][3]=1;

            R=1;
         }
         *DESis-=2*eA;
      }

      if((y+a<0)||(y+a>Ly-1)){}                                                    //Comprobacion de fronteras
      else if((P[x][y+a][0]==1)&&(P[x][y+a][1]==P[x][y][1])&&(P[x][y+a][2]==1)&&(  //Enlace en direccion contraria
      (((R==2)||(R==3))&&(C==1)&&(((P[x][y][3]==0)||(P[x][y][3]==2))&&((P[x][y+a][3]==0)||(P[x][y+a][3]==1))))||
                        ((C==0)&&(((P[x][y][3]==1)||(P[x][y][3]==3))&&((P[x][y+a][3]==2)||(P[x][y+a][3]==3))))))
      {
         E[x][y]+=EP;
         E[x][y+a]+=EP;
         if(C==1)
         {
            if(P[x][y][3]==0) P[x][y][3]=1;
            else if(P[x][y][3]==2) P[x][y][3]=3;
            if(P[x][y+a][3]==0) P[x][y+a][3]=2;
            else if(P[x][y+a][3]==1) P[x][y+a][3]=3;
         }
         if(C==0)
         {
            if(P[x][y][3]==1) P[x][y][3]=0;
            else if(P[x][y][3]==3) P[x][y][3]=2;
            if(P[x][y+a][3]==2) P[x][y+a][3]=0;
            else if(P[x][y+a][3]==3) P[x][y+a][3]=1;

            if(R==1) R=3;
            else R=2;
         }
         *DESis-=2*eA;
      }
   }
   else if((P[x][y][1]==2)||(P[x][y][1]==4))
   {
      if((x-a<0)||(x-a>Lx-1)){}                                                  //Comprobacion de fronteras
      else if((P[x-a][y][0]==1)&&(P[x-a][y][1]==P[x][y][1])&&(P[x][y][2]==1)&&(  //Enlace en direccion de crecimiento
      (((R==1)||(R==3)||(R==4))&&(C==1)&&(((P[x-a][y][3]==0)||(P[x-a][y][3]==2))&&((P[x][y][3]==0)||(P[x][y][3]==1))))||
                                ((C==0)&&(((P[x-a][y][3]==1)||(P[x-a][y][3]==3))&&((P[x][y][3]==2)||(P[x][y][3]==3))))))
      {
         E[x-a][y]+=EP;
         E[x][y]+=EP;
         if(C==1)
         {
            if(P[x-a][y][3]==0) P[x-a][y][3]=1;
            else if(P[x-a][y][3]==2) P[x-a][y][3]=3;
            if(P[x][y][3]==0) P[x][y][3]=2;
            else if(P[x][y][3]==1) P[x][y][3]=3;
         }
         if(C==0)
         {
            if(P[x-a][y][3]==1) P[x-a][y][3]=0;
            else if(P[x-a][y][3]==3) P[x-a][y][3]=2;
            if(P[x][y][3]==2) P[x][y][3]=0;
            else if(P[x][y][3]==3) P[x][y][3]=1;

            R=1;
         }
         *DESis-=2*eA;
      }

      if((x+a<0)||(x+a>Lx-1)){}                                                    //Comprobacion de fronteras
      else if((P[x+a][y][0]==1)&&(P[x+a][y][1]==P[x][y][1])&&(P[x+a][y][2]==1)&&(  //Enlace en direccion contraria
      (((R==2)||(R==3))&&(C==1)&&(((P[x][y][3]==0)||(P[x][y][3]==2))&&((P[x+a][y][3]==0)||(P[x+a][y][3]==1))))||
                        ((C==0)&&(((P[x][y][3]==1)||(P[x][y][3]==3))&&((P[x+a][y][3]==2)||(P[x+a][y][3]==3))))))
      {
         E[x][y]+=EP;
         E[x+a][y]+=EP;
         if(C==1)
         {
            if(P[x][y][3]==0) P[x][y][3]=1;
            else if(P[x][y][3]==2) P[x][y][3]=3;
            if(P[x+a][y][3]==0) P[x+a][y][3]=2;
            else if(P[x+a][y][3]==1) P[x+a][y][3]=3;
         }
         if(C==0)
         {
            if(P[x][y][3]==1) P[x][y][3]=0;
            else if(P[x][y][3]==3) P[x][y][3]=2;
            if(P[x+a][y][3]==2) P[x+a][y][3]=0;
            else if(P[x+a][y][3]==3) P[x+a][y][3]=1;

            if(R==1) R=3;
            else R=2;
         }
         *DESis-=2*eA;
      }
   }
   return R;
}

//Enlaces Actina-CLinker
double EnAL(int Lx,int Ly,double eL,int x,int y,int*** P,double** E,int C)
{
   double EG=0,DESis=0;

   if(C==1) EG=-eL;
   else if(C==0) EG=+eL;

//Enlaces con C-Linkers
   if((P[x][y][1]==1)||(P[x][y][1]==3))                            //Si liga horizontalmente:
   {
      if((x==0)||(x==1)){}
      else if((P[x-1][y][0]==2)&&(P[x-1][y][1]==1)&&(              //Hay una CL en la direccion adecuada?
            ((C==0)&&((P[x-1][y][3]==1)||(P[x-1][y][3]==5)))||     //Hay un enlace que destruir?
            ((C==1)&&((P[x-1][y][3]==0)||(P[x-1][y][3]==3)))))     //Hay posibilidad de establecer un enlace?
      {
         if((P[x-1][y][3]==3)&&(P[x][y][1]!=P[x-2][y][1]))         //Si ya hay un enlace solo se puede unir con una antiparalela
         {
            E[x][y]+=EG;
            E[x-1][y]+=EG;
            P[x-1][y][3]=5;
            DESis-=2*eL;
         }
         else if(P[x-1][y][3]!=3)
         {
            E[x][y]+=EG;
            E[x-1][y]+=EG;
            if(P[x-1][y][3]==1) P[x-1][y][3]=0;
            else if(P[x-1][y][3]==5) P[x-1][y][3]=3;
            else if(P[x-1][y][3]==0) P[x-1][y][3]=1;
            DESis-=2*eL;
         }
      }

      if((x==Lx-1)||(x==Lx-2)){}
      else if((P[x+1][y][0]==2)&&(P[x+1][y][1]==1)&&(
      ((C==0)&&((P[x+1][y][3]==3)||(P[x+1][y][3]==5)))||
      ((C==1)&&((P[x+1][y][3]==0)||(P[x+1][y][3]==1)))))
      {
         if((P[x+1][y][3]==1)&&(P[x][y][1]!=P[x+2][y][1]))
         {
            E[x][y]+=EG;
            E[x+1][y]+=EG;
            P[x+1][y][3]=5;
            DESis-=2*eL;
         }
         else if(P[x+1][y][3]!=1)
         {
            E[x][y]+=EG;
            E[x+1][y]+=EG;
            if(P[x+1][y][3]==3) P[x+1][y][3]=0;
            else if(P[x+1][y][3]==5) P[x+1][y][3]=1;
            else if(P[x+1][y][3]==0) P[x+1][y][3]=3;
            DESis-=2*eL;
         }
      }
   }
   else if((P[x][y][1]==2)||(P[x][y][1]==4))                       //Si liga verticalmente:
   {
      if((y==0)||(y==1)){}
      else if((P[x][y-1][0]==2)&&(P[x][y-1][1]==2)&&(
      ((C==0)&&((P[x][y-1][3]==2)||(P[x][y-1][3]==6)))||
      ((C==1)&&((P[x][y-1][3]==0)||(P[x][y-1][3]==4)))))
      {
         if((P[x][y-1][3]==4)&&(P[x][y][1]!=P[x][y-2][1]))
         {
            E[x][y]+=EG;
            E[x][y-1]+=EG;
            P[x][y-1][3]=6;
            DESis-=2*eL;
         }
         else if(P[x][y-1][3]!=4)
         {
            E[x][y]+=EG;
            E[x][y-1]+=EG;
            if(P[x][y-1][3]==2) P[x][y-1][3]=0;
            else if(P[x][y-1][3]==6) P[x][y-1][3]=4;
            else if(P[x][y-1][3]==0) P[x][y-1][3]=2;
            DESis-=2*eL;
         }
      }

      if((y==Ly-1)||(y==Ly-2)){}
      else if((P[x][y+1][0]==2)&&(P[x][y+1][1]==2)&&(
      ((C==0)&&((P[x][y+1][3]==4)||(P[x][y+1][3]==6)))||
      ((C==1)&&((P[x][y+1][3]==0)||(P[x][y+1][3]==2)))))
      {
         if((P[x][y+1][3]==2)&&(P[x][y][1]!=P[x][y+2][1]))
         {
            E[x][y]+=EG;
            E[x][y+1]+=EG;
            P[x][y+1][3]=6;
            DESis-=2*eL;
         }
         else if(P[x][y+1][3]!=2)
         {
            E[x][y]+=EG;
            E[x][y+1]+=EG;
            if(P[x][y+1][3]==4) P[x][y+1][3]=0;
            else if(P[x][y+1][3]==6) P[x][y+1][3]=2;
            else if(P[x][y+1][3]==0) P[x][y+1][3]=4;
            DESis-=2*eL;
         }
      }
   }
   return DESis;
}

//Enlaces Actina-Myosina
double EnAM(int Lx,int Ly,double eM,int x,int y,int*** P,double** E,int** FM,int C)
{
   double EZ=0,DESis=0;

   if(C==1) EZ=-eM;
   else if(C==0) EZ=+eM;

   if((P[x][y][1]==1)||(P[x][y][1]==3))                                             //Direccion X
   {
      if((x==0)||(x==1)){}
      else if((P[x-1][y][0]==3)&&(P[x-1][y][1]==1)&&(P[x-1][y][3]!=C)&&                //Comprobacion estado Myosina
             (P[x-2][y][0]==1)&&(                                                      //Comprobacion otra actina
             ((P[x][y][1]==1)&&(P[x-2][y][1]==3))||                                    //Combinaciones de direciones
             ((P[x][y][1]==3)&&(P[x-2][y][1]==1))))
      {
         E[x][y]+=EZ;
         E[x-1][y]+=EZ;
         P[x-1][y][3]=C;
         E[x-2][y]+=EZ;

         if(C==0) P[x-1][y][2]=0;
         else if(C==1) P[x-1][y][2]=FM[P[x-1][y][4]][2];

         DESis-=3*eM;
      }

      if((x==Lx-1)||(x==Lx-2)){}
      else if((P[x+1][y][0]==3)&&(P[x+1][y][1]==1)&&(P[x+1][y][3]!=C)&&                //Comprobacion estado Myosina
             (P[x+2][y][0]==1)&&(                                                      //Comprobacion otra actina
             ((P[x][y][1]==1)&&(P[x+2][y][1]==3))||                                    //Combinaciones de direciones
             ((P[x][y][1]==3)&&(P[x+2][y][1]==1))))
      {
         E[x][y]+=EZ;
         E[x+1][y]+=EZ;
         P[x+1][y][3]=C;
         E[x+2][y]+=EZ;

         if(C==0) P[x+1][y][2]=0;
         else if(C==1) P[x+1][y][2]=FM[P[x+1][y][4]][2];

         DESis-=3*eM;
      }
   }
   else if((P[x][y][1]==2)||(P[x][y][1]==4))                                                 //Direccion Y
   {
      if((y==0)||(y==1)){}
      else if((P[x][y-1][0]==3)&&(P[x][y-1][1]==2)&&(P[x][y-1][3]!=C)&&                //Comprobacion estado Myosina
               (P[x][y-2][0]==1)&&(                                                      //Comprobacion otra actina
               ((P[x][y][1]==2)&&(P[x][y-2][1]==4))||                                    //Combinaciones de direciones
               ((P[x][y][1]==4)&&(P[x][y-2][1]==2))))
      {
         E[x][y]+=EZ;
         E[x][y-1]+=EZ;
         P[x][y-1][3]=C;
         E[x][y-2]+=EZ;

         if(C==0) P[x][y-1][2]=0;
         else if(C==1) P[x][y-1][2]=FM[P[x][y-1][4]][2];

         DESis-=3*eM;
      }

      if((y==Ly-1)||(y==Ly-2)){}
      else if((P[x][y+1][0]==3)&&(P[x][y+1][1]==2)&&(P[x][y+1][3]!=C)&&                //Comprobacion estado Myosina
               (P[x][y+2][0]==1)&&(                                                      //Comprobacion otra actina
               ((P[x][y][1]==2)&&(P[x][y+2][1]==4))||                                    //Combinaciones de direciones
               ((P[x][y][1]==4)&&(P[x][y+2][1]==2))))
      {
         E[x][y]+=EZ;
         E[x][y+1]+=EZ;
         P[x][y+1][3]=C;
         E[x][y+2]+=EZ;

         if(C==0) P[x][y+1][2]=0;
         else if(C==1) P[x][y+1][2]=FM[P[x][y+1][4]][2];

         DESis-=3*eM;
      }
   }
   return DESis;
}

//Enlaces CLinker-Actina
double EnLA(int Lx,int Ly,double eL,int x,int y,int*** P,double** E,int C,int R,double* DESis)
{
   double EG=0;

   if(C==1) EG=-eL;
   else if(C==0) EG=+eL;

//Si la particula es una C-Linker:
   if(P[x][y][1]==1)                               // Que enlaza en el eje X
   {
      if((x==0)||(x==Lx-1)){}
      else
      {
         if((P[x-1][y][0]==1)&&((P[x-1][y][1]==1)||(P[x-1][y][1]==3))&&(
               ((C==0)&&((P[x][y][3]==3)||(P[x][y][3]==5)))||
               ((C==1)&&((P[x][y][3]==0)||(P[x][y][3]==1))&&((R==1)||(R==3)||(R==4)))))
         {
            if((P[x][y][3]==1)&&(P[x-1][y][1]!=P[x+1][y][1]))
            {
               E[x][y]+=EG;
               E[x-1][y]+=EG;
               P[x][y][3]=5;
               *DESis-=2*eL;
            }
            else if(P[x][y][3]!=1)
            {
               E[x][y]+=EG;
               E[x-1][y]+=EG;
               if(P[x][y][3]==3)
               {
                  P[x][y][3]=0;
                  R=1;
               }
               else if(P[x][y][3]==5)
               {
                  P[x][y][3]=1;
                  R=1;
               }
               else if(P[x][y][3]==0) P[x][y][3]=3;
               *DESis-=2*eL;
            }
         }
         if((P[x+1][y][0]==1)&&((P[x+1][y][1]==1)||(P[x+1][y][1]==3))&&(
               ((C==0)&&((P[x][y][3]==1)||(P[x][y][3]==5)))||
               ((C==1)&&((P[x][y][3]==0)||(P[x][y][3]==3))&&((R==2)||(R==3)||(R==4)))))
         {
            if((P[x][y][3]==3)&&(P[x-1][y][1]!=P[x+1][y][1]))
            {
               E[x][y]+=EG;
               E[x+1][y]+=EG;
               P[x][y][3]=5;
               *DESis-=2*eL;
            }
            else if(P[x][y][3]!=3)
            {
               E[x][y]+=EG;
               E[x+1][y]+=EG;
               if(P[x][y][3]==1)
               {
                  P[x][y][3]=0;
                  if(R==1) R=3;
                  else R=2;
               }
               else if(P[x][y][3]==5)
               {
                  P[x][y][3]=3;
                  if(R==1) R=3;
                  else R=2;
               }
               else if(P[x][y][3]==0) P[x][y][3]=1;
               *DESis-=2*eL;
            }
         }
      }
   }
   else if(P[x][y][1]==2)                               // Que enlaza en el eje Y
   {
      if((y==0)||(y==Ly-1)){}
      else
      {
         if((P[x][y-1][0]==1)&&((P[x][y-1][1]==2)||(P[x][y-1][1]==4))&&(
            ((C==0)&&((P[x][y][3]==4)||(P[x][y][3]==6)))||
            ((C==1)&&((P[x][y][3]==0)||(P[x][y][3]==2))&&((R==1)||(R==3)||(R==4)))))
         {
            if((P[x][y][3]==2)&&(P[x][y-1][1]!=P[x][y+1][1]))
            {
               E[x][y]+=EG;
               E[x][y-1]+=EG;
               P[x][y][3]=6;
               *DESis-=2*eL;
            }
            else if(P[x][y][3]!=2)
            {
               E[x][y]+=EG;
               E[x][y-1]+=EG;
               if(P[x][y][3]==4)
               {
                  P[x][y][3]=0;
                  R=1;
               }
               else if(P[x][y][3]==6)
               {
                  P[x][y][3]=2;
                  R=1;
               }
               else if(P[x][y][3]==0) P[x][y][3]=4;
               *DESis-=2*eL;
            }
         }
        if((P[x][y+1][0]==1)&&((P[x][y+1][1]==2)||(P[x][y+1][1]==4))&&(
               ((C==0)&&((P[x][y][3]==2)||(P[x][y][3]==6)))||
               ((C==1)&&((P[x][y][3]==0)||(P[x][y][3]==4))&&((R==2)||(R==3)||(R==4)))))
         {
            if((P[x][y][3]==4)&&(P[x][y-1][1]!=P[x][y+1][1]))
            {
               E[x][y]+=EG;
               E[x][y+1]+=EG;
               P[x][y][3]=6;
               *DESis-=2*eL;
            }
            else if(P[x][y][3]!=4)
            {
               E[x][y]+=EG;
               E[x][y+1]+=EG;
               if(P[x][y][3]==2)
               {
                  P[x][y][3]=0;
                  if(R==1) R=3;
                  else R=2;
               }
               else if(P[x][y][3]==6)
               {
                  P[x][y][3]=4;
                  R=1;
               }
               else if(P[x][y][3]==0) P[x][y][3]=2;
               *DESis-=2*eL;
            }
         }
      }
   }
   return R;
}

//Enlaces Myosina-Actina
double EnMA(int Lx,int Ly,double eM,int x,int y,int*** P,double** E,int** FM,int C)
{
   double EZ=0,DESis=0;

   if(C==1) EZ=-eM;
   else if(C==0) EZ=+eM;

   if(P[x][y][3]!=C)
   {
      if(P[x][y][1]==1)
      {
         if((x==0)||(x==Lx-1)){}
         else if((P[x-1][y][0]==1)&&(P[x+1][y][0]==1)&&(                         //Comprobacion existencia actina adyacente a la myosina
                ((P[x-1][y][1]==1)&&(P[x+1][y][1]==3))||                         //Comprobacion para una combinacion
                ((P[x-1][y][1]==3)&&(P[x+1][y][1]==1))))                         //La otra combinacion
         {
            E[x-1][y]+=EZ;
            E[x][y]+=EZ;
            P[x][y][3]=C;
            E[x+1][y]+=EZ;

            if(C==0) P[x][y][2]=0;
            else if(C==1) P[x][y][2]=FM[P[x][y][4]][2];

            DESis-=3*eM;
         }
      }
      if(P[x][y][1]==2)
      {
         if((y==0)||(y==Ly-1)){}
         else if((P[x][y-1][0]==1)&&(P[x][y+1][0]==1)&&(                         //Comprobacion existencia actina adyacente a la myosina
                ((P[x][y-1][1]==2)&&(P[x][y+1][1]==4))||                         //Comprobacion para una combinacion
                ((P[x][y-1][1]==4)&&(P[x][y+1][1]==2))))                         //La otra combinacion
         {
            E[x][y-1]+=EZ;
            E[x][y]+=EZ;
            P[x][y][3]=C;
            E[x][y+1]+=EZ;

            if(C==0) P[x][y][2]=0;
            else if(C==1) P[x][y][2]=FM[P[x][y][4]][2];

            DESis-=3*eM;
         }
      }
   }
   return DESis;
}

//Envejecimiento del filamento de actina
void HATP(int Lx,int Ly,double eA,double eL,double eM,int*** P,double** E,int** FM,int x,int y,long* HA)
{
   double DESis=0;

   if(P[x][y][1]==1)
   {
      (*HA)++;
      En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,0,0,&DESis);
      P[x][y+1][2]=0;
      En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,1,0,&DESis);
   }
   else if(P[x][y][1]==2)
   {
      (*HA)++;
      En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,0,0,&DESis);
      P[x+1][y][2]=0;
      En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,1,0,&DESis);
   }
   else if(P[x][y][1]==3)
   {
      (*HA)++;
      En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,0,0,&DESis);
      P[x][y-1][2]=0;
      En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,1,0,&DESis);
   }
   else if(P[x][y][1]==4)
   {
      (*HA)++;
      En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,0,0,&DESis);
      P[x-1][y][2]=0;
      En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,1,0,&DESis);
   }
}

//Entrada o salida de particulas mediante un reservorio
void ResV(int Lx,int Ly,double eA,double eL,double eM,double mu,int*** P,double** E,int** FM,int x,int y,long* RI,long* RO,int* N,int* NX,int tpro,int OldT,int OldS)
{
   double p=0;
   int R=0;
   double DESis=0;

   if(P[x][y][0]==0)                 //Si no hay particula se prueba a introducir una
   {
      if(*NX!=0)                     //Entra solo si hay proteinas en el systema
      {
         P[x][y][0]=tpro;
         if(tpro==1) P[x][y][1]=randomInt(4)+1;        //Eleccion aleatoria de la direccion
         else P[x][y][1]=randomInt(2)+1;
         if(tpro==1) P[x][y][2]=1;                            //La actina entra con ATP
         else P[x][y][2]=0;

         En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,1,4,&DESis);   //Calculo de la energia al entrar

         p=ramdomFloat();

         if(p<(exp(-DESis+mu)))                      //Eleccion de si la proteina entra o no
         {
            *RI=*RI+1;
            (*N)++;
            (*NX)--;
         }
         else                                        //Si no entra se recupera la configuracion inicial
         {
            En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,0,0,&DESis);

            P[x][y][0]=0;
            P[x][y][1]=0;
            P[x][y][2]=0;
         }
      }
   }
   else if(P[x][y][0]!=0)                           //Si hay una particula se prueba a sacarla
   {
      R=En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,0,0,&DESis);

      p=ramdomFloat();

      if(p<(exp(-mu+DESis)))       //Si se saca
      {
         (*RO)++;
         (*N)--;
         (*NX)++;
         P[x][y][0]=0;
         P[x][y][1]=0;
         P[x][y][2]=0;
      }
      else
      {
         En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,1,R,&DESis);
      }
   }
}

//Autoregulacion del potencial quimico
inline double RegMu(double muE,int N,int NX,int D)
{
   double mu=0;

   if(NX!=0)
   {
      if(N!=0) mu=muE-log((double)(D+1)*N/(NX+N));
      else mu=muE-log((double)(D+1)*0.001/(NX+0.001));
   }
   return mu;
}



//Registro del crecimiento de los filamentos
void FilNew(int Lx,int Ly,int*** P,int x,int y,int** FA,int** FL,int** FM,double* TA,int t,double TBreak)
{
   int nf=0,nl=0,nm=0;

   if((P[x][y][0]==1)&&(P[x][y][3]==2)) //El filamento crece (+Actina en el lado debil)
   {
      int a=0,b=0,c=0,d=0;
                                                      //Generalizacion para acortar codigo
      if((P[x][y][1]==1)||(P[x][y][1]==2)) a=1;          //Direccion positiva
      else if((P[x][y][1]==3)||(P[x][y][1]==4)) a=-1;    //Direccion negativa

      if(((P[x][y][1]==1)||(P[x][y][1]==3))&&(P[x][y-a][3]==3))  //Eje Y
      {
         if(P[x][y-a][4]!=0)                             //Filamento ya existente
         {
            nf=P[x][y-a][4];                             //Reconocimiento
            P[x][y][4]=nf;                               //Addicion al filamento

            FA[nf][3]++;

            FA[0][1]++;

            for(b=-1;b<=1;b+=2)                          //Etiquetaje de CLinkers o Myosinas adyacentes
            {
               if((x+b>-1)&&(x+b<Lx))
               {
                  if((P[x+b][y][0]==2)&&(P[x+b][y][3]==5)&&(P[x+(2*b)][y][4]!=0))  //CLinker
                  {
                     for(nl=1;nl<=FL[0][0];nl++) if(FL[nl][0]==0) break;   //Eleccion etiqueta
                     if(nl>FL[0][0]) FL[0][0]=nl;
                     P[x+b][y][4]=nl;                                      //Etiquetage

                     FL[nl][0]=nf;                                         //Etiquetas de los filamentos adyacentes
                     FL[nl][1]=P[x+(2*b)][y][4];

                     FA[0][2]++;
                  }
                  else if((P[x+b][y][0]==3)&&(P[x+b][y][1]==1)&&(P[x+b][y][3]==1)&&(P[x+(2*b)][y][4]!=0))   //Myosina
                  {
                     for(nm=1;nm<=FM[0][0];nm++) if(FM[nm][0]==0) break;   //Eleccion etiqueta
                     if(nm>FM[0][0]) FM[0][0]=nm;
                     P[x+b][y][4]=nm;                                      //Etiquetage

                     FM[nm][0]=nf;                                         //Etiquetas de los filamentos adyacentes
                     FM[nm][1]=P[x+(2*b)][y][4];
                     FM[nm][2]=P[x+b][y][2];                               //Numero de ATP en la myosina
                     FM[nm][3]=t;                                          //Tiempo de la ultima absorcion

                     FA[0][3]++;
                  }
               }
            }
         }
         else                                                        //Filamento nuevo
         {
            for(nf=1;nf<=FA[0][0];nf++) if(FA[nf][3]==0) break;      //Eleccion etiqueta filamento
            if(nf>FA[0][0]) FA[0][0]=nf;
            TA[0]=FA[0][0];

            FA[nf][0]=x;                                             //Creacion filamento
            FA[nf][1]=y-(2*a);
            FA[nf][2]=P[x][y][1];
            FA[nf][3]=3;
            TA[nf]=0;
            FA[0][1]+=3;

            for(c=0;c<=2;c++)                                        //Etiquetaje
            {
               d=y-(c*a);
               P[x][d][4]=nf;                                           //Actinas

               for(b=-1;b<=1;b+=2)
               {
                  if((x+b>-1)&&(x+b<Lx))
                  {
                     if((P[x+b][d][0]==2)&&(P[x+b][d][3]==5)&&(P[x+(2*b)][d][4]!=0)) //CLinker
                     {
                        for(nl=1;nl<=FL[0][0];nl++) if(FL[nl][0]==0) break;
                        if(nl>FL[0][0]) FL[0][0]=nl;
                        P[x+b][d][4]=nl;

                        FL[nl][0]=nf;
                        FL[nl][1]=P[x+(2*b)][d][4];

                        FA[0][2]++;
                     }
                     else if((P[x+b][d][0]==3)&&(P[x+b][d][1]==1)&&(P[x+b][d][3]==1)&&(P[x+(2*b)][d][4]!=0)) //Myosin
                     {
                        for(nm=1;nm<=FM[0][0];nm++) if(FM[nm][0]==0) break;
                        if(nm>FM[0][0]) FM[0][0]=nm;
                        P[x+b][d][4]=nm;

                        FM[nm][0]=nf;
                        FM[nm][1]=P[x+(2*b)][d][4];
                        FM[nm][2]=P[x+b][d][2];
                        FM[nm][3]=t;

                        FA[0][3]++;
                     }
                  }
               }
            }
         }
      }
      else if(((P[x][y][1]==2)||(P[x][y][1]==4))&&(P[x-a][y][3]==3))   //Eje X
      {
         if(P[x-a][y][4]!=0)
         {
            nf=P[x-a][y][4];
            P[x][y][4]=nf;

            FA[nf][3]++;

            FA[0][1]++;

            for(b=-1;b<=1;b+=2)
            {
               if((y+b>-1)&&(y+b<Ly))
               {
                  if((P[x][y+b][0]==2)&&(P[x][y+b][3]==6)&&(P[x][y+(2*b)][4]!=0))
                  {
                     for(nl=1;nl<=FL[0][0];nl++) if(FL[nl][0]==0) break;
                     if(nl>FL[0][0]) FL[0][0]=nl;
                     P[x][y+b][4]=nl;

                     FL[nl][0]=nf;
                     FL[nl][1]=P[x][y+(2*b)][4];

                     FA[0][2]++;
                  }
                  else if((P[x][y+b][0]==3)&&(P[x][y+b][1]==2)&&(P[x][y+b][3]==1)&&(P[x][y+(2*b)][4]!=0))
                  {
                     for(nm=1;nm<=FM[0][0];nm++) if(FM[nm][0]==0) break;
                     if(nm>FM[0][0]) FM[0][0]=nm;
                     P[x][y+b][4]=nm;

                     FM[nm][0]=nf;
                     FM[nm][1]=P[x][y+(2*b)][4];
                     FM[nm][2]=P[x][y+b][2];
                     FM[nm][3]=t;

                     FA[0][3]++;
                  }
               }
            }
         }
         else
         {
            for(nf=1;nf<=FA[0][0];nf++) if(FA[nf][3]==0) break;
            if(nf>FA[0][0]) FA[0][0]=nf;
            TA[0]=FA[0][0];

            FA[nf][0]=x-(2*a);
            FA[nf][1]=y;
            FA[nf][2]=P[x][y][1];
            FA[nf][3]=3;
            TA[nf]=0;
            FA[0][1]+=3;

            for(c=0;c<=2;c++)
            {
               d=x-(c*a);
               P[d][y][4]=nf;

               for(b=-1;b<=1;b+=2)
               {
                  if((y+b>-1)&&(y+b<Ly))
                  {
                     if((P[d][y+b][0]==2)&&(P[d][y+b][3]==6)&&(P[d][y+(2*b)][4]!=0))
                     {
                        for(nl=1;nl<=FL[0][0];nl++) if(FL[nl][0]==0) break;
                        if(nl>FL[0][0]) FL[0][0]=nl;
                        P[d][y+b][4]=nl;

                        FL[nl][0]=nf;
                        FL[nl][1]=P[d][y+(2*b)][4];

                        FA[0][2]++;
                     }
                     else if((P[d][y+b][0]==3)&&(P[d][y+b][1]==2)&&(P[d][y+b][3]==1)&&(P[d][y+(2*b)][4]!=0))
                     {
                        for(nm=1;nm<=FM[0][0];nm++) if(FM[nm][0]==0) break;
                        if(nm>FM[0][0]) FM[0][0]=nm;
                        P[d][y+b][4]=nm;

                        FM[nm][0]=nf;
                        FM[nm][1]=P[d][y+(2*b)][4];
                        FM[nm][2]=P[d][y+b][2];
                        FM[nm][3]=t;

                        FA[0][3]++;
                     }
                  }
               }
            }
         }
      }
   }
   else if(P[x][y][0]==2) //Se enlaza una CL a un filamento
   {
      if((P[x][y][3]==5)&&(P[x-1][y][4]!=0)&&(P[x+1][y][4]!=0)) //Eje X
      {
         for(nl=1;nl<=FL[0][0];nl++) if(FL[nl][0]==0) break;
         if(nl>FL[0][0]) FL[0][0]=nl;
         P[x][y][4]=nl;

         FL[nl][0]=P[x-1][y][4];
         FL[nl][1]=P[x+1][y][4];

         FA[0][2]++;
      }
      else if((P[x][y][3]==6)&&(P[x][y-1][4]!=0)&&(P[x][y+1][4]!=0)) //Eje Y
      {
         for(nl=1;nl<=FL[0][0];nl++) if(FL[nl][0]==0) break;
         if(nl>FL[0][0]) FL[0][0]=nl;
         P[x][y][4]=nl;

         FL[nl][0]=P[x][y-1][4];
         FL[nl][1]=P[x][y+1][4];

         FA[0][2]++;
      }
   }
   else if(P[x][y][0]==3)  //Se enlaza una myosina a un filamento
   {
      if((P[x][y][3]==1)&&(P[x][y][1]==1)&&(P[x-1][y][4]!=0)&&(P[x+1][y][4]!=0)) //Eje X
      {
         for(nm=1;nm<=FM[0][0];nm++) if(FM[nm][0]==0) break;
         if(nm>FM[0][0]) FM[0][0]=nm;
         P[x][y][4]=nm;

         FM[nm][0]=P[x-1][y][4];
         FM[nm][1]=P[x+1][y][4];
         FM[nm][2]=P[x][y][2];
         FM[nm][3]=t;

         FA[0][3]++;
      }
      else if((P[x][y][3]==1)&&(P[x][y][1]==2)&&(P[x][y-1][4]!=0)&&(P[x][y+1][4]!=0)) //Eje Y
      {
         for(nm=1;nm<=FM[0][0];nm++) if(FM[nm][0]==0) break;
         if(nm>FM[0][0]) FM[0][0]=nm;
         P[x][y][4]=nm;

         FM[nm][0]=P[x][y-1][4];
         FM[nm][1]=P[x][y+1][4];
         FM[nm][2]=P[x][y][2];
         FM[nm][3]=t;

         FA[0][3]++;
      }
   }
}

//Registro de la rotura de los filamentos
void FilOld(int Lx,int Ly,int OldT,int OldS,int*** P,int x,int y,int** FA,int** FL,int** FM,double* TA)
{
   int nf=0,nl=0,nm=0;
   int e=0;

   if(OldT==1) {       //Se libera una actina
      int a=0,b=0,c=0,d=0;
      nf=P[x][y][4];

      if((FA[nf][2]==1)||(FA[nf][2]==2)) a=1;
      else if((FA[nf][2]==3)||(FA[nf][2]==4)) a=-1;

      if(OldS==1){    //Se libera una actina el extremo fuerte
         FA[nf][3]--;   //Reducion L
         if(FA[nf][3]>=3) {  //Si se conserva el filamento se borran las etiquetas de esan actina solo
            P[x][y][4]=0;
            FA[0][1]--;
            if((FA[nf][2]==1)||(FA[nf][2]==3)){       // Direccion Y
               FA[nf][1]=y+a;  // ??? cambio inicio porque ???
               for(b=-1;b<=1;b+=2){                      // Laterales
                  if((x+b>-1)&&(x+b<Lx)){
                     if((P[x+b][y][0]==2)&&(P[x+b][y][1]==1)&&(P[x+b][y][4]!=0)) {  //MTYPE = CRL (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                        nl=P[x+b][y][4];                    // Id CRL
                        if(nl==FL[0][0]) FL[0][0]--;        // Total CRL num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta, no se actuliza FA[0][NCRL]
                        P[x+b][y][4]=0;                     // CRL fuera de listas
                        for(e=0;e<=1;e++) FL[nl][e]=0;      // CRL suelto : no se elimina ???
                        FA[0][2]--;                         // Total CRL FA[0][NCRL]--
                        }
                     else if((P[x+b][y][0]==3)&&(P[x+b][y][1]==1)&&(P[x+b][y][4]!=0)) { //MTYPE = MYO (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                        nm=P[x+b][y][4];                    // Id MYO
                        if(nm==FM[0][0]) FM[0][0]--;        // Total MYO num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta
                        P[x+b][y][4]=0;                     // MYO fuera de listas
                        P[x+b][y][2]=0;                     // MYO sin ATP
                        for(e=0;e<=3;e++) FM[nm][e]=0;      // MYO suelto : no se elimina ???
                        FA[0][3]--;                         // Total MYO FA[0][NMYO]--
                        }
                     }
                  }
               }
            else if((FA[nf][2]==2)||(FA[nf][2]==4)) { // Direccion X
               FA[nf][0]=x+a;  // ??? cambio inicio porque ???
               for(b=-1;b<=1;b+=2) {  // Laterales
                  if((y+b>-1)&&(y+b<Ly)) {
                     if((P[x][y+b][0]==2)&&(P[x][y+b][1]==2)&&(P[x][y+b][4]!=0)){ //MTYPE = CRL (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                        nl=P[x][y+b][4];                // Id CRL
                        if(nl==FL[0][0]) FL[0][0]--;    // Total CRL num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta, no se actuliza FA[0][NCRL]
                        P[x][y+b][4]=0;                 // CRL fuera de listas
                        for(e=0;e<=1;e++) FL[nl][e]=0;  // CRL suelto : no se elimina ???
                        FA[0][2]--;                     // Total CRL FA[0][NCRL]--
                        }
                     else if((P[x][y+b][0]==3)&&(P[x][y+b][1]==2)&&(P[x][y+b][4]!=0)){ //MTYPE = MYO (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                        nm=P[x][y+b][4];               // Id MYO
                        if(nm==FM[0][0]) FM[0][0]--;   // Total MYO num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta
                        P[x][y+b][4]=0;                // MYO fuera de listas
                        P[x][y+b][2]=0;                // MYO sin ATP
                        for(e=0;e<=3;e++) FM[nm][e]=0; // MYO suelto : no se elimina ???
                        FA[0][3]--;                    // Total MYO FA[0][NMYO]--
                        }
                     }
                  }
               }
            }
         else {     //Si se elimina el filamento se borran todas las actinas que lo formaban
            if((FA[nf][2]==1)||(FA[nf][2]==3)) {      // Direccion Y
               for(c=0;c<=2;c++) {
                  d=y+(c*a);
                  P[x][d][4]=0;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((x+b>-1)&&(x+b<Lx))
                     {
                        if((P[x+b][d][0]==2)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nl=P[x+b][d][4];
                           if(nl==FL[0][0]) FL[0][0]--;
                           P[x+b][d][4]=0;

                           for(e=0;e<=1;e++) FL[nl][e]=0;
                           FA[0][2]--;
                        }
                        else if((P[x+b][d][0]==3)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nm=P[x+b][d][4];
                           if(nm==FM[0][0]) FM[0][0]--;
                           P[x+b][d][4]=0;
                           P[x+b][d][2]=0;

                           for(e=0;e<=3;e++) FM[nm][e]=0;
                           FA[0][3]--;
                        }
                     }
                  }
               }
            }
            else if((FA[nf][2]==2)||(FA[nf][2]==4)) { // Direccion X
               for(c=0;c<=2;c++) {
                  d=x+(c*a);
                  P[d][y][4]=0;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((y+b>-1)&&(y+b<Ly))
                     {
                        if((P[d][y+b][0]==2)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nl=P[d][y+b][4];
                           if(nl==FL[0][0]) FL[0][0]--;
                           P[d][y+b][4]=0;

                           for(e=0;e<=1;e++) FL[nl][e]=0;
                           FA[0][2]--;
                        }
                        else if((P[d][y+b][0]==3)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nm=P[d][y+b][4];
                           if(nm==FM[0][0]) FM[0][0]--;
                           P[d][y+b][4]=0;
                           P[d][y+b][2]=0;

                           for(e=0;e<=3;e++) FM[nm][e]=0;
                           FA[0][3]--;
                        }
                     }
                  }
               }
            }

            for(e=0;e<=3;e++) FA[nf][e]=0;
            if(nf==FA[0][0]) FA[0][0]--;
            TA[0]=FA[0][0];
            TA[nf]=0;
            FA[0][1]-=3;
         }
      }
      else if(OldS==2)                             //Se libera una actina del extremo debil
      {
         FA[nf][3]--;

         if(FA[nf][3]>=3)
         {
            P[x][y][4]=0;
            FA[0][1]--;

            if((FA[nf][2]==1)||(FA[nf][2]==3))
            {
               for(b=-1;b<=1;b+=2)
               {
                  if((x+b>-1)&&(x+b<Lx))
                  {
                     if((P[x+b][y][0]==2)&&(P[x+b][y][1]==1)&&(P[x+b][y][4]!=0))
                     {
                        nl=P[x+b][y][4];
                        if(nl==FL[0][0]) FL[0][0]--;
                        P[x+b][y][4]=0;

                        for(e=0;e<=1;e++) FL[nl][e]=0;
                        FA[0][2]--;
                     }
                     else if((P[x+b][y][0]==3)&&(P[x+b][y][1]==1)&&(P[x+b][y][4]!=0))
                     {
                        nm=P[x+b][y][4];
                        if(nm==FM[0][0]) FM[0][0]--;
                        P[x+b][y][4]=0;
                        P[x+b][y][2]=0;

                        for(e=0;e<=3;e++) FM[nm][e]=0;
                        FA[0][3]--;
                     }
                  }
               }
            }
            else if((FA[nf][2]==2)||(FA[nf][2]==4))
            {
               for(b=-1;b<=1;b+=2)
               {
                  if((y+b>-1)&&(y+b<Ly))
                  {
                     if((P[x][y+b][0]==2)&&(P[x][y+b][1]==2)&&(P[x][y+b][4]!=0))
                     {
                        nl=P[x][y+b][4];
                        if(nl==FL[0][0]) FL[0][0]--;
                        P[x][y+b][4]=0;

                        for(e=0;e<=1;e++) FL[nl][e]=0;
                        FA[0][2]--;
                     }
                     else if((P[x][y+b][0]==3)&&(P[x][y+b][1]==2)&&(P[x][y+b][4]!=0))
                     {
                        nm=P[x][y+b][4];
                        if(nm==FM[0][0]) FM[0][0]--;
                        P[x][y+b][4]=0;
                        P[x][y+b][2]=0;

                        for(e=0;e<=3;e++) FM[nm][e]=0;
                        FA[0][3]--;
                     }
                  }
               }
            }
         }
         else
         {
            if((FA[nf][2]==1)||(FA[nf][2]==3))
            {
               for(c=0;c<=2;c++)
               {
                  d=y-(c*a);
                  P[x][d][4]=0;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((x+b>-1)&&(x+b<Lx))
                     {
                        if((P[x+b][d][0]==2)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nl=P[x+b][d][4];
                           if(nl==FL[0][0]) FL[0][0]--;
                           P[x+b][d][4]=0;

                           for(e=0;e<=1;e++) FL[nl][e]=0;
                           FA[0][2]--;
                        }
                        else if((P[x+b][d][0]==3)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nm=P[x+b][d][4];
                           if(nm==FM[0][0]) FM[0][0]--;
                           P[x+b][d][4]=0;
                           P[x+b][d][2]=0;

                           for(e=0;e<=3;e++) FM[nm][e]=0;
                           FA[0][3]--;
                        }
                     }
                  }
               }
            }
            else if((FA[nf][2]==2)||(FA[nf][2]==4))
            {
               for(c=0;c<=2;c++)
               {
                  d=x-(c*a);
                  P[d][y][4]=0;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((y+b>-1)&&(y+b<Ly))
                     {
                        if((P[d][y+b][0]==2)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nl=P[d][y+b][4];
                           if(nl==FL[0][0]) FL[0][0]--;
                           P[d][y+b][4]=0;

                           for(e=0;e<=1;e++) FL[nl][e]=0;
                           FA[0][2]--;
                        }
                        else if((P[d][y+b][0]==3)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nm=P[d][y+b][4];
                           if(nm==FM[0][0]) FM[0][0]--;
                           P[d][y+b][4]=0;
                           P[d][y+b][2]=0;

                           for(e=0;e<=3;e++) FM[nm][e]=0;
                           FA[0][3]--;
                        }
                     }
                  }
               }
            }

            for(e=0;e<=3;e++) FA[nf][e]=0;
            if(nf==FA[0][0]) FA[0][0]--;
            TA[0]=FA[0][0];
            TA[nf]=0;
            FA[0][1]-=3;
         }
      }
      else if(OldS==3)                          //Se libera una actina central
      {
         int MemL=FA[nf][3];

         FA[0][1]-=MemL;
         TA[nf]=0;                              //Se elimina la tension
         P[x][y][4]=0;

         if((FA[nf][2]==1)||(FA[nf][2]==3))     //Eje Y
         {
            for(b=-1;b<=1;b+=2)                 //Borrado CL o Myo adyacentes
            {
               if((x+b>-1)&&(x+b<Lx))
               {
                  if((P[x+b][y][0]==2)&&(P[x+b][y][1]==1)&&(P[x+b][y][4]!=0))
                  {
                     nl=P[x+b][y][4];
                     if(nl==FL[0][0]) FL[0][0]--;
                     P[x+b][y][4]=0;

                     for(e=0;e<=1;e++) FL[nl][e]=0;
                     FA[0][2]--;
                  }
                  else if((P[x+b][y][0]==3)&&(P[x+b][y][1]==1)&&(P[x+b][y][4]!=0))
                  {
                     nm=P[x+b][y][4];
                     if(nm==FM[0][0]) FM[0][0]--;
                     P[x+b][y][4]=0;
                     P[x+b][y][2]=0;

                     for(e=0;e<=3;e++) FM[nm][e]=0;
                     FA[0][3]--;
                  }
               }
            }

            int MemY=FA[nf][1];

            if(abs(y-MemY)>=3)          //Segmento debil es filamento
            {
               FA[nf][3]=abs(y-MemY);   //Se conserva la etiqueta
               FA[0][1]+=FA[nf][3];

               for(c=1;c<=FA[nf][3];c++)  //Se elimina la tension
               {
                  d=y-(a*c);
                  for(b=-1;b<=1;b+=2)
                  {
                     if((x+b>-1)&&(x+b<Lx))
                     {
                        if((P[x+b][d][0]==3)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nm=P[x+b][d][4];
                           P[x+b][d][2]=0;
                           FM[nm][2]=0;
                        }
                     }
                  }
               }
            }
            else  //Si el segmento debil no llega a filamento se eliminan todas las etiquetas y se suprime el filamento
            {
               for(e=0;e<=3;e++) FA[nf][e]=0;
               if(nf==FA[0][0]) FA[0][0]--;
               TA[0]=FA[0][0];

               for(c=1;c<=abs(y-MemY);c++)
               {
                  d=y-(c*a);
                  P[x][d][4]=0;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((x+b>-1)&&(x+b<Lx))
                     {
                        if((P[x+b][d][0]==2)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nl=P[x+b][d][4];
                           if(nl==FL[0][0]) FL[0][0]--;
                           P[x+b][d][4]=0;

                           for(e=0;e<=1;e++) FL[nl][e]=0;
                           FA[0][2]--;
                        }
                        else if((P[x+b][d][0]==3)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nm=P[x+b][d][4];
                           if(nm==FM[0][0]) FM[0][0]--;
                           P[x+b][d][4]=0;
                           P[x+b][d][2]=0;

                           for(e=0;e<=3;e++) FM[nm][e]=0;
                           FA[0][3]--;
                        }
                     }
                  }
               }
            }

            if((MemL-abs(y-MemY)-1)>=3)                        //Si el segmento fuerte es filamento se crea un filamento nuevo
            {
               for(nf=1;nf<=FA[0][0];nf++) if(FA[nf][3]==0) break;
               if(nf>FA[0][0]) FA[0][0]=nf;
               TA[0]=FA[0][0];

               FA[nf][0]=x;
               FA[nf][1]=y+a;
               FA[nf][2]=P[x][y+a][1];
               FA[nf][3]=MemL-abs(y-MemY)-1;
               TA[nf]=0;
               FA[0][1]+=FA[nf][3];

               for(c=1;c<=FA[nf][3];c++)
               {
                  d=y+(c*a);
                  P[x][d][4]=nf;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((x+b>-1)&&(x+b<Lx))
                     {
                        if((P[x+b][d][0]==2)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nl=P[x+b][d][4];
                           FL[nl][0]=nf;
                           FL[nl][1]=P[x+(2*b)][d][4];

                           FA[0][2]++;
                        }
                        else if((P[x+b][d][0]==3)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nm=P[x+b][d][4];
                           FM[nm][0]=nf;
                           FM[nm][1]=P[x+(2*b)][d][4];

                           P[x+b][d][2]=0;
                           FM[nm][2]=0;

                           FA[0][3]++;
                        }
                     }
                  }
               }
            }
            else  //Si el segmento fuerte no llega a filamento se eliminan todas las etiquetas y se suprime el filamento
            {
               for(c=1;c<=(MemL-abs(y-MemY)-1);c++)
               {
                  d=y+(c*a);
                  P[x][d][4]=0;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((x+b>-1)&&(x+b<Lx))
                     {
                        if((P[x+b][d][0]==2)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nl=P[x+b][d][4];
                           if(nl==FL[0][0]) FL[0][0]--;
                           P[x+b][d][4]=0;

                           for(e=0;e<=1;e++) FL[nl][e]=0;
                           FA[0][2]--;
                        }
                        else if((P[x+b][d][0]==3)&&(P[x+b][d][1]==1)&&(P[x+b][d][4]!=0))
                        {
                           nm=P[x+b][d][4];
                           if(nm==FM[0][0]) FM[0][0]--;
                           P[x+b][d][4]=0;
                           P[x+b][d][2]=0;

                           for(e=0;e<=3;e++) FM[nm][e]=0;
                           FA[0][3]--;
                        }
                     }
                  }
               }
            }
         }
         else if((FA[nf][2]==2)||(FA[nf][2]==4)) //Eje X
         {
            for(b=-1;b<=1;b+=2)
            {
               if((y+b>-1)&&(y+b<Ly))
               {
                  if((P[x][y+b][0]==2)&&(P[x][y+b][1]==2)&&(P[x][y+b][4]!=0))
                  {
                     nl=P[x][y+b][4];
                     if(nl==FL[0][0]) FL[0][0]--;
                     P[x][y+b][4]=0;

                     for(e=0;e<=1;e++) FL[nl][e]=0;
                     FA[0][2]--;
                  }
                  else if((P[x][y+b][0]==3)&&(P[x][y+b][1]==2)&&(P[x][y+b][4]!=0))
                  {
                     nm=P[x][y+b][4];
                     if(nm==FM[0][0]) FM[0][0]--;
                     P[x][y+b][4]=0;
                     P[x][y+b][2]=0;

                     for(e=0;e<=3;e++) FM[nm][e]=0;
                     FA[0][3]--;
                  }
               }
            }

            int MemX=FA[nf][0];

            if(abs(x-MemX)>=3)
            {
               FA[nf][3]=abs(x-MemX);
               FA[0][1]+=FA[nf][3];

               for(c=1;c<=FA[nf][3];c++)
               {
                  d=x-(a*c);

                  for(b=-1;b<=1;b+=2)
                  {
                     if((y+b>-1)&&(y+b<Ly))
                     {
                        if((P[d][y+b][0]==3)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nm=P[d][y+b][4];
                           P[x][y+b][2]=0;
                           FM[nm][2]=0;
                        }
                     }
                  }
               }
            }
            else
            {
               for(e=0;e<=3;e++) FA[nf][e]=0;
               if(nf==FA[0][0]) FA[0][0]--;
               TA[0]=FA[0][0];

               for(c=1;c<=abs(x-MemX);c++)
               {
                  d=x-(c*a);
                  P[d][y][4]=0;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((y+b>-1)&&(y+b<Ly))
                     {
                        if((P[d][y+b][0]==2)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nl=P[d][y+b][4];
                           if(nl==FL[0][0]) FL[0][0]--;
                           P[d][y+b][4]=0;

                           for(e=0;e<=1;e++) FL[nl][e]=0;
                           FA[0][2]--;
                        }
                        else if((P[d][y+b][0]==3)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nm=P[d][y+b][4];
                           if(nm==FM[0][0]) FM[0][0]--;
                           P[d][y+b][4]=0;
                           P[d][y+b][2]=0;

                           for(e=0;e<=3;e++) FM[nm][e]=0;
                           FA[0][3]--;
                        }
                     }
                  }
               }
            }

            if((MemL-abs(x-MemX)-1)>=3)
            {
               for(nf=1;nf<=FA[0][0];nf++) if(FA[nf][3]==0) break;
               if(nf>FA[0][0]) FA[0][0]=nf;
               TA[0]=FA[0][0];

               FA[nf][0]=x+a;
               FA[nf][1]=y;
               FA[nf][2]=P[x+a][y][1];
               FA[nf][3]=MemL-abs(x-MemX)-1;
               TA[nf]=0;
               FA[0][1]+=FA[nf][3];

               for(c=1;c<=FA[nf][3];c++)
               {
                  d=x+(c*a);
                  P[d][y][4]=nf;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((y+b>-1)&&(y+b<Ly))
                     {
                        if((P[d][y+b][0]==2)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nl=P[d][y+b][4];
                           FL[nl][0]=nf;
                           FL[nl][1]=P[d][y+(2*b)][4];

                           FA[0][2]++;
                        }
                        else if((P[d][y+b][0]==3)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nm=P[d][y+b][4];
                           FM[nm][0]=nf;
                           FM[nm][1]=P[d][y+(2*b)][4];

                           P[x][y+b][2]=0;
                           FM[nm][2]=0;

                           FA[0][3]++;
                        }
                     }
                  }
               }
            }
            else
            {
               for(c=1;c<=(MemL-abs(x-MemX)-1);c++)
               {
                  d=x+(c*a);
                  P[d][y][4]=0;

                  for(b=-1;b<=1;b+=2)
                  {
                     if((y+b>-1)&&(y+b<Ly))
                     {
                        if((P[d][y+b][0]==2)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nl=P[d][y+b][4];
                           if(nl==FL[0][0]) FL[0][0]--;
                           P[d][y+b][4]=0;

                           for(e=0;e<=1;e++) FL[nl][e]=0;
                           FA[0][2]--;
                        }
                        else if((P[d][y+b][0]==3)&&(P[d][y+b][1]==2)&&(P[d][y+b][4]!=0))
                        {
                           nm=P[d][y+b][4];
                           if(nm==FM[0][0]) FM[0][0]--;
                           P[d][y+b][4]=0;
                           P[d][y+b][2]=0;

                           for(e=0;e<=3;e++) FM[nm][e]=0;
                           FA[0][3]--;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   else if(OldT==2)  { //Una CL se suelta de un filamento
   
      nl=P[x][y][4];
      if(nl==FL[0][0]) FL[0][0]--;
      P[x][y][4]=0;

      for(e=0;e<=1;e++) FL[nl][e]=0;
      FA[0][2]--;
   }
   else if(OldT==3) {  //Una Myo se suelta de un filamento
      nm=P[x][y][4];
      if(nm==FM[0][0]) FM[0][0]--;
      P[x][y][4]=0;
      P[x][y][2]=0;

      for(e=0;e<=3;e++) FM[nm][e]=0;
      FA[0][3]--;
   }
}

//Registro de vecinos mediante las etiquetas de CL y Myo
inline void Neighbor(int Lx, int Ly,int*** P,int** FA,int** FL,int** FM,int** NH)
{
   int x=0,y=0,n=0,nf1=0,nf2=0;
   int Here=0;

   NH[0][0]=FA[0][0];

   for(x=0;x<Lx;x++)
   {
      for(y=0;y<Ly;y++)
      {
         if(((P[x][y][0]==2)||(P[x][y][0]==3))&&(P[x][y][4]!=0)) //Hay una CL o Myo que une filamentos
         {
            n=P[x][y][4];
            if(P[x][y][0]==2)
            {
               nf1=FL[n][0];
               nf2=FL[n][1];
            }
            else if(P[x][y][0]==3)
            {
               nf1=FM[n][0];
               nf2=FM[n][1];
            }

            //Se registran los filamentos adyacentes evitando repeticiones
            Here=0;
            for(n=1;n<=NH[nf1][0];n++) if(NH[nf1][n]==nf2) Here=1;
            if(Here==0)
            {
               NH[nf1][0]++;
               n=NH[nf1][0];
               NH[nf1][n]=nf2;
            }

            Here=0;
            for(n=1;n<=NH[nf2][0];n++) if(NH[nf2][n]==nf1) Here=1;
            if(Here==0)
            {
               NH[nf2][0]++;
               n=NH[nf2][0];
               NH[nf2][n]=nf1;
            }
         }
      }
   }
}

int DireccionXY(int Dir, int * x, int * y){
   switch(Dir) {
      case 1 :  *x= 0; *y= 1; break;
      case 2 :  *x= 1; *y= 0; break;
      case 3 :  *x= 0; *y=-1; break;
      case 4 :  *x=-1; *y= 0; break;
      default : 
         printf("Error en el codigod e direccion: %d\n",Dir);
         return 1;
   }
   return 0;
}

int ThereIsLateral(char D, int Lx, int Ly,int *** P ,int x, int y, int * lx, int * ly){
   int dx,dy=0;
   
   //check default left
   switch(P[x][y][DIR]) {
      case 1 :  dx=-1; dy= 0; break;
      case 2 :  dx= 0; dy= 1; break;
      case 3 :  dx= 1; dy= 0; break;
      case 4 :  dx= 0; dy=-1; break;
      default : 
         printf("Error en el codigod e direccion: %d\n",P[x][y][DIR]);
         return 0; //FALSE
   }
   //If right
   if (D == 'R'){
      dx*=-1;
      dy*=-1;
      }
   
   
   if ((x+dx > Lx) ||(y+dy > Ly) || (x+dx < 0) ||(y+dy < 0) || // Not Inbound
       P[dx][dy][MTYPE]==ACT ||                                // Not ACT
       P[dx][dy][MTYPE]==NONE ||                               // Not empty
       (((P[x][y][DIR]+P[dx][dy][DIR])%2) !=0 ))               // Not Parallel
      return 0; //FALSE
   
   return 0; //TRUE 
}

//~ int EtiquetaFA(int Lx,int Ly,int*** P,int nz,int** FA){
   //~ int x0=FA[nz][0],y0=FA[nz][1];
   //~ int dx,dy;
   //~ int d;
   //~ int L=FA[nz][3];
   
   //~ if (DireccionXY(FA[nz][2],&dx,&dy)
      //~ return 1;
   
   //~ for(d=0 ; d<L && xi>0 && xi<Lx && yi>0  && yi<Ly ; d++, xi=x0+d*dx , yi=y0+d*dy){              
      //~ for(b=-1;b<=1;b+=2){
         //~ if((x+b>-1)&&(x+b<Lx)) {
            //~ if((P[x+b][z][MTYPE]==MY0)&&(P[x+b][z][DIR]==1)&&(P[x+b][z][TAG]!=0)) {
            //~ nf=P[x+b][z][TAG];
            //~ P[x+b][z][ATP]=0;
            //~ FM[nf][2]=0;
            //~ }
         //~ }
      //~ }
//~ }

int RemoveMolecule(int x, int y, int*** P,int** FA,double* TA,int** FL,int** FM,int** NH){
   int nf; //TAG id
   
   switch(P[x][y][MTYPE]) {
      case ACT : 
            printf("ERROR: Removing ACT not implemented");
            break;
      
      case CRL : 
            if (P[x][y][TAG] != 0){        //MTYPE = CRL (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
               nf=P[x][y][TAG];              // Id CRL
               if(nf==FL[0][0]) FL[0][0]--;  // Total CRL num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta, no se actuliza FA[0][NCRL]
               P[x][y][TAG]=0;               // CRL fuera de listas
               FA[GLOBAL][NCRL]--;           // Total CRL FA[0][NCRL]--
               FL[nf][0]=0 ;                 // CRL suelto : no se elimina ???
               FL[nf][1]=0 ;                 // CRL suelto : no se elimina ???
               }
            break;
            
      case MYO :  
            if (P[x][y][TAG] != 0) {     //MTYPE = MYO (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
               nf=P[x][y][TAG];              // Id MYO
               if(nf==FM[0][0]) FM[0][0]--;  // Total MYO num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta
               P[x][y][4]=0;                 // MYO fuera de listas
               P[x][y][2]=0;                 // MYO sin ATP
               FA[0][3]--;                   // Total MYO FA[0][NMYO]--
               FM[nf][0]=0;                  // MYO suelto : no se elimina ???
               FM[nf][1]=0;                  // MYO suelto : no se elimina ???
               }
            break;
      default :
         printf("Unrecognized Molecule %d\n",P[x][y][MTYPE]);
         return 1;
   }
   return 0;   
}
int RemoveFilLateral(int Lx,int Ly,int nz,int*** P,int** FA,double* TA,int** FL,int** FM,int** NH){
   //Estado original del filamento
   int x0=FA[nz][0];  //Origen X
   int y0=FA[nz][1];  //Origen Y
   int MD=FA[nz][2];   //Direccion
   int ML=FA[nz][3];   //Longitud
   
   int dx,dy;          //Direccion desplazamiento (coordenadas)
   int d,xi,yi;
   int lx=0,ly=0;
   
   if (DireccionXY(MD,&dx,&dy))
      return 1;

#ifdef CHECK      
   if ((x0+ML*dx > Lx) ||(y0+ML*dy > Ly) || (x0+ML*dx < 0) ||(y0+ML*dy < 0)){
      printf("Error (%s : %d ): Filament lenght out of bound\n",__FILE__, __LINE__);
      return 1;
      }
#endif

   for(d=0, xi=x0 , yi=y0; d<ML ; d++, xi=x0+d*dx , yi=y0+d*dy){

#ifdef CHECK      
      if (P[xi][yi][MTYPE] != ACT ){
         printf("Error (%s : %d ): Next element %d is not ACT\n",__FILE__, __LINE__, d);
         return 1;
         }
      if (P[xi][yi][DIR] != MD ){
         printf("Error (%s : %d ): Next element %d have not the same direction\n",__FILE__, __LINE__, d);
         return 1;
         }
#endif

      if (ThereIsLateral('L',Lx,Ly,P,xi,yi,&lx,&ly)&&(P[lx][ly][TAG]!=0) ){
         RemoveMolecule(lx,ly,P,FA,TA,FL,FM,NH);
         }
      
      if (ThereIsLateral('R',Lx,Ly,P,xi,yi,&lx,&ly)&&(P[lx][ly][TAG]!=0) ){
         RemoveMolecule(lx,ly,P,FA,TA,FL,FM,NH);
         }
      }
      
   return 0;
}

//Rotura filamento
void FilBreak(int Lx,int Ly,double eA,double eL,double eM,double fTDis,int nz,int*** P,double** E,int** FA,double* TA,int** FL,int** FM,int** NH)
{
   int Mx0=FA[nz][0],My0=FA[nz][1],MD=FA[nz][2],ML=FA[nz][3];  //Estado original del filamento
   double MT=TA[nz];
   int x=0,y=0,z=0,R=0,nf=0;
   int a=0,b=0,c=0,d=0;
   double DESis=0;

//Direccion
   if((MD==1)||(MD==2)) a=1;
   else if((MD==3)||(MD==4)) a=-1;

//Punto de rotura en el centro del filamento  
   if((MD==1)||(MD==3)){
      x=Mx0;
      y=My0+(ML/2)*a;
      }
   else if((MD==2)||(MD==4)){
      x=Mx0+(ML/2)*a;
      y=My0;
      }

//Rotura
   R=En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,0,0,&DESis);  
   P[x][y][2]=0;                                
   En(Lx,Ly,eA,eL,eM,x,y,P,E,FM,1,R,&DESis);    

   if((MD==1)||(MD==3)){//Direccion Y
      R=En(Lx,Ly,eA,eL,eM,x,y+a,P,E,FM,0,0,&DESis); //Siempre se eliminan los dos enlaces de la proteina elegida ---> por que dos veces?
      P[x][y+a][2]=0;
      En(Lx,Ly,eA,eL,eM,x,y+a,P,E,FM,1,R,&DESis);
      if(ML>3){                                      //Si la longitud del filamento es mayor que 3 se rompe tambien otro enlace anterior
         R=En(Lx,Ly,eA,eL,eM,x,y-a,P,E,FM,0,0,&DESis);
         P[x][y-a][2]=0;
         En(Lx,Ly,eA,eL,eM,x,y-a,P,E,FM,1,R,&DESis);
         }
      }
   else if((MD==2)||(MD==4)){//Direccion X
      R=En(Lx,Ly,eA,eL,eM,x+a,y,P,E,FM,0,0,&DESis);
      P[x+a][y][2]=0;
      En(Lx,Ly,eA,eL,eM,x+a,y,P,E,FM,1,R,&DESis);
      if(ML>3){
         R=En(Lx,Ly,eA,eL,eM,x-a,y,P,E,FM,0,0,&DESis);
         P[x-a][y][2]=0;
         En(Lx,Ly,eA,eL,eM,x-a,y,P,E,FM,1,R,&DESis);
         }
      }

//Transmision de tension. Si tiene vecinos la tension se transmite equitativamente entre ellos
   if (NH[nz][0]!=0){ 
      MT*=fTDis;
      MT/=NH[nz][0];
      for(z=1;z<=NH[nz][0];z++) {
         nf=NH[nz][z];
         TA[nf]+=MT;
         }
      }
   TA[nz]=0;

//Correccion etiquetaje
   FA[0][1]-=ML;                       //Correccion sumatorio longitud
   if((MD==1)||(MD==3)){                  //Direccion Y
      if(abs(y-My0)>=4){                  //Existe filamento en el extremo debil despues de eliminar 3 enlaces de Actinas
         FA[nz][3]=abs(y-My0)-1;                //Longitud del filamento restante
         FA[0][1]+=FA[nz][3];                   //Sumatorio longitud 
         
         //Se consumen todos los ATP de la myosinas en el filamento
         for(d=0;d<FA[nz][3];d++){              //Correccion etiquetas filamento
            z=My0+d*a;
            for(b=-1;b<=1;b+=2){                // miramos laterales
               if((x+b>-1)&&(x+b<Lx)){
                  if((P[x+b][z][0]==3)&&(P[x+b][z][1]==1)&&(P[x+b][z][4]!=0)){  //MTYPE = MYO (CONNECTED, pero no sabemos a quien) , Solo se toma en cuanta una direccion ???
                     nf=P[x+b][z][4];  // Id MYO
                     P[x+b][z][2]=0;   // MYO sin ATP
                     FM[nf][2]=0;      // MYO sin ATP --> el unico sitio donde se actuliza FM
                     }
                  }
               }
            }

         P[x][y-a][4]=0;                        //Borrado etiquetas de la actina hidrolizada extra --> no se ha hecho antes ???

         for(b=-1;b<=1;b+=2){    /// Esto no es redundante ??? 
            if((x+b>-1)&&(x+b<Lx)){
               if((P[x+b][y-a][0]==2)&&(P[x+b][y-a][1]==1)&&(P[x+b][y-a][4]!=0)){ //MTYPE = CRL (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                  nf=P[x+b][y-a][4];             // Id CRL
                  if(nf==FL[0][0]) FL[0][0]--;   // Total CRL num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta, no se actuliza FA[0][NCRL]
                  P[x+b][y-a][4]=0;              // CRL fuera de listas
                  for(c=0;c<=1;c++) FL[nf][c]=0; // CRL suelto : no se elimina ???
                  FA[0][2]--;                    // Total CRL FA[0][NCRL]--
                  }
               
               else if((P[x+b][y-a][0]==3)&&(P[x+b][y-a][1]==1)&&(P[x+b][y-a][4]!=0)) { //MTYPE = MYO (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                  nf=P[x+b][y-a][4];            // Id MYO
                  if(nf==FM[0][0]) FM[0][0]--;  // Total MYO num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta
                  P[x+b][y-a][4]=0;             // MYO fuera de listas
                  P[x+b][y-a][2]=0;             // MYO sin ATP

                  for(c=0;c<=1;c++) FM[nf][c]=0;// MYO suelto : no se elimina ???
                  FA[0][3]--;                   // Total MYO FA[0][NMYO]--
                  }
               }
            }
         /// No se actualiza el tamano total de filamentos ni el tamno del filamento !!!
         } 
      else { 
         // NO Existe filamento en el extremo debil despues de eliminar 3 enlaces de Actinas
         for(c=0;c<=3;c++) FA[nz][c]=0;  // Eliminamos filamento ACT
         for(d=0;d<abs(y-My0);d++){
            z=My0+d*a;     //Siguente ACT
            P[x][z][4]=0;  //Siguente ACT, TAG=0
            for(b=-1;b<=1;b+=2){ //Laterales
               if((x+b>-1)&&(x+b<Lx)) {
                  if((P[x+b][z][0]==2)&&(P[x+b][z][1]==1)&&(P[x+b][z][4]!=0)) { //MTYPE = CRL (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                     nf=P[x+b][z][4];              // Id CRL
                     if(nf==FL[0][0]) FL[0][0]--;  // Total CRL num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta, no se actuliza FA[0][NCRL]
                     P[x+b][z][4]=0;               // CRL fuera de listas
                     for(c=0;c<=1;c++) FL[nf][c]=0;// CRL suelto : no se elimina ???
                     FA[0][2]--;                   // Total CRL FA[0][NCRL]--
                     }
                  else if((P[x+b][z][0]==3)&&(P[x+b][z][1]==1)&&(P[x+b][z][4]!=0)){ //MTYPE = MYO (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                     nf=P[x+b][z][4];              // Id MYO
                     if(nf==FM[0][0]) FM[0][0]--;  // Total MYO num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta, no se actuliza FA[0][NMYO]
                     P[x+b][z][4]=0;               // MYO anterior fuera de listas
                     P[x+b][z][2]=0;               // MYO anterior sin ATP 
                     for(c=0;c<=3;c++) FM[nf][c]=0;// MYO suelto : no se elimina ???
                     FA[0][3]--;                   // Total MYO FA[0][NMYO]--
                     }
                  }
               }
            }
         }

      // Tratemos el resto del filamento del extremo fuerte
      if(ML-abs(y-My0)>=4){                  // Si es filamento
         P[x][y][4]=0;                       // Eliminamos Actina seleccionada
         for(b=-1;b<=1;b+=2){                // Laterales
            if((x+b>-1)&&(x+b<Lx)){   
               if((P[x+b][y][0]==2)&&(P[x+b][y][1]==1)&&(P[x+b][y][4]!=0)){ //MTYPE = CRL (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                  nf=P[x+b][y][4];                 // Id CRL
                  if(nf==FL[0][0]) FL[0][0]--;     // Total CRL num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta, no se actuliza FA[0][NCRL]
                  P[x+b][y][4]=0;                  // CRL fuera de listas
                  for(c=0;c<=1;c++) FL[nf][c]=0;   // CRL suelto : no se elimina ???
                  FA[0][2]--;                      // Total CRL FA[0][NCRL]--
                  }
               else if((P[x+b][y][0]==3)&&(P[x+b][y][1]==1)&&(P[x+b][y][4]!=0)) { //MTYPE = MYO (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                  nf=P[x+b][y][4];                 // Id MYO
                  if(nf==FM[0][0]) FM[0][0]--;     // Total MYO num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta, no se actuliza FA[0][NMYO]
                  P[x+b][y][4]=0;                  // MYO anterior fuera de listas
                  P[x+b][y][2]=0;                  // MYO anterior sin ATP 
                  for(c=0;c<=3;c++) FM[nf][c]=0;   // MYO suelto : no se elimina ???
                  FA[0][3]--;                      // Total MYO FA[0][NMYO]--
               }
            }
         }

         for(nz=1;nz<=FA[0][0];nz++) if(FA[nz][3]==0) break;  // Id nuevo filamento
         if(nz>FA[0][0]) FA[0][0]=nz;                         // Actualizamos total filamentos
         TA[0]=FA[0][0];

         FA[nz][0]=x;                 // Punto init x del nuevo filamento
         FA[nz][1]=y+a;               // Punto init y del nuevo filamento
         FA[nz][2]=P[x][y+a][1];      // Direccion del nuevo filamento es la misma
         FA[nz][3]=ML-abs(y-My0)-1;   // Tamano del nuevo filamento
         FA[0][1]+=FA[nz][3];         // Insertamos Tamano del nuevo filamento

         for(d=1;d<ML-abs(y-My0);d++) {  // ERROR en el for, fataba -1
            z=y+d*a;
            /* SEG FAULT !!! z=-1*/
            P[x][z][4]=nz;            // Insertamos nuevo Id del filamento
            for(b=-1;b<=1;b+=2){      // Laterales
               if((x+b>-1)&&(x+b<Lx)) {
                  if((P[x+b][z][0]==2)&&(P[x+b][z][1]==1)&&(P[x+b][z][4]!=0)){ //MTYPE = CRL (NOT CONNECTED) + mismo error direccion ???
                     nf=P[x+b][z][4];              // Id CRL
                     FL[nf][0]=nz;                 // Insertamos Id Filamento en CRL  ???
                     FL[nf][1]=P[x+(2*b)][z][4];   // Actualizamos segundo enlace CRL ???
                     }
                  else if((P[x+b][z][0]==3)&&(P[x+b][z][1]==1)&&(P[x+b][z][4]!=0)) { //MTYPE = MYO (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                     nf=P[x+b][z][4];              // Id MYO
                     FM[nf][0]=nz;                 // Insertamos Id Filamento en MYO  ???
                     FM[nf][1]=P[x+(2*b)][z][4];   // Actualizamos segundo enlace MYO ???
                     FM[nf][2]=0;                  // Quitamos ATP MYO ???
                     P[x+b][z][2]=0;               // Quitamos ATP MYO en P ???
                     }
                  }
               }
            }
         }
      else { // Si no es filamento
         for(d=0;d<ML-abs(y-My0);d++) {  // ERROR distancias el el for
            z=y+d*a;
            P[x][z][4]=0;                    // Eliminamos Actina
            for(b=-1;b<=1;b+=2) {            // Laterales
               if((x+b>-1)&&(x+b<Lx)) {
                  if((P[x+b][z][0]==2)&&(P[x+b][z][1]==1)&&(P[x+b][z][4]!=0)) { //MTYPE = CRL (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                     nf=P[x+b][z][4];                 // Id CRL
                     if(nf==FL[0][0]) FL[0][0]--;     // Total CRL num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta
                     P[x+b][z][4]=0;                  // CRL fuera de listas
                     for(c=0;c<=1;c++) FL[nf][c]=0;   // CRL suelto : no se elimina ???
                     FA[0][2]--;                      // FA[0][NCRL]--
                     }
                  else if((P[x+b][z][0]==3)&&(P[x+b][z][1]==1)&&(P[x+b][z][4]!=0)) { //MTYPE = MYO (CONNECTED, pero no sabemos a quien) + mismo error direccion ???
                     nf=P[x+b][z][4];                // Id MYO
                     if(nf==FM[0][0]) FM[0][0]--;    // Total MYO num conn -- ---> Porque hay que comparar ??? se puede perder la cuenta
                     P[x+b][z][4]=0;                 // MYO fuera de listas
                     P[x+b][z][2]=0;                 // MYO ATP=0
                     for(c=0;c<=3;c++) FM[nf][c]=0;  // MYO suelto : no se elimina ???
                     FA[0][3]--;                     // FA[0][NMYO]--
                  }
               }
            }
         }
      }
   }
   else if((MD==2)||(MD==4)) //Direccion X
   {
      if(abs(x-Mx0)>=4)
      {
         FA[nz][3]=abs(x-Mx0)-1;
         FA[0][1]+=FA[nz][3];

         for(d=0;d<FA[nz][3];d++)
         {
            z=Mx0+d*a;

            for(b=-1;b<=1;b+=2)
            {
               if((y+b>-1)&&(y+b<Ly))
               {
                  if((P[z][y+b][0]==3)&&(P[z][y+b][1]==2)&&(P[z][y+b][4]!=0))
                  {
                     nf=P[z][y+b][4];
                     FM[nf][2]=0;
                     P[z][y+b][2]=0;
                  }
               }
            }
         }

         P[x-a][y][4]=0;

         for(b=-1;b<=1;b+=2)
         {
            if((y+b>-1)&&(y+b<Ly))
            {
               if((P[x-a][y+b][0]==2)&&(P[x-a][y+b][1]==2)&&(P[x-a][y+b][4]!=0))
               {
                  nf=P[x-a][y+b][4];
                  if(nf==FL[0][0]) FL[0][0]--;
                  P[x-a][y+b][4]=0;

                  for(c=0;c<=1;c++) FL[nf][c]=0;
                  FA[0][2]--;
               }
               else if((P[x-a][y+b][0]==3)&&(P[x-a][y+b][1]==2)&&(P[x-a][y+b][4]!=0))
               {
                  nf=P[x-a][y+b][4];
                  if(nf==FM[0][0]) FM[0][0]--;
                  P[x-a][y+b][4]=0;
                  P[x-a][y+b][2]=0;

                  for(c=0;c<=1;c++) FM[nf][c]=0;
                  FA[0][3]--;
               }
            }
         }
      }
      else
      {
         for(c=0;c<=3;c++) FA[nz][c]=0;

         for(d=0;d<abs(x-Mx0);d++)
         {
            z=Mx0+d*a;
            P[z][y][4]=0;

            for(b=-1;b<=1;b+=2)
            {
               if((y+b>-1)&&(y+b<Ly))
               {
                  if((P[z][y+b][0]==2)&&(P[z][y+b][1]==2)&&(P[z][y+b][4]!=0))
                  {
                     nf=P[z][y+b][4];
                     if(nf==FL[0][0]) FL[0][0]--;
                     P[z][y+b][4]=0;

                     for(c=0;c<=1;c++) FL[nf][c]=0;
                     FA[0][2]--;
                  }
                  else if((P[z][y+b][0]==3)&&(P[z][y+b][1]==2)&&(P[z][y+b][4]!=0))
                  {
                     nf=P[z][y+b][4];
                     if(nf==FM[0][0]) FM[0][0]--;
                     P[z][y+b][4]=0;
                     P[z][y+b][2]=0;

                     for(c=0;c<=3;c++) FM[nf][c]=0;
                     FA[0][3]--;
                  }
               }
            }
         }
      }

      if(ML-abs(x-Mx0)>=4) // ERROR dist nuevo fil
      {
         P[x][y][4]=0;

         for(b=-1;b<=1;b+=2)
         {
            if((y+b>-1)&&(y+b<Ly))
            {
               if((P[x][y+b][0]==2)&&(P[x][y+b][1]==2)&&(P[x][y+b][4]!=0))
               {
                  nf=P[x][y+b][4];
                  if(nf==FL[0][0]) FL[0][0]--;
                  P[x][y+b][4]=0;

                  for(c=0;c<=1;c++) FL[nf][c]=0;
                  FA[0][2]--;
               }
               else if((P[x][y+b][0]==3)&&(P[x][y+b][1]==2)&&(P[x][y+b][4]!=0))
               {
                  nf=P[x][y+b][4];
                  if(nf==FM[0][0]) FM[0][0]--;
                  P[x][y+b][4]=0;
                  P[x][y+b][2]=0;

                  for(c=0;c<=3;c++) FM[nf][c]=0;
                  FA[0][3]--;
               }
            }
         }

         for(nz=1;nz<=FA[0][0];nz++) if(FA[nz][3]==0) break;
         if(nz>FA[0][0]) FA[0][0]=nz;
         TA[0]=FA[0][0];

         FA[nz][0]=x+a;
         FA[nz][1]=y;
         FA[nz][2]=P[x+a][y][1];
         FA[nz][3]=ML-abs(x-Mx0)-1;
         FA[0][1]+=FA[nz][3];

         for(d=1;d<ML-abs(x-Mx0);d++) // ERROR dist nuevo fil
         {
            z=x+d*a;
            P[z][y][4]=nz;

            for(b=-1;b<=1;b+=2)
            {
               if((y+b>-1)&&(y+b<Ly))
               {
                  if((P[z][y+b][0]==2)&&(P[z][y+b][1]==2)&&(P[z][y+b][4]!=0))
                  {
                     nf=P[z][y+b][4];
                     FL[nf][0]=nz;
                     FL[nf][1]=P[z][y+(2*b)][4];
                  }
                  else if((P[z][y+b][0]==3)&&(P[z][y+b][1]==2)&&(P[z][y+b][4]!=0))
                  {
                     nf=P[z][y+b][4];
                     FM[nf][0]=nz;
                     FM[nf][1]=P[z][y+(2*b)][4];
                     FM[nf][2]=0;
                     P[z][y+b][2]=0;
                  }
               }
            }
         }
      }
      else
      {
         for(d=0;d<ML-abs(x-Mx0);d++)// ERROR dist nuevo fil
         {
            z=x+d*a;
            P[z][y][4]=0;

            for(b=-1;b<=1;b+=2)
            {
               if((y+b>-1)&&(y+b<Ly))
               {
                  if((P[z][y+b][0]==2)&&(P[z][y+b][1]==2)&&(P[z][y+b][4]!=0))
                  {
                     nf=P[z][y+b][4];
                     if(nf==FL[0][0]) FL[0][0]--;
                     P[z][y+b][4]=0;

                     for(c=0;c<=1;c++) FL[nf][c]=0;
                     FA[0][2]--;
                  }
                  else if((P[z][y+b][0]==3)&&(P[z][y+b][1]==2)&&(P[z][y+b][4]!=0))
                  {
                     nf=P[z][y+b][4];
                     if(nf==FM[0][0]) FM[0][0]--;
                     P[z][y+b][4]=0;
                     P[z][y+b][2]=0;

                     for(c=0;c<=3;c++) FM[nf][c]=0;
                     FA[0][3]--;
                  }
               }
            }
         }
      }
   }
}

//Registro de la lista de Bundles
void ControlBundle(int* DP,int* CDP,int Lx, int Ly,int*** P,int** FA,int** FL,int** FM,int** NH,int** BA,double * TA,int NFA,int Cheese,int TBreak)
{
   if(*CDP==Cheese)
   {
      (*DP)++;
      *CDP=0;

      FILE* docBA;         //Archivos ActinaBundle
      FILE* docNBA;        //Archivos Actina
      FILE* docCL;         //Archivos C-Linker
      FILE* docMY;         //Archivos Myosin
      FILE* docATA;        //Archivos ActinaBundle
      FILE* docABA;        //Archivos tamano bundle all
      FILE* docDBA;        //Archivos tamano bundle
      FILE* docTFA;        //Archivos averagetension filamentos
      FILE* docDFA;        //Archivos control filamentos
      FILE* docBundles;    //Archivos control vecinos

      char dBA[]="DataA_0000_0000.dat";       //Cadenas de texto para nombrar archivos
      char dNBA[]="DataNBA_0000.dat";
      char dCL[]="DataCL_0000.dat";
      char dMY[]="DataMY_0000.dat";
      char dATA[]="DataATA_0000_0000.dat";
      char dDBA[]="DataBA_0000.dat";
      char dTFA[]="DataTA_0000.dat";
      char dDFA[]="DataFA_0000.dat";
      char dBundles[]="DataBundles_0000.dat";

      int x=0,y=0,n=0,m=0,k=0,l=0;
      int nba=0;
      int Here=0;
      int time=(*DP)*Cheese;

//Control Filamentos
      sprintf(dDFA,"DataFA_%d.dat",*DP);
      docDFA=fopen(dDFA,"a");
      fprintf(docDFA,"Fil:%d\nA:%d\nCL:%d\nMyo:%d\n",FA[0][0],FA[0][1],FA[0][2],FA[0][3]);
      fprintf(docDFA,"\nLista Filamentos\n");

      fprintf(docDFA,"nf\t[x0][y0]\tD\tL\tT\n");
      for(x=1;x<NFA;x++) fprintf(docDFA,"%d\t[%d][%d]\t%d\t%d\t%.1lf\n",x,FA[x][0],FA[x][1],FA[x][2],FA[x][3],TA[x]);
      fclose(docDFA);

//Control Bundles
      for(x=0;x<(NFA/3);x++) for(y=0;y<NFA;y++) BA[x][y]=0;
      for(x=0;x<NFA;x++) for(y=0;y<NFA;y++) NH[x][y]=0;

      Neighbor(Lx,Ly,P,FA,FL,FM,NH);

      for(x=1;x<=NH[0][0];x++)
      {
         if(NH[x][0]!=0)
         {
            Here=0;
            for(n=1;n<=BA[0][0];n++) for(m=1;m<=BA[n][0];m++) if(x==BA[n][m]) Here=1;
            if(Here==0)
            {
               BA[0][0]++;
               nba=BA[0][0];
               BA[nba][0]=1;
               BA[nba][1]=x;
               BA[0][nba]=FA[x][3];
               for(l=1;l<=BA[nba][0];l++)
               {
                  for(k=1;k<=NH[BA[nba][l]][0];k++)
                  {
                     Here=0;
                     for(m=1;m<=BA[nba][0];m++) if(NH[BA[nba][l]][k]==BA[nba][m]) Here=1;
                     if(Here==0)
                     {
                        BA[nba][0]++;
                        BA[nba][BA[nba][0]]=NH[BA[nba][l]][k];
                        BA[0][nba]+=FA[NH[BA[nba][l]][k]][3];
                     }
                  }
               }
            }
         }
      }

//Impresion lista Bundles
      sprintf(dBundles,"DataBundles_%d.dat",*DP);
      docBundles=fopen(dBundles,"a");
      for(x=1;x<(NFA/3);x++)
      {
         fprintf(docBundles,"nf:%d\n",x);
         fprintf(docBundles,"Lista de Vecinos:\n");
         for(y=1;y<NFA;y++)
         {
            fprintf(docBundles,"%d\t",BA[x][y]);
            if(BA[x][y]==0)
            {
               fprintf(docBundles,"\n\n");
               break;
            }
         }
      }
      fclose(docBundles);

//Impresion tamano Bundles
      int* NaBu;
      int conB=0,nBu=0;
      double TFM=0;

      NaBu=(int*)malloc((BA[0][0]+1)*sizeof(int));
      for(x=0;x<=BA[0][0];x++) NaBu[x]=0;

      sprintf(dDBA,"DataBA_%d.dat",*DP);
      docDBA=fopen(dDBA,"a");

      for(n=1;n<=BA[0][0];n++)fprintf(docDBA,"%i %i \n",n,BA[0][n]);
      for(n=1;n<=9;n++)
      {
         sprintf(dTFA,"DataTA_%d.dat",n);
         docTFA=fopen(dTFA,"a");
         if(BA[n][0]!=0)
         {
            for(m=1;m<=BA[n][0];m++) TFM+=TA[BA[n][m]];
            TFM/=BA[n][0];
         }
         fprintf(docTFA,"%i %lf \n",time,TFM);
         fclose(docTFA);
      }
      fclose(docDBA);

      docABA=fopen("DataABA.dat","a");
      fprintf(docABA,"%i %i %i %i %i %i %i %i %i %i \n",time,BA[0][1],BA[0][2],BA[0][3],BA[0][4],BA[0][5],BA[0][6],BA[0][7],BA[0][8],BA[0][9]);
      fclose(docABA);

//Impresion foto Bundles
      for(x=0;x<Lx;x++)
      {
         for(y=0;y<Ly;y++)
         {
            if((P[x][y][0]==1)&&(P[x][y][4]!=0))
            {
               conB=0;
               for(n=1;n<=BA[0][0];n++)
               {
                  for(m=1;m<=BA[n][0];m++)
                  {
                     if(BA[n][m]==P[x][y][4])
                     {
                        if(NaBu[n]==0)
                        {
                           nBu++;
                           NaBu[n]=nBu;
                        }
                        conB=1;
                        break;
                     }
                  }
                  if(conB==1) break;
               }

               if(conB==1)
               {
                  sprintf(dBA,"DataA_%d_%d.dat",*DP,NaBu[n]);
                  docBA=fopen(dBA,"a");
                  fprintf(docBA,"%i %i \n",x,y);
                  fclose(docBA);
               }
               else
               {
                  sprintf(dNBA,"DataNBA_%d.dat",*DP);
                  docNBA=fopen(dNBA,"a");
                  fprintf(docNBA,"%i %i \n",x,y);
                  fclose(docNBA);

               }
            }
            else if(P[x][y][0]==1)
            {
               sprintf(dNBA,"DataNBA_%d.dat",*DP);
               docNBA=fopen(dNBA,"a");
               fprintf(docNBA,"%i %i \n",x,y);
               fclose(docNBA);

            }
            else if(P[x][y][0]==2)
            {
               sprintf(dCL,"DataCL_%d.dat",*DP);
               docCL=fopen(dCL,"a");
               fprintf(docCL,"%i %i \n",x,y);
               fclose(docCL);
            }
            else if(P[x][y][0]==3)
            {
               sprintf(dMY,"DataMY_%d.dat",*DP);
               docMY=fopen(dMY,"a");
               fprintf(docMY,"%i %i \n",x,y);
               fclose(docMY);
            }
         }
      }
      free(NaBu);

      int Tint=TBreak/5;
//Impresion foto tension
      for(x=0;x<Lx;x++)
      {
         for(y=0;y<Ly;y++)
         {
            if(P[x][y][0]==1)
            {
               if(P[x][y][4]!=0)
               {
                  if(TA[P[x][y][4]]<=Tint)
                  {
                     sprintf(dATA,"DataA_%d_%dT.dat",*DP,1);
                     docATA=fopen(dATA,"a");
                     fprintf(docATA,"%i %i \n",x,y);
                     fclose(docATA);
                  }
                  else if(TA[P[x][y][4]]<=2*Tint)
                  {
                     sprintf(dATA,"DataA_%d_%dT.dat",*DP,2);
                     docATA=fopen(dATA,"a");
                     fprintf(docATA,"%i %i \n",x,y);
                     fclose(docATA);
                  }
                  else if(TA[P[x][y][4]]<=3*Tint)
                  {
                     sprintf(dATA,"DataA_%d_%dT.dat",*DP,3);
                     docATA=fopen(dATA,"a");
                     fprintf(docATA,"%i %i \n",x,y);
                     fclose(docATA);
                  }
                  else if(TA[P[x][y][4]]<=4*Tint)
                  {
                     sprintf(dATA,"DataA_%d_%dT.dat",*DP,4);
                     docATA=fopen(dATA,"a");
                     fprintf(docATA,"%i %i \n",x,y);
                     fclose(docATA);
                  }
                  else if(TA[P[x][y][4]]<=5*Tint)
                  {
                     sprintf(dATA,"DataA_%d_%dT.dat",*DP,5);
                     docATA=fopen(dATA,"a");
                     fprintf(docATA,"%i %i \n",x,y);
                     fclose(docATA);
                  }
               }
            }
         }
      }
   }
}


