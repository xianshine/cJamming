#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float distance(float *  r1,float * r2, float box) {
   int i;
   float dr[3],dist2;
   

   dist2=0;
   for (i=0;i<3;i++) {
      dr[i]=fabs(r1[i]-r2[i]);
      if (dr[i]>0.5*box) dr[i]=box-dr[i];
      dist2+=pow(dr[i],2);
   }
   return sqrt(dist2);
}

float pairP(float r, float sigma, float eps) {
   float energy;
   if (r<sigma) energy=0.5*eps*pow((1-r/sigma),2);
   else
      energy=0;

   return energy;
}

int main(int argc, char ** argv)
{
   int i1, i2, i,j, it1,it2;

/*
   printf( "argc = %d\n", argc );
   for( i = 0; i < argc; ++i ) {
        printf( "argv[ %d ] = %s\n", i, argv[ i ] );
   }
*/
   
   int nB;
   float xlo,xhi,xyz[100000][3],box;
   char iname[20],oname[20];
   FILE *fi,*fo;

   nB=atoi(argv[1]);
   xlo=atof(argv[2]);
   xhi=atof(argv[3]);
   strcpy(iname,argv[4]);
   strcpy(oname,argv[5]);
   box=xhi-xlo;
   //printf("%d %f %f\n",nB,xlo,xhi);

   //float matDist[nB][nB];
   float dist;

   if ((fi=fopen(iname, "r"))==NULL)
   {
      printf ("Can't open file %s!\n",iname);
      exit(0);
   }
   else {
      j=0;
      //printf("opened file %s!\n",iname);
   }
      

   for (i1=0;i1<nB;i1++) {
      fscanf (fi, "%f %f %f", &xyz[i1][0],&xyz[i1][1],&xyz[i1][2]);
      //printf("%d %d %f %f %f\n", j,types[i1],xyz[i1][0],xyz[i1][1],xyz[i1][2]);
   }

/*
   for (i1=0;i1<nB;i1++) {
      if (i1%100==0) printf("Now at %d \n",i1);
      for (i2=i1+1;i2<nB;i2++) {
         dist=distance(xyz[i1],xyz[i2],box);
         matDist[i1][i2]=dist;
         matDist[i2][i1]=dist;
      }
   }      
*/      
   if ((fo=fopen(oname, "w"))==NULL)
   {
      printf ("Can't open file %s!\n",oname);
      exit(0);
   }
   else {
      j=0;
      //printf("opened file %s!\n",oname);
   }
   for (i1=0;i1<nB;i1++) {
      for (i2=0;i2<nB;i2++) {
         dist=distance(xyz[i1],xyz[i2],box);
         fprintf(fo, "%15.10f ",dist);
      }
      fprintf(fo,"\n");
   }
      fclose(fi);
      fclose(fo);
}
