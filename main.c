#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "inc.h"
/*-------------------- computation of one iteration*/

void OneStep (PointProb** myMatrix, double** Pij, double** old_Pij,
	      DomainPtr Dom,  MPI_Comm comm2D);
/*--------------------normalize to get prob. vector*/
int normliz(double** x, int ni, int nj, MPI_Comm comm2D);
/*--------------------allReduce calculates inf for convergence*/
double ComputErr(double** x, double** y, int ni, int nj, MPI_Comm comm2D);
/*--------------------generates matrix  probabs.. */
void get_prob (int i, int j, int ni, int nj, int qi, int qj, int *dims,
	       double *prob) ;
/*--------------------auxilary function*/

int main (int argc, char **argv){
 int ni, nj, westNB, eastNB, northNB, southNB;
 int myid, procs, itmax=10000 ;
 DomainPtr Dom; 
 MPI_Comm comm2D;
/*-------------------- main matrix : stores local probabs + matrix itself. */
 PointProb** myMatrix;
/*--------------------probability that particle is at this point*/
/*-------------------- AND previous value of Pij*/
 double **Pij, **old_Pij;
 int dims[2], coords[2], periods[2];
 int qi, qj;   /* virtual process dimensions  */  
 int i,j,k, ierr = 0;
 int it;   /* iteration count */
 double locprob[5];
 /* double t, tot;  DEBUG */
 char str1[20], str2[20];
 FILE *fo;  
 
 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &procs);
 MPI_Comm_rank(MPI_COMM_WORLD, &myid);
/*-------------------- Read proc-mesh dimensions */ 
 if (myid == 0) {
   if (argc == 2) {
     qi = atoi(argv[1]);
     if (procs % qi){
       printf(" ERROR: qi = %d does not divide procs = %d \n",qi,procs);
       ierr =1;
     }
     qj = procs/qi;
     dims[0] =  qi; 
     dims[1] =  qj;
   } else {
     printf(" enter value of qi after main.ex \n") ;
     ierr =2;
   }
 }
 //
 MPI_Bcast(&ierr, 1, MPI_INT, 0, MPI_COMM_WORLD);
 if (ierr>0) {
   return(ierr);
   MPI_Finalize();
   exit(1);
 }
/*-------------------- broadcast proc-mesh dimensions */ 
 MPI_Bcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
 qi = dims[0] ; 
 qj = dims[1];
 if (procs != qi*qj) {
   if (myid == 0) 
     printf(" ERROR: qi*qj does not equal nprocs - aborting \n");
   exit(1);
 }
/*-------------------- hard-wired dimensions of mesh */
 ni =  10;
 nj =  10;
/*-------------------- output files. Need to be different for each PE */
 strcpy(str2,"OUT/out_");
 sprintf(str1,"%d",myid);
 strcat(str2,str1);
 fo = fopen(str2, "w");
/*-------------------- CREATE cartesian communicator */
 periods[0]=periods[1]=1;
/*--------------------2D topoology created*/
 MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&comm2D);
 /*-------------------- my rank in 2-D form */
 MPI_Cart_coords(comm2D,myid,2,coords);
 qi = coords[0];
 qj = coords[1];
 printf(" (qi, qj) (%d %d ) == myid %d\n",qi,qj,myid);
/*-------------------- rank of the west neighbor of the node 
                       standard array  indexing */
 coords[0]=qi;
 coords[1]=qj-1;
 MPI_Cart_rank(comm2D,coords,&westNB);
/*-------------------- rank of the east neighbor of the node*/
 coords[0]=qi; 
 coords[1]=qj+1; 
 MPI_Cart_rank(comm2D,coords,&eastNB);
/*-------------------- rank of the south neighbor of the node*/
 coords[0]=qi+1;
 coords[1]=qj; 
 MPI_Cart_rank(comm2D,coords,&southNB);
/*-------------------- rank of the north neighbor of the node*/
 coords[0]=qi-1;
 coords[1]=qj; 
 MPI_Cart_rank(comm2D,coords,&northNB);
 /*-------------------- create Domain Dom */
 Dom = (DomainPtr) malloc(sizeof(myDomain));
/*-------------------- create ghost cell arrays */
 Dom->south = (double*)malloc(sizeof(double)*nj); 
 Dom->north = (double*)malloc(sizeof(double)*nj); 
 Dom->east  = (double*)malloc(sizeof(double)*ni); 
 Dom->west  = (double*)malloc(sizeof(double)*ni); 
/*-------------------- assign other fields */
 Dom->ni = ni;
 Dom->nj = nj;
 Dom->qi = qi;
 Dom->qj = qj;
 Dom->westNB = westNB;
 Dom->eastNB = eastNB;
 Dom->southNB = southNB;
 Dom->northNB = northNB;
/*----------------------------------------------------------------------*/
/*-------------------- create matrix and local probabilities */  
 myMatrix=(PointProb**)malloc(sizeof(PointProb*)*ni);  
 Pij=(double**)malloc(sizeof(double*)*ni);  
 old_Pij=(double**)malloc(sizeof(double*)*ni);
/*-------------------- allocate for each row */
 for(i=0;i<ni;i++) {
   myMatrix[i]=(PointProb*)malloc(sizeof(PointProb)*nj);
   Pij[i] = (double*) malloc(sizeof(double)*nj);
   old_Pij[i] = (double*) malloc(sizeof(double)*nj);
   //   
   for(j=0;j<nj;j++){
     get_prob(i, j, ni, nj, qi, qj, dims, locprob);
     for (k=0; k<5; k++) 
       myMatrix[i][j].prob[k] = locprob[k] ;
     Pij[i][j]     = 0.0; 
     old_Pij[i][j] = 0.0; 
   }
 }
/*---------------------------------------------------------------------*/
/*-------------------- set initial probability */
 if (qi == 0 && qj == 0) 
   old_Pij[0][0] = 1.0; 
/*-------------------- all info has been set now */
 int count = 0;
 
 /* ----- MPI_Wtime */
 double starttime, endtime;
 starttime = MPI_Wtime();

 for (it = 0; it < itmax; it++) {
   count = it;

   //if (qi == 0 && qj == 0) 
   //printf("iteration %d:\n",count);
/*-------------------- Do one iteration  */
   OneStep (myMatrix, Pij, old_Pij, Dom, comm2D);
/*-------------------- this is optional - can scale at end */
   normliz(Pij, ni, nj, comm2D);
/*-------------------- check convergence at least 8 digits or so*/   
   if (ComputErr(Pij, old_Pij, ni, nj, comm2D) < 1.e-09)
     break;
/*------------------ else copy and repeat */
   for(i=0; i<ni; i++) 
     for (j=0;j<nj;j++)
       old_Pij[i][j] = Pij[i][j]; 
 }
/*-------------------- Write data to files */
 for(i=0; i<ni; i++) {
   for (j=0;j<nj;j++)
     fprintf(fo,"%e ",Pij[i][j]) ;
   fprintf(fo,"\n"); 
 }

 printf("Stopped after %d iterations\n",count);
 endtime   = MPI_Wtime();
 printf("proc %d took %f seconds\n",myid, endtime-starttime);
 /*-------------------- free arrays */
 for (i=0; i< ni; i++){
   free (myMatrix[i]);
   free (Pij[i]);
   free (old_Pij[i]);
 }
 free (myMatrix);
 free (Pij);
 free (old_Pij);  
 free (Dom->south) ;
 free (Dom->north) ;
 free (Dom->east) ;
 free (Dom->west);
 free (Dom);
 MPI_Finalize();
 return(0);
}   

