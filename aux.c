#include <mpi.h>
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "inc.h"

void OneStep (PointProb** myMatrix, double** Pij, double** old_Pij,
	      DomainPtr Dom,  MPI_Comm comm2D){
/*-------------------- Does one iteration */
  int i,j;
  int ni =Dom->ni, nj=Dom->nj;
  //int qi = Dom->qi, qj = Dom->qj;
  int westNB = Dom->westNB, eastNB = Dom->eastNB, southNB = Dom->southNB, northNB = Dom->northNB;


  // Create sending buffer of ghost cell arrays
  double* toSouth = (double*)malloc(sizeof(double)*nj);
  double* toNorth = (double*)malloc(sizeof(double)*nj);
  double* toEast = (double*)malloc(sizeof(double)*ni);
  double* toWest = (double*)malloc(sizeof(double)*ni);

  // Initialize Pij
  for(int i = 0; i < ni; i++) {
  	for (int j = 0; j < nj; j++) {
		Pij[i][j] = 0;
	}
  }

  // Initialize boundary buffer
  for (int i = 0; i < ni; i++) {
        Dom->east[i] = 0;
        Dom->west[i] = 0;
	toEast[i] = 0;
	toWest[i] = 0;
  }
  for (int j = 0; j < nj; j++) {
        Dom->north[j] = 0;
        Dom->south[j] = 0;
	toNorth[j] = 0;
	toSouth[j] = 0;
  }


  // Main Loop
  for (i=0; i<ni; i++) {
  	for (j=0; j<nj; j++) {
		double pc = myMatrix[i][j].prob[0];
	       	double pw = myMatrix[i][j].prob[1];
		double pe = myMatrix[i][j].prob[2];
	        double pn = myMatrix[i][j].prob[3];
	       	double ps = myMatrix[i][j].prob[4];
		Pij[i][j] += pc * old_Pij[i][j];

		// NW corner
		if (i == 0 && j == 0) {
			toWest[i] += pw * old_Pij[i][j];
                        Pij[i][j+1] += pe * old_Pij[i][j];
                        toNorth[j] += pn * old_Pij[i][j];
                        Pij[i+1][j] += ps * old_Pij[i][j];
		}
		
		// NE corner
		else if (i == 0 && j == nj-1) {
			Pij[i][j-1] += pw * old_Pij[i][j];
                        toEast[i] += pe * old_Pij[i][j];
                        toNorth[j] += pn * old_Pij[i][j];
                        Pij[i+1][j] += ps * old_Pij[i][j];
		}

		// SW corner 
		else if (i == ni-1 && j == 0) {
                        toWest[i] += pw * old_Pij[i][j];
                        Pij[i][j+1] += pe * old_Pij[i][j];
                        Pij[i-1][j] += pn * old_Pij[i][j];
                        toSouth[j] += ps * old_Pij[i][j];
                }

		// SE corner
		else if (i == ni-1 && j == nj-1) {
                        Pij[i][j-1] += pw * old_Pij[i][j];
                        toEast[i] += pe * old_Pij[i][j];
                        Pij[i-1][j] += pn * old_Pij[i][j];
                        toSouth[j] += ps * old_Pij[i][j];
                }

		// North buffer, storing data to be sent to north processor
		else if (i == 0) {
			Pij[i][j-1] += pw * old_Pij[i][j];
                        Pij[i][j+1] += pe * old_Pij[i][j];
                        toNorth[j] += pn * old_Pij[i][j];
                        Pij[i+1][j] += ps * old_Pij[i][j];
		}

		// South buffer, storing data to be sent to south processor
		else if (i == ni-1) {
                        Pij[i][j-1] += pw * old_Pij[i][j];
                        Pij[i][j+1] += pe * old_Pij[i][j];
                        Pij[i-1][j] += pn * old_Pij[i][j];
                        toSouth[j] += ps * old_Pij[i][j];
                }

		// West buffer, storing data to be sent to west processor
		else if (j == 0) {
                        toWest[i] += pw * old_Pij[i][j];
                        Pij[i][j+1] += pe * old_Pij[i][j];
                        Pij[i-1][j] += pn * old_Pij[i][j];
                        Pij[i+1][j] += ps * old_Pij[i][j];
                }

		// East buffer, storing data to be sent to east processor
		else if (j == nj-1) {
                        Pij[i][j-1] += pw * old_Pij[i][j];
                        toEast[i] += pe * old_Pij[i][j];
                        Pij[i-1][j] += pn * old_Pij[i][j];
                        Pij[i+1][j] += ps * old_Pij[i][j];
                }

		// Inner part without ghost cells
		//else if (i > 0 && i < ni-1 && j > 0 && j < nj-1) 
		else {
                	Pij[i][j-1] += pw * old_Pij[i][j];
                	Pij[i][j+1] += pe * old_Pij[i][j];
                	Pij[i-1][j] += pn * old_Pij[i][j];
                	Pij[i+1][j] += ps * old_Pij[i][j];
		}
				
	}
  }


  // Exchanging boundary points (ghost cell)
  MPI_Send(toWest, ni, MPI_DOUBLE, westNB, 1, comm2D);
  MPI_Send(toEast, ni, MPI_DOUBLE, eastNB, 2, comm2D);
  MPI_Send(toNorth, nj, MPI_DOUBLE, northNB, 3, comm2D);
  MPI_Send(toSouth, nj, MPI_DOUBLE, southNB, 4, comm2D);
  
  // Receiving
  MPI_Recv(Dom->west, ni, MPI_DOUBLE, westNB, 2, comm2D, MPI_STATUS_IGNORE);
  MPI_Recv(Dom->east, ni, MPI_DOUBLE, eastNB, 1, comm2D, MPI_STATUS_IGNORE);
  MPI_Recv(Dom->north, nj, MPI_DOUBLE, northNB, 4, comm2D, MPI_STATUS_IGNORE);
  MPI_Recv(Dom->south, nj, MPI_DOUBLE, southNB, 3, comm2D, MPI_STATUS_IGNORE);

  //if (qi == 0 && qj == 1)
  //      printf(" point received from proc %d  is %f \n", eastNB, Dom->east[3]);

  // Handling data receiving from north & south processor
  for (j = 0; j < nj; j++) {
  	Pij[0][j] += Dom->north[j];
	Pij[ni-1][j] += Dom->south[j];
  }

  // Handling data receiving from east & west processor
  for (i = 0; i < ni; i++) {
        Pij[i][nj-1] += Dom->east[i];
        Pij[i][0] += Dom->west[i];
  }
  
  /*-------------------- free arrays */
  free(toNorth);
  free(toSouth);
  free(toEast);
  free(toWest);
  
}


double ComputErr(double** x, double** y, int ni, int nj, MPI_Comm comm2D){
  double max=0;
  double temp;
  int i,j;
  /*local computation of max{abs(u_k+1(i)-u_k(i))}*/
  
  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      temp=fabs(x[i][j]-y[i][j]);
      if(temp>max)
        max=temp;
    }
  }
/*--------------------allreduce the results of max at root node*/
  MPI_Allreduce(&max,&temp,1,MPI_DOUBLE,MPI_MAX,comm2D);
  return temp;
}

int normliz(double** x, int ni, int nj, MPI_Comm comm2D){
  double temp=0.0;
  int i,j;
  /*local computation of max{abs(u_k+1(i)-u_k(i))}*/
  for(i=0;i<ni;i++)
    for(j=0;j<nj;j++)
      temp+=x[i][j];
/*--------------------allreduce the results of max at root node*/
  MPI_Allreduce(&temp,&temp,1,MPI_DOUBLE,MPI_SUM,comm2D);
  if (temp == 0.0)
    return 1;
  for(i=0;i<ni;i++)
    for(j=0;j<nj;j++)
      x[i][j] /=temp;   
  return 0;
}


void get_prob(int i, int j, int ni, int nj, int pi, int pj,
		 int *dims, double *prob) {
/* get probability stencil */
/* ------------------------*/
/* i, j = local coordinates of point in each cell 
   ni, nj = dimensions of cell [assumed to be all the same size]
   pi, pj = process id in cartesian topology 
   dims   = array of size to containing the max number of 
            processes in each direction of cartesian mesh 
------------------------------------------------------------*/
  double tw, te, tn, ts;
  int ii, jj, niT, njT;
  /*-------------------- global indices */
  ii = pi*ni + i;
  jj = pj*nj + j;
  /*-------------------- global matrix dimensions */
  niT = dims[0]*ni;
  njT = dims[1]*nj;  
  // T1: pl = .15  pr = .3   == all multiplied by .25 times fun. of i,j
  // T2  p = p = .2
  te = 0.25* 0.30 * (double)(njT - jj+1)/ (double) njT;
  tw = 0.25* 0.15 * (double)(njT + jj-1)/ (double) njT;
  ts = 0.25* 0.20 * (double)(niT - ii+1)/ (double) niT;
  tn = 0.25* 0.20 * (double)(niT + ii-1)/ (double) niT;
  /*-------------------- save trans. probabilities */
  prob[1] = tw;      // west
  prob[2] = te;      // east
  prob[3] = tn;      // north 
  prob[4] = ts;      // south
  prob[0] = 1.0 - (te+tw+tn+ts);
}
