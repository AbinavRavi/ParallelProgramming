#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "heat.h"
#include "helper.h"

void init_parallel(int iproc,
               int jproc,
               int imax,
               int jmax,
               int *myrank,
               int *il,
               int *ir,
               int *jb,
               int *jt,
               int *l_rank,
               int *r_rank,
               int *b_rank,
               int *t_rank,
               int *omg_i,
               int *omg_j,
               int num_proc)
{

  *omg_i = (*myrank % iproc) + 1;
  *omg_j = ((*myrank+1-*omg_i)/iproc)+1;

  *il = (*omg_i-1)*(imax/iproc) + 1;
  *ir = (*omg_i!=iproc)?((*omg_i)*(imax/iproc)):imax;

  *jb = (*omg_j-1)*(jmax/jproc) + 1;
  *jt = (*omg_j!=jproc)?((*omg_j)*(jmax/jproc)):jmax;

  if(*il == 1)      *l_rank = MPI_PROC_NULL;
  else              *l_rank = *myrank - 1;

  if(*ir == imax)   *r_rank = MPI_PROC_NULL;
  else              *r_rank = *myrank + 1;


  if(*jb == 1)      *b_rank = MPI_PROC_NULL;
  else              *b_rank = *myrank - iproc;

  if(*jt == jmax)   *t_rank = MPI_PROC_NULL;
  else              *t_rank = *myrank + iproc;


//printf("Thread_id: %d omg_ij: %d%d \nil: %d, ir: %d, jb: %d, jt: %d \n",*myrank,*omg_i,*omg_j, *il,*ir,*jb,*jt);

//printf("l_rank: %d, r_rank: %d, b_rank: %d, t_rank: %d \n \n", *l_rank,*r_rank,*b_rank,*t_rank);
}

void pressure_comm(double *P, int il, int ir,int jb, int jt, int n,
          				int l_rank,int r_rank,int b_rank, int t_rank,
									double *bufSend,double *bufRecv, MPI_Request *request1, 
									MPI_Request *request2, MPI_Status status, int chunk, MPI_Comm comm)

{
	int x_dim = ir - il + 1;
	int y_dim = jt - jb + 1;

  //Send to left &  recieve from right
	if (l_rank != MPI_PROC_NULL)
	{
		for(int j=1; j<=y_dim; j++)
		{
			bufSend[j-1] = P[map(1, j ,n+2)];
		}
		MPI_Send( bufSend, y_dim, MPI_DOUBLE, l_rank, 1, comm );
	}
	if (r_rank != MPI_PROC_NULL) // Receive from right
	{
		MPI_Recv( bufRecv, y_dim, MPI_DOUBLE, r_rank, 1, comm, &status );
		for(int j=1; j<=y_dim; j++)
		{
			P[map(x_dim+1, j ,n+2)] = bufRecv[j-1];
		}
  }
	//MPI_Wait(request1, &status);
	//MPI_Wait(request2, &status);
  // send to right & recieve from left
  if (r_rank != MPI_PROC_NULL)
  {
		for(int j=1; j<=y_dim; j++)
		{
			bufSend[j-1] = P[map(x_dim, j ,n+2)];
		}
		MPI_Send( bufSend, y_dim, MPI_DOUBLE, r_rank, 1, comm);
	}
	if (l_rank != MPI_PROC_NULL)
	{
		MPI_Recv(bufRecv, y_dim, MPI_DOUBLE, l_rank, 1, comm, &status);
		for(int j=1; j<=y_dim; j++)
		{
			P[map(0, j ,n+2)] = bufRecv[j-1];
		}
	}
	//MPI_Wait(request1, &status);
	//MPI_Wait(request2, &status);

  ///send to bottom recieve from top
	if (t_rank != MPI_PROC_NULL)
	{
		for(int i=1; i<=x_dim; i++)
		{
			bufSend[i-1] = P[map(i, y_dim ,n+2)];
		}
		MPI_Send( bufSend, x_dim, MPI_DOUBLE, t_rank, 1, comm );
	}
	if (b_rank != MPI_PROC_NULL)
	{
		MPI_Recv( bufRecv, x_dim, MPI_DOUBLE, b_rank, 1, comm, &status);
		for(int i=1; i<=x_dim; i++)
		{
			P[map(i, 0 ,n+2)] =bufRecv[i-1];
		}
  }
	//MPI_Wait(request1, &status);
	//MPI_Wait(request2, &status);

    ///send to bottom recieve from top
  if (b_rank != MPI_PROC_NULL)
	{
		for(int i=1; i<=x_dim; i++)
		{
			bufSend[i-1] = P[map(i, 1 ,n+2)];
		}
		MPI_Send( bufSend, x_dim, MPI_DOUBLE, b_rank, 1, comm );
	}

	if (t_rank != MPI_PROC_NULL) // Receive from the top
	{
		MPI_Recv( bufRecv, x_dim, MPI_DOUBLE, t_rank, 1, comm, &status );
		for(int i=1; i<=x_dim; i++)
		{
			P[map(i, y_dim+1 ,n+2)] = bufRecv[i-1];
		}
	}
	//MPI_Wait(request1, &status);
	//MPI_Wait(request2, &status);

}


double jacobi(double *h_new, double *h_old, int niters, int energy_intensity, int n, int iter_energy,  const int nsources, int sources[nsources][2], int rank, int size, int px, int py, MPI_Comm comm, int output)
{
    h_old = (double *)calloc(1, (n + 2) * (n + 2) * sizeof(double)); // extended with halos of width 1
    h_new = (double *)calloc(1, (n + 2) * (n + 2) * sizeof(double)); // extended with halos of width 1
    double *tmp;
		
		int il, ir, jb, jt, l_rank, r_rank, b_rank, t_rank, omg_i, omg_j, chunk;
		
		init_parallel(px, py, n+2, n+2, &rank, &il, &ir, &jb, &jt, 
									&l_rank, &r_rank, &b_rank, &t_rank, &omg_i, &omg_j, size);

		MPI_Request *request1 = NULL;
		MPI_Request *request2 = NULL;
		MPI_Status status;

		int x_dim = ir - il + 1;
		int y_dim = jt - jb + 1;
 		chunk = 0;

		double *bufSend = (double *)calloc(0, (n+2)* sizeof(double));
		double *bufRecv = (double *)calloc(0, (n+2)* sizeof(double));

		pressure_comm(h_old, il, ir, jb, jt, n, l_rank, r_rank, b_rank, t_rank, 
									bufSend, bufRecv, request1, request2, status, chunk, comm);
    
		for (int iter = 0; iter < niters; ++iter)
    {
      for (int j = 1; j < y_dim+1; ++j)
      {
        for (int i = 1; i < x_dim+1; ++i)
        {
          h_new[map(i, j ,n+2)] = h_old[map(i, j, n+2)] / 2.0 + (h_old[map(i - 1, j, n+2)] + h_old[map(i + 1, j, n+2)] + h_old[map(i, j - 1, n+2)] + h_old[map(i, j + 1, n+2)]) / 4.0 / 2.0;
        }
      }
      if (iter < iter_energy)
      {
        for (int i = 0; i < nsources; ++i)
        {
            h_new[map(sources[i][0], sources[i][1], n+2)] += energy_intensity; // heat rate
        }
      }

      tmp = h_new; // swap arrays
      h_new = h_old;
      h_old = tmp;
    	
		}
    if (output) printarr(h_new, n, rank);

    return calculate_total_heat(h_new, n);
}
