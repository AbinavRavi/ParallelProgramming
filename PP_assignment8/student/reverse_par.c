#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "helper.h"

void reverse(char *str, int strlen)
{
        // parallelize this function and make sure to call reverse_str()
        // on each processor to reverse the substring.
//printf("Before reverse");
        int np, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

int chunk;
int extra_chunk;
chunk=strlen/(np);
extra_chunk=strlen%np;



    if (rank == 0)
    {int x=0;
	for(int i=0;i<np;i++)
	{
if(extra_chunk==0){	MPI_Send(&str[(x)*(chunk+1)+(i-x)*chunk], chunk, MPI_CHAR , i, 0,MPI_COMM_WORLD);}

else{	MPI_Send(&str[i*(chunk+1)], chunk+1, MPI_CHAR , i, 0,MPI_COMM_WORLD);extra_chunk=extra_chunk-1;x=i+1;}

         }
     }

extra_chunk=strlen%np;

if(0<=rank && rank<extra_chunk)
{

char buffer[chunk+1];
MPI_Recv(&buffer, chunk+1, MPI_CHAR ,0, 0,MPI_COMM_WORLD , MPI_STATUS_IGNORE);

reverse_str(buffer, chunk+1);
MPI_Send(&buffer, chunk+1, MPI_CHAR , 0, 0,MPI_COMM_WORLD);


}
else
{
char buffer[chunk];
MPI_Recv(&buffer,chunk , MPI_CHAR ,0, 0,MPI_COMM_WORLD , MPI_STATUS_IGNORE);

reverse_str(buffer, chunk);
MPI_Send(&buffer, chunk, MPI_CHAR , 0, 0,MPI_COMM_WORLD);


}
if(rank==0)
{
	for(int i=0;i<np-extra_chunk;i++)
	{
	MPI_Recv(&str[chunk*i], chunk, MPI_CHAR ,np-1-i, 0,MPI_COMM_WORLD , MPI_STATUS_IGNORE);
	}
	for(int j=0;j<extra_chunk;j++)
	{
	MPI_Recv(&str[chunk*(np-extra_chunk)+j*(chunk+1)], chunk+1, MPI_CHAR ,extra_chunk-1-j, 0,MPI_COMM_WORLD , MPI_STATUS_IGNORE);
	}

}



  }
