using namespace std;

#include <iostream>

#include <assert.h>

#include <string>



#include <mpi.h>

int main(int argc, char **argv)
{
	int size, myrank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int P,Q;
	int array[2];
	
  	if(myrank == 0)
  	{

     		P = atoi(argv[1]); Q = atoi(argv[2]);
     		array[0] = P;
     		array[1] = Q;
  	}

  	MPI_Bcast(array, 2, MPI_INT, 0, MPI_COMM_WORLD);

  	if(myrank != 0)
  	{
    		P = array[0];
    		Q = array[1];
  	}
	if(size != P * Q)

        {
                printf("Bad Size\n");
                MPI_Finalize();
                return 0;
        }
	char subcommunicator;
	int color;
    	int key;
   
	int color1 = myrank / Q;
        int color2 = myrank % Q;
        MPI_Comm newComm;
        MPI_Comm_split(MPI_COMM_WORLD,color1,myrank,&newComm);
        MPI_Comm newComm1;
        MPI_Comm_split(MPI_COMM_WORLD,color2,myrank,&newComm1);


	printf("[MPI process %d] I am grid (%d,%d) and am ready to grid out.\n",myrank,color1,color2);
	

	

	MPI_Finalize();


}
