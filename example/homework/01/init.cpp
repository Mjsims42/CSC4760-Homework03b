class LinearDistribution
{
public:

  LinearDistribution(int _P, int _M) : the_P(_P), the_M(_M)
                                       {nominal = the_M/the_P;
					extra   = the_M%the_P;
					factor1 = extra*(nominal+1);}
  virtual ~LinearDistribution() {}
  
  void global_to_local(int I, int &p, int &i) const
              {p = (I < factor1) ? I/(nominal+1) : extra+((I-factor1)/nominal);
		i = I - ((p < extra) ? p*(nominal+1) : (factor1+(p-extra)*nominal));}
  int local_to_global(int p, int i) const
  {return i + ((p < extra) ? p*(nominal+1) : (factor1+(p-extra)*nominal));}
					      
  int m(int p) const {return (p < extra) ? (nominal+1) : nominal;}
  
  int M() const {return the_M;}
  int P() const {return the_P;}

protected:
  int the_M, the_P, nominal, extra, factor1;
};


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

	int P,Q,M;
	int array[3];

  	if(myrank == 0)
  	{
     		P = atoi(argv[1]); Q = atoi(argv[2]); M = atoi(argv[3]);
     		array[0] = P;
     		array[1] = Q;
		array[2] = M;
		printf("P: %d, Q: %d, M: %d\n", P, Q, M);
  	}
	
  	MPI_Bcast(array, 3, MPI_INT, 0, MPI_COMM_WORLD);

  	if(myrank != 0)
  	{
    		P = array[0];
    		Q = array[1];
		M = array[2];
  	}
	
	if(size != P * Q)
        {
                printf("Bad Size\n");
                MPI_Finalize();
                return 0;
        }
	char subcommunicator;
	
	int color1 = myrank / Q;

	int color2 = myrank % Q;


 	
      	
	MPI_Comm newComm;
	MPI_Comm_split(MPI_COMM_WORLD,color1,myrank,&newComm);
	
	MPI_Comm newComm1;
	MPI_Comm_split(MPI_COMM_WORLD,color2,myrank,&newComm1);
	// Get my rank in the new communicator
	int *X = new int[M];
	int *count = new int[Q];
	int *displs = new int[Q];
	int z0;
	int z1;
	int factor = ((M % P) * ((M/P) + 1));
	if(M % P != 0)
	{
		z0 = 1;
	}
	else
	{
		z0 = 0;
	}
	if((((color1 + 1) * ((M/P) + 1)) <= factor) && M % P != 0)
	{
		z1 = 0;
	}
	else
	{
		z1 = 1;
	}
	int X_Length = ((M/P) + z0 - z1);
	MPI_Request request;
	if(color2 == 0)
	{
		if(myrank == 0)
		{
			printf("I am process 0,0");
			for(int i = 0; i < M; i++)
			{
				X[i] = i + 1;
			}
			printf("for loop size of %d/%d for rank 0: \n", M,P);
        		for(int i = 0; i < M; i++)
        		{
                		printf("rank 0: %i \n",X[i]);
        		}
			printf("We are going to scatter rank: %d\n", myrank);
			for(int i = 0; i < Q; i++) 
        		{
                		count[i] = Q;
                		displs[i] = i * Q;
        		}
        		MPI_Scatterv(X,count,displs,MPI_INT,X,Q,MPI_INT,0, newComm1);	
		}
		else
		{	printf("\n\nrank %d size is %d\n\n",myrank,X_Length);	
			MPI_Scatterv(NULL, NULL, NULL, MPI_INT,X,Q, MPI_INT,0, newComm1);
		}
	delete[] count;
        delete[] displs;
	}

	if(color2 == 0)
	{
		printf("after Scatterv for rank: %d\n\n",myrank);
		for(int i =0; i < Q; i++)
        	{
                	printf("rank: %d count: %d value: %d\n",myrank,i,X[i]);
        	}
	}
	printf("lets Bcast horizontally rank: %d\n", myrank);
	MPI_Bcast(X, X_Length, MPI_INT, 0, newComm);
	for(int i = 0; i < (M/Q) + 1; i++)
        {
        printf("\nProcess(%d,%d) has at count: %d Value: %d\n\n",color1,color2,i,X[i]);
        }
	int *Y = new int[M/Q];
	int nominal1 = M/P;
	int nominal2 = M/Q;
	int extra1 = M%P;
	int extra2 = M%Q;
	for(int i = 0; i < (M/Q) + 1; i++)
	{
        	Y[i] = 0;
        }
	for(int i = 0; i < X_Length; i++)
	{
		int I = i + ((color1 < extra1) ? (nominal1+1)*color1 : 
				(extra1*(nominal1 + 1) + (color1-extra1)*nominal1));

		int qhat = (I < extra2*(nominal2+1)) ? I/(nominal2+1) : 
			(extra2 + (I-extra2*(nominal2 + 1))/nominal2);
		int jhat = I - ((qhat < extra2) ? (nominal2+1)*qhat :
				(extra2*(nominal2+1) + (qhat - extra2)*nominal2));
		if(qhat == color2)
		{
			Y[jhat] = X[i];
		}
	
	} 
	for(int i = 0; i < (M/Q) + 1; i++)
        {        
        printf("\nProcess(%d,%d) has at count: %d Value: %d\n\n",color1,color2,i,Y[i]);
        }
	
	int *Z = new int[M/Q];
	printf("lets reduceALL rank: %d\n", myrank);
	MPI_Allreduce(Y, Z, (M/Q) + 1, MPI_INT, MPI_SUM, newComm1);
	for(int i = 0; i < (M/Q) + 1; i++)
	{	
	printf("\nProcess %d has at count: %d Value: %d\n\n",myrank,i,Z[i]);
	}




	delete[] Y;
	delete[] X;

	MPI_Finalize();
}
