#include <mpi.h>
#include <iostream>
#include <stdio.h>
#define PING_PONG_LIMIT 10
using namespace std;
int main (int argc, char* argv[])
       {


          if (argc < 2) {
                std::cerr<<"usage: mpirun "<<argv[0]<<" <value>"<<std::endl;
                return -1;
                  }


                MPI_Init(NULL, NULL);


                int world_size;
                MPI_Comm_size(MPI_COMM_WORLD, &world_size);

                int world_rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

                if(world_size!=2)
                {
                        std::cerr<<"Please enter no of processors as 2";
                        return -1;
                }



                int count=0;

                        int processor1_rev;
                        

        int k;
        int Variable1=atoi(argv[1]);
        int recv;
    if (world_rank == 0)
                {
        MPI_Send(&Variable1, 1, MPI_INT, 1, 0,
                 MPI_COMM_WORLD);
                        MPI_Recv(&recv, 1, MPI_INT, 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                   printf("%d received variable %d from %d\n",
               world_rank, recv, 1);
                 }
 else
        {

        MPI_Recv(&k, 1, MPI_INT, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        k+=2;

MPI_Send(&k, 1, MPI_INT, 0, 0,MPI_COMM_WORLD);
           }
MPI_Finalize();
  return 0;
}

