#include <mpi.h>
#include <math.h>
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
 #endif

  int check2DHeat(double** H, long n, long rank, long P, long k);

 #ifdef __cplusplus
}
#endif


double genH0(long r, long c, long n) {
  double val = (double)(c == (n/2));
  return val;
}



int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <N> <K>"<<std::endl;
    return -1;
  }
  MPI_Init(&argc,&argv);
  long n, K;
  n = atol(argv[1]);
  K = atol(argv[2]);

  int rank,np;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&np);

   int p = sqrt(np);
   long division= n/p; 
   
   int r_division = rank/p,c_division = rank%p;
   double**  H = new double*[division];
   double**  G = new double*[division];
   for(long i=0;i<division;i++)
     {
     H[i] = new double[division];
     G[i] = new double[division];
     
   }
   
   long r_start = (r_division*division),c_start = (c_division*division);
   long r_end = r_start+division,c_end = c_start+division;

  for (long r = r_start,rset=0; r<r_end; r++,rset++) {
    for (long c= c_start,cset=0; c<c_end; c++,cset++) {
       H[rset][cset] = genH0(r, c,n);
    }
  } 
  
  double *textp = new double[division];
  double *tm = new double[division];
  double *rightgt = new double[division];
  double *leftft = new double[division];
  double *tempp1 = new double[division];
  double *bmtl = new double[division];
  double *rightgtl = new double[division];
  double *leftftl = new double[division];
  
  int top,bottom,left,right;
  left =c_division?rank-1:-1;
  right = (c_division == (p-1))?-1:rank+1;
  top = rank-p;
  bottom = rank+p;
  MPI_Status sts[4];
  MPI_Request rqts[8];
  long count = 0;
  double start = MPI_Wtime();
  for (long it = 0; it<K; it++) 
  {
     for(long i =0,iteration = 0;i < division;i++)
       {
   leftft[iteration] = H[i][0];
   rightgt[iteration] = H[i][division-1];
   textp[iteration]  = H[0][i];
   tm[iteration] = H[division-1][i];
   leftftl[iteration] = H[i][0];
   rightgtl[iteration] = H[i][division-1];
   tempp1[iteration]  = H[0][i];
   bmtl[iteration] = H[division-1][i];
   iteration++;
       }
     count = 0;
     if(top>=0){
       MPI_Isend(tempp1,division,MPI_DOUBLE,top,0,MPI_COMM_WORLD,&rqts[count]);
       count++;
     }
     if(right != -1){
       MPI_Isend(rightgtl,division,MPI_DOUBLE,right,1,MPI_COMM_WORLD,&rqts[count]);
       count++;
     }
     if(bottom < np){
       MPI_Isend(bmtl,division,MPI_DOUBLE,bottom,2,MPI_COMM_WORLD,&rqts[count]);
       count++;
     }
     if(left != -1){
       MPI_Isend(leftftl,division,MPI_DOUBLE,left,3,MPI_COMM_WORLD,&rqts[count]);
       count++;
     }
     for(long i=1;(i+1)<division;i++)
       {
   for(long j=1;(j+1)<division;j++)
     {
       G[i][j] = (H[i][j] + H[i-1][j] + H[i+1][j] + H[i][j-1] + H[i][j+1])/5;
     }
       }
     if(top >= 0){
       MPI_Recv(textp,division,MPI_DOUBLE,top,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      
     }
     if(right != -1){
       MPI_Recv(rightgt,division,MPI_DOUBLE,right,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      
     }
     if(bottom < np){
       MPI_Recv(tm,division,MPI_DOUBLE,bottom,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       
     }
     if(left != -1){
       MPI_Recv(leftft,division,MPI_DOUBLE,left,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       
     }
     G[0][0] = (H[0][0] + textp[0] + leftft[0] + H[0][1] + H[1][0])/5;
     G[0][division-1] = (H[0][division-1] + textp[division-1] + rightgt[0]+H[0][division-2]+H[1][division-1])/5;
     G[division-1][0] = (H[division-1][0] + tm[0] + leftft[division-1]+H[division-1][1]+H[division-2][0])/5;
     G[division-1][division-1] = (H[division-1][division-1] + tm[division-1] +
          rightgt[division-1]+H[division-1][division-2]+H[division-2][division-1])/5;
     for(long i=1,j=1;i<division-1;i++,j++)
   {
      G[0][i] = (H[0][i]+H[1][i]+textp[i]+H[0][i-1]+H[0][i+1])/5;
      G[division-1][i] = (H[division-1][i]+H[division-2][i]+tm[i]+H[division-1][i-1]+H[division-1][i+1])/5;
      G[j][0] = (H[j][0]+leftft[j]+H[j][1]+H[j-1][0]+H[j+1][0])/5;
      G[j][division-1] = (H[j][division-1]+H[j-1][division-1]+H[j+1][division]+rightgt[j]+H[j][division-2])/5;
   }
      H = G;
      check2DHeat(H,n,rank,np,it);
      MPI_Waitall(count,rqts,sts);
             
   }
  if(rank == 0)
    {
       double end = MPI_Wtime();
       cerr<<end-start<<endl; 
    }
 for(long i=0;i<division;i++)
   delete[] H[i];
  delete[] H;
  delete[] textp;
  delete[] tm;
  delete[] rightgt;
  delete[] leftft;
  MPI_Finalize();

  return 0;
}
