/*  Shefali
 *  Sharma
 *  sharma92
 */

#ifndef A1_HPP
#define A1_HPP

#include <vector>
#include <mpi.h>
#include <unistd.h>

int connected_components(std::vector<signed char>& A, int n, int q, const char* out, MPI_Comm comm) {
    
    int localmax;
    MPI_Comm comm_new;
    
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    //std::vector<int> P(0);
    int b = n/q;
    //P.resize(q * q, 0);
    std::vector<int> M_Matrix(b*b,0);
    std::vector<int> M_helper(b*b,0);
    std::vector<int> P(b*b, 0);
    std::vector<int> R(b*b, 0);
    std::vector<int> T(b*b, 0);
    std::vector<int> P_Prime(b*b, 0);
    std::vector<int> MaxPT(b*b, 0);
    //M_Matrix.resize(q*q,0);
    
    //printf("q = %d \n", q);
    //int ne = ((2.0 * q * b * b) / (n * n - n));
    int col = rank % q;
    int row = rank / q;
    
    //if (col == row) for (int i = 0; i < b; ++i) A[i * b + i] = 1;
    //if(col<=row){
    for(int i = 0; i<b; i++){
        for(int j = 0; j<b; j++){
            if(A[i * b + j] == 1){
                //printf("A[ %d ] = %d \n", (i * b + j), A[i * b + j]);
                //printf("   row = %d, col = %d \n", row, col);
                P[i * b + j] = i + (row*b);
            }
            else{
                P[i * b + j] = -1;
            }
        }
    }
    
    //}
    
    //MPI_Barrier(comm);
    /*    sleep(rank);
     printf("   row = %d, col = %d \n", row, col);
     
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(P[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     
     printf("\n");
     */
    //    printf("..\n");
    MPI_Barrier(comm);
    
    //printf("Calculating local max . . .\n");
    for(int i = 0; i < b; ++i){
        for(int j = 0; j < b; j++){
            if(M_Matrix[0 + j] <= P[i * b + j]){
                M_Matrix[0 + j] = P[i * b + j];
                //M_Matrix[i * b + j] = 0;
                
            }
            //else{
            //    M_Matrix[i * b + j] = 0;
            //}
        }
        
    }
    for(int i = 0; i < b; i++){
        for(int j = 0; j < b; j++){
            M_Matrix[i * b + j] = M_Matrix[0 + j];
        }
    }
    
    /*    sleep(rank);
     printf("   row = %d, col = %d \n", row, col);
     
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(M_Matrix[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     
     printf("\n");
     */
    MPI_Barrier(comm);
    
    MPI_Comm_split(comm, col,rank, &comm_new);
    MPI_Allreduce(M_Matrix.data(), P.data(), (b * b), MPI_INT, MPI_MAX, comm_new);
    MPI_Comm_free(&comm_new);
    MPI_Barrier(comm);
    
    sleep(rank);
    printf("\n P = \n");
    for (int i = 0; i < b; ++i) {
        
        for (int j = 0; j < b; ++j) {
            std::cout << static_cast<int>(P[i * b + j]) << " ";
        }
        std::cout << std::endl;
    }
    
    MPI_Barrier(comm);
    
    for (int i = 0; i < b; ++i) {
        
        for (int j = 0; j < b; ++j) {
            if(static_cast<int>(A[i * b + j])==1)
            {
                M_helper[i * b + j] = P[i * b + j];
            }
            else
            {
                M_helper[i * b + j] = 0;
            }
        }
    }
    
    /*    sleep(rank);
     printf("\n M_Helper = \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(M_helper[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     */
    MPI_Barrier(comm);
    for (int i = 0; i < b; ++i) {
        localmax=0;
        for (int j = 0; j < b; ++j){
            if(localmax < M_helper[i * b + j])
            {
                localmax = M_helper[ i * b + j];
            }
        }
        for (int j = 0; j < b; ++j)
        {
            M_helper[i * b + j] = localmax;
        }
        
    }
    
    /*    sleep(rank);
     printf("\n M_helper after local max=  \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(M_helper[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     */
    MPI_Barrier(comm);
    
    MPI_Comm_split(comm, row, rank, &comm_new);
    MPI_Allreduce(M_helper.data(), R.data(), (b * b), MPI_INT, MPI_MAX, comm_new);
    MPI_Comm_free(&comm_new);
    MPI_Barrier(comm);
    
    /*    sleep(rank);
     printf("\nR = \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(R[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     */
    MPI_Barrier(comm);
    
    for (int i = 0; i < b; ++i) {
        for (int j = 0; j < b; ++j) {
            if(R[i * b + j] == j + b*col){
                M_Matrix[i * b + j] = P[i * b + j];
            }
            else{
                M_Matrix[i * b + j] = 0;
            }
        }
    }
    MPI_Barrier(comm);
    /*    sleep(rank);
     //M_Matrix.resize(b*b, 0);
     printf("\nM_Matrix = \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(M_Matrix[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     MPI_Barrier(comm);
     */
    for (int i = 0; i < b; ++i) {
        localmax=0;
        for (int j = 0; j < b; ++j){
            if(localmax < M_Matrix[i * b + j])
            {
                localmax = M_Matrix[ i * b + j];
            }
        }
        for (int j = 0; j < b; ++j)
        {
            M_Matrix[i * b + j] = localmax;
        }
        
    }
    /*
     printf("\nM_Matrix After local max= \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(M_Matrix[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     */
    MPI_Barrier(comm);
    
    MPI_Comm_split(comm, row, rank, &comm_new);
    MPI_Allreduce(M_Matrix.data(), P_Prime.data(), (b * b), MPI_INT, MPI_MAX, comm_new);
    MPI_Comm_free(&comm_new);
    MPI_Barrier(comm);
    
    sleep(rank);
    printf("\nP_Prime = \n");
    for (int i = 0; i < b; ++i) {
        
        for (int j = 0; j < b; ++j) {
            std::cout << static_cast<int>(P_Prime[i * b + j]) << " ";
        }
        std::cout << std::endl;
    }
    
    for (int i = 0; i < b; ++i) {
        for (int j = 0; j < b; ++j) {
            if(P[i * b + j] == i + b*row){
                M_Matrix[i * b + j] = P_Prime[i * b + j];
            }
            else{
                M_Matrix[i * b + j] = 0;
            }
        }
    }
    
    MPI_Barrier(comm);
    sleep(rank);
    //M_Matrix.resize(b*b, 0);
    /*printf("\nM_Matrix = \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(M_Matrix[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     */
    MPI_Barrier(comm);
    
    for (int i = 0; i < b; ++i) {
        localmax=0;
        for (int j = 0; j < b; ++j){
            if(localmax < M_Matrix[i * b + j])
            {
                localmax = M_Matrix[ i * b + j];
            }
        }
        for (int j = 0; j < b; ++j)
        {
            M_Matrix[i * b + j] = localmax;
        }
        
    }
    
    /*    printf("\nM_Matrix After local max= \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(M_Matrix[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     */
    MPI_Barrier(comm);
    
    MPI_Comm_split(comm, row, rank, &comm_new);
    MPI_Allreduce(M_Matrix.data(), T.data(), (b * b), MPI_INT, MPI_MAX, comm_new);
    MPI_Comm_free(&comm_new);
    MPI_Barrier(comm);
    
    /*    sleep(rank);
     printf("\nT = \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(T[i * b + j]) << " ";
     }
     std::cout << std::endl;
     }
     */
    for (int i = 0; i < b; ++i) {
        
        for (int j = 0; j < b; ++j) {
            if(T[i * b + j]<P_Prime[i * b + j]){
                MaxPT[i * b + j] = P_Prime[i * b + j];
            }
            else{
                MaxPT[i * b + j] = T[i * b + j];
            }
        }
    }
    
    sleep(rank);
    printf("\nMaxPT = \n");
    for (int i = 0; i < b; ++i) {
        
        for (int j = 0; j < b; ++j) {
            std::cout << static_cast<int>(MaxPT[i * b + j]) << " ";
        }
        std::cout << std::endl;
    }
    
    
    //-------------------------------------------------------------------------
    /*
     
     std::vector<int> temp(b*b, 0);
     std::vector<int> temp_helper(b*b, 0);
     std::vector<int> temp_M_Matrix(b*b, 0);
     
     sleep(rank);
     printf("\n P = \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(P[i * b + j]) << " ";
     temp[i * b + j] = P[i * b + j];
     }
     std::cout << std::endl;
     }
     
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     if(static_cast<int>(A[i * b + j])==1)
     {
     temp_helper[i * b + j] = P[i * b + j];
     }
     else
     {
     temp_helper[i * b + j] = 0;
     }
     }
     }
     
     for (int i = 0; i < b; ++i) {
     localmax=0;
     for (int j = 0; j < b; ++j)
     {
     if(localmax<temp_helper[j*b+i]){
     localmax = temp_helper[j*b+i];
     }
     }
     for (int j = 0; j < b; ++j)
     {
     temp_helper[j*b+i] = localmax;
     }
     
     }
     
     MPI_Barrier(comm);
     
     MPI_Comm_split(comm, col, rank, &comm_new);
     MPI_Allreduce(temp_helper.data(), temp_M_Matrix.data(), (b * b), MPI_INT, MPI_MAX, comm_new);
     MPI_Comm_free(&comm_new);
     MPI_Barrier(comm);
     
     sleep(rank);
     
     printf("\n P After row Allreduce= \n");
     for (int i = 0; i < b; ++i) {
     
     for (int j = 0; j < b; ++j) {
     std::cout << static_cast<int>(temp_M_Matrix[i * b + j]) << " ";
     //temp[i * b + j] = P[i * b + j];
     }
     std::cout << std::endl;
     }
     
     */
    return -1;
} // connected_components

#endif // A1_HPP

