#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#define VET_SIZE 100000

void bubbleSort(int n, int *vetor)
{
    int c = 0, d, troca, trocou = 1;
    while ((c < (n - 1)) & trocou)
    {
        trocou = 0;
        for (d = 0; d < n - c - 1; d++)
            if (vetor[d] > vetor[d + 1])
            {
                troca = vetor[d];
                vetor[d] = vetor[d + 1];
                vetor[d + 1] = troca;
                trocou = 1;
            }
        c++;
    }
}

int main(int argc, char **argv){
    MPI_Status status; // Message status
    double t1, t2; // Count exectuion time
    int my_rank; // Process ID
    int proc_n; // Number of process
    int aux; // Store auxiliary values

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);

    t1 = (my_rank == 0)? MPI_Wtime(): 0; // Time counter on process 0
    #define LAST (proc_n - 1) // Last process

    int i; // for counters
    unsigned char end = 0; // Control the main loop
    unsigned char k; // Process counter
    int psize = VET_SIZE/proc_n; // tamanho do vetor do processo
    int pTochange = VET_SIZE/(proc_n*3); // Size to change date with others
    int left ,right; // Auxiliar vector for change operation
    unsigned char vet_ctrl[proc_n]; // Control vector
    int *vector = (int *)malloc((pTochange + psize) * sizeof(int));
    left = (my_rank !=0) ? my_rank - 1:0 ; // Vizinho esq
    right = (my_rank < LAST) ? my_rank + 1: LAST; // vizinho dir

    memset(vet_ctrl, 0, sizeof(vet_ctrl));

    for (i = 0; i < psize; i++)
    {
        vector[i] = (proc_n - my_rank) * psize - i;
    }
    while(!end) {
        // ordenacao local
        bubbleSort(psize, vector);
        // manda pra direita se n for o ultmo
        if (my_rank != LAST){
            MPI_Send(&vector[psize - 1], 1, MPI_INT, right, 0, MPI_COMM_WORLD);
        }
        // recebe da esquerda se n for 0
        if (my_rank != 0){
            MPI_Recv(&aux, 1, MPI_INT, left, 0, MPI_COMM_WORLD, &status);
            if (aux < vector[0]){
                vet_ctrl[my_rank] = 1;
            }
            else{
                vet_ctrl[my_rank] = 0;
            }
        }
        if (my_rank == 1){
            MPI_Send(vector, 1, MPI_INT, left, 0, MPI_COMM_WORLD);
        }
        if (my_rank == 0){
            MPI_Recv(&aux, 1, MPI_INT, right, 0, MPI_COMM_WORLD, &status);
            if(aux > vector[psize-1] ){
                vet_ctrl[my_rank] = 1;
            }else{
                vet_ctrl[my_rank] = 0;
            }
        }
        for (i = 0; i < proc_n; i++)
        {
            MPI_Bcast(&vet_ctrl[i], 1, MPI_UNSIGNED_CHAR, i, MPI_COMM_WORLD);
        }
        k = 0;
        for (i = 0; i < (proc_n*2)-2; i++){
            if (vet_ctrl[i] == 1) {
                k++;
            }
        }
        if (k == proc_n){
            end = 1;
            break;
        }
        if (my_rank != 0) {
            MPI_Send(vector, pTochange, MPI_INT, left, 0, MPI_COMM_WORLD);
            // send to left
        }
        if (my_rank != LAST){
            MPI_Recv(&vector[psize], pTochange, MPI_INT, right, 0, MPI_COMM_WORLD, &status);
            bubbleSort(psize,vector + pTochange);
            MPI_Send(&vector[psize], pTochange, MPI_INT, right, 0, MPI_COMM_WORLD);
        }
        if (my_rank !=0){
            MPI_Recv(vector, pTochange, MPI_INT, left, 0, MPI_COMM_WORLD, &status);
        }
    }
    t2 = (my_rank == 0)?MPI_Wtime():0;
    if (my_rank == 0)
    {
        printf("Elapsed: %.4f s\n\n", t2 - t1);
    }
    free(vector);
    MPI_Finalize();
    return 0;
}