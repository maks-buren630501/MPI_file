#include <iostream>
#include<string>
#include<string.h>
#include<mpi.h>
#include<chrono>
#include<fstream>

using namespace std;


int* getMatrix(int n, int m)
{
    return new int[n*m];
}

void showMatrix(int *matrix, int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            cout<<matrix[i*m +j]<<" ";
        }
        cout<<endl;
    }
}

void initMatrix(int *matrix,int m, int n)
{
    for(int i = 0; i < m*n; i++)
    {
        matrix[i] = std::rand() %20 +1;
    }
}

void initMatrixZero(int *matrix,int m, int n)
{
    for(int i = 0; i < m*n; i++)
    {
        matrix[i] = 0;
    }
}

void multMatrix2x2(int *matrixA, int *matrixB, int *matrixC)
{
    for(int i = 0; i < 2; i++)
    {
        for(int j =0; j < 2; j++)
        {
            for(int k =0; k < 2; k++)
            {
                matrixC[i*2 + j] += matrixA[i*2 + k] * matrixB[k*2 + j];
            }
        }
    }
}

void multMatrix(int *matrixA, int *matrixB, int *matrixC, int m1, int n1,int m2, int n2)
{
    for(int i = 0; i < m1; i++)
    {
        for(int j = 0; j < n2; j++)
        {
            for(int k =0; k < n1; k++)
            {
                matrixC[i*n2 + j] += matrixA[i*n1 + k] * matrixB[k*n2 + j];
            }
        }
    }
}

void copyMatrix(int *matrixA, int* matrixB, int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            matrixA[i*n +j] = matrixB[i*n +j];
        }
    }
}

void addMat(int *matR,int *matC,int offset,int size)
{
    for(int i = 0; i < size; i++)
    {
        matR[offset + i] = matC[i];
    }
}

bool compareMatrix(int *matA,int *matB,int m,int n)
{
    int fl = 0;
    for(int i = 0; i < m*n; i++)
    {
            if(matA[i] != matB[i])
            {
                fl = 1;
            }
    }
    if(fl == 0)
    {
        return true;
    }
    else
    {
        return false;
    }

}

void matrixToFile(const char* filePath, int *matrix, int m, int n)
{
    fstream file;
    file.open(filePath,std::fstream::out|std::fstream::binary);
    file.write(reinterpret_cast<const char *>(matrix),m*n*sizeof(int));
    file.close();
}

void matrixFromFile(const char* filePath, int *matrix, int m, int n)
{
    fstream file;
    file.open(filePath,std::fstream::in|std::fstream::binary);
    file.read((char*)matrix,m*n*sizeof(int));
    file.close();
}



int main(int argc, char* argv[])
{

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int m = 100;
    int n = 100;
    int r = m/(size-1);

    MPI_File aFile;
    MPI_File bFile;
    MPI_File cFile;


    if(rank == 0)
    {
        MPI_Status status;
        int *matA = getMatrix(m,n);
        initMatrix(matA,m,n);
        matrixToFile("matA.bin",matA,m,n);
        cout<<"A matrix"<<endl;
        showMatrix(matA,m,n);
        cout<<endl;

        int *matB = getMatrix(n,m);
        initMatrix(matB,n,m);
        matrixToFile("matB.bin",matB,n,m);
        cout<<"B matrix"<<endl;
        showMatrix(matB,n,m);
        cout<<endl;

        fstream file;
        file.open("matC.bin",std::fstream::out|std::fstream::binary);
        file.close();

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_File_open(MPI_COMM_WORLD,"matA.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&aFile);
        MPI_File_open(MPI_COMM_WORLD,"matB.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&bFile);
        MPI_File_open(MPI_COMM_WORLD,"matC.bin",MPI_MODE_RDWR,MPI_INFO_NULL,&cFile);

        MPI_File_close(&aFile);
        MPI_File_close(&bFile);
        MPI_File_seek(cFile,0,MPI_SEEK_SET);

        int *matC = new int[m*m];
        initMatrixZero(matC,m,m);
        MPI_File_read(cFile,matC,m*m,MPI_INT,&status);
        cout<<"final matrix"<<endl;
        showMatrix(matC,m,m);
        MPI_File_close(&cFile);

        //std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
        //std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
        //cout <<"time MPI = "<<time_span2.count()<<"seconds"<<endl;

    }
    else
    {
        MPI_Status status;
        MPI_Barrier(MPI_COMM_WORLD);
        int *matA = new int[r*n];
        int *matB = new  int[n*m];

        MPI_File_open(MPI_COMM_WORLD,"matA.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&aFile);
        MPI_File_open(MPI_COMM_WORLD,"matB.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&bFile);

        MPI_File_seek(aFile,(rank-1)*r*n*sizeof(int),MPI_SEEK_SET);
        MPI_File_read(aFile,matA,r*n,MPI_INT,&status);
        MPI_File_read(bFile,matB,n*m,MPI_INT,&status);
        int *matC = new int[r*m];
        initMatrixZero(matC,r,m);
        multMatrix(matA,matB,matC,r,n,n,m);



        MPI_File_open(MPI_COMM_WORLD,"matC.bin",MPI_MODE_RDWR,MPI_INFO_NULL,&cFile);
        MPI_File_seek(cFile,(rank-1)*r*m*sizeof(int),MPI_SEEK_SET);
        MPI_File_write(cFile,matC,r*m,MPI_INT,&status);

        MPI_File_close(&aFile);
        MPI_File_close(&bFile);
        MPI_File_close(&cFile);



    }

    MPI_Finalize();


    return 0;
}
