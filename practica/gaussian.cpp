// 
// This program is an implementation of parallel gaussian elimination.
//

#include <iostream>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <ctime>
#include <vector>
#include <iomanip>
#include <sys/time.h>

using namespace std;

// Sorts the input row into chunks to be scattered two all the processors.
void sortByProcess(vector<double> list1, double* list2, int count, int size);

// Swaps two rows.
void swap(double** list, int count, int row1, int row2);

void testProcesses(int coreNum, int argc, char * argv[]);

int rank_1; // size;

//Convertimos el main en un proceso
int runProcess(int size, int argc, char * argv[]){
  double sTime, eTime, rTime;
  ifstream inFile;
  int num_rows = 3200;
  int num_cols = 3200;
  int cur_control = 0;
  double * send_buffer = NULL;
  double * recv_buffer = NULL;
  double ** data = NULL;
  double determinant;
  vector<double> file_buffer;

  struct timeval timeBegin, timeEnd;
  long int totalTime;

  // If the input file is not given, print message and exit.
  if(argc < 2)
  {
    cout << "No input file given." << endl;
    MPI_Finalize();
    return 0;
  }
  // If the root node (0), then open the input file and read in the
  // number of rows.
  if(!rank_1)
  {
    inFile.open(argv[1]);
    inFile >> num_rows;
    file_buffer.resize(num_rows);
  }

  send_buffer = new double[num_rows];

  // Broadcasts the number of rows to each processor.
  MPI_Bcast (&num_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  num_cols = num_rows / size;
  
  // Allocate the memory on each processor.
  data = new double*[num_cols];
  for(int i = 0; i < num_cols; i++)
    data[i] = new double[num_rows];
  for(int i = 0; i < num_cols; i++)
  {
    for(int j = 0; j < num_rows; j++)
      data[i][j] = 0;
  }
  recv_buffer = new double[num_cols];
  
  // Scatter the data.
  for(int i = 0; i < num_rows; i++)
  {
    if(!rank_1)
    {
      for(int j = 0; j < num_rows; j++) {
        inFile >> file_buffer[j];
	//cout << setw(2) << file_buffer[j] << "\t";
      }
      //cout << endl;
      sortByProcess(file_buffer, send_buffer, num_rows, size);
    }
    // Scatters the data so that each process gets the next value for their columns.
    MPI_Scatter(send_buffer, num_cols, MPI_DOUBLE, recv_buffer, num_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(int j = 0; j < num_cols; j++)
    {
      data[j][i] = recv_buffer[j];
    }
  }
  delete [] recv_buffer;
  delete [] send_buffer;
  // Begin timing.
  MPI_Barrier(MPI_COMM_WORLD);
  //sTime = MPI_Wtime();
  
  gettimeofday(&timeBegin, 0);

  send_buffer = new double[num_rows];
  int cur_row = 0;
  int swaps = 0;
  double det_val = 1;
  int cur_index = 0;
  for(int i = 0; i < num_rows; i++)
  { 
    // Find the row to swap with.
    int rowSwap;
    if(cur_control == rank_1)
    {
      rowSwap = cur_row;
      double max = data[cur_index][cur_row];
      // Find the row to swap with.
      for(int j = cur_row + 1; j < num_rows; j++)
      {
        if(data[cur_index][j] > max)
        {
          rowSwap = j;
          max = data[cur_index][j];
        }
      }
    }
   
    // Find out if you need to swap and then act accordingly.
    MPI_Bcast(&rowSwap, 1, MPI_INT, cur_control, MPI_COMM_WORLD);
    if(rowSwap != cur_row)
    {
      swap(data, num_cols, cur_row, rowSwap);
      swaps++;
    }
    
    if(cur_control == rank_1)
    {
      // Generate the coefficients.
      for(int j = cur_row; j < num_rows; j++)
        send_buffer[j] = data[cur_index][j] / data[cur_index][cur_row];
    }
    // Send and recv the coefficients.
    MPI_Bcast(send_buffer, num_rows, MPI_DOUBLE, cur_control, MPI_COMM_WORLD);
    // Apply the coefficients to the data.
    for(int j = 0; j < num_cols; j++)
    {
      for(int k = cur_row + 1; k < num_rows; k++)
      {
        data[j][k] -= data[j][cur_row] * send_buffer[k];
      }
    }
    
    // Update the determinant value.
    if(cur_control == rank_1)
    {
      det_val = det_val * data[cur_index][cur_row];
      cur_index++;
    }
      
    // Increment the row that we are looking at
    // and increment the counter that tells each process where
    // to recv from. The counter resets to zero to give us a 
    // "round robin" communication pattern. Probably not very efficient,
    // but it will do for now.
    cur_control++;
    if(cur_control == size)
      cur_control = 0;
    cur_row++;
  }
  
  // Reduce all the determinant values from each process
  // with a multiplication operation.
  // Personally I really like the method I used to find the determinant:
  //   1. Each process just keeps multiplying the pivot value into the product.
  //   2. The reduce does a multiply on all of the individual products.
  // So there really is no extra work to find the determinant.
  MPI_Reduce(&det_val, &determinant, 1, MPI_DOUBLE, MPI_PROD, 0, MPI_COMM_WORLD);  
  // If we did an odd number of row swaps, negate the determinant.
  if(swaps % 2)
    determinant = -determinant;
  
  // End timing.
  MPI_Barrier(MPI_COMM_WORLD);
  //eTime = MPI_Wtime();
  //rTime = eTime - sTime;
  //cout << endl; 
  // Print matrix
  /*for (int i = 0; i < num_rows; ++i) {
	  for (int j = 0; j < num_rows; ++j) {
		  cout << data[i][j] << "\t";
	  }
	  cout << endl;
  }*/


    //Implementación nuetra de medir el tiempo
    gettimeofday(&timeEnd, 0);
    totalTime = timeEnd.tv_usec - timeBegin.tv_usec;
    totalTime += timeEnd.tv_sec - timeBegin.tv_sec;
    cout <<"Número de cores: " << rank_1 << ". Tiempo medido: "<<((long double) totalTime )/ 1000000.0 << endl;

  // If root node, output the runtime.
  if(!rank_1)
  {
    //cout << "Run Time: " << rTime << endl;
    cout << "Determinant value: " << determinant << endl;
  }
  
  // A bit of house cleaning.
  delete [] send_buffer;
  for(int i = 0; i < num_cols; i++)
    delete [] data[i];
  delete [] data;
  
  // Finalize and exit. 
  //MPI_Finalize();
} 

int main(int argc, char * argv[]){
    testProcesses(8, argc, argv);
}


void testProcesses(int coreNum, int argc, char * argv[]){
	MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank_1);
  	MPI_Comm_size(MPI_COMM_WORLD, &coreNum);

//for(int i = 1; i <= coreNum; i = 2*i){

        runProcess(coreNum, argc, argv);
   // }
	MPI_Finalize();
    //cout << "Nro de core: " << rank_1<< endl;
}


void sortByProcess(vector<double> list2, double* list1, int count, int size)
{
  int index = 0;
  for(int i = 0; i < size; i++)
  {
    for(int j = i; j < count; j += size)
    {
      list1[index] = list2[j];
      index++;
    }
  }
  return;
}

void swap(double** list, int count, int row1, int row2)
{
  double temp;
  if(row1 == row2)
    return;
  for(int i = 0; i < count; i++)
  {
    temp = list[i][row1];
    list[i][row1] = list[i][row2];
    list[i][row2] = temp;
  }
  return;
}

