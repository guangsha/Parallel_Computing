////////////////////////////////////////////////////
// Cyclic Jacobi Algorithm for diagonalization 
// of symmetric matrices (parallel)                         
////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "omp.h"

using namespace std;

struct RotationParameter
{
  bool applied;
  int p;
  int q;
  double SinTheta;
  double CosTheta;
};

void RotateMore( double **matrix, int size, int index, RotationParameter *RotationParameterList )
{
  int p, q;
  double SinTheta, CosTheta;
  double *ColumnP = new double[size];
    
  p = RotationParameterList[index].p;
  q = RotationParameterList[index].q;
  SinTheta = RotationParameterList[index].SinTheta;
  CosTheta = RotationParameterList[index].CosTheta;
    
  for ( int i = 0; i < size; i++ )
    {
      if ( ( i != p ) && ( i != q ) )
        {
	  ColumnP[i] = CosTheta * matrix[i][p] + SinTheta * matrix[i][q];
	  matrix[i][q] = CosTheta * matrix[i][q] - SinTheta * matrix[i][p];
        }
    }
    
  for ( int i = 0; i < size; i++ )
    {
      if ( ( i != p ) && ( i != q ) )
        {
	  matrix[i][p] = ColumnP[i];
        }
    }
    
  delete [] ColumnP;
}

void rotate( int p, int q, double **matrix, int size, double epsilon, int index,
	     RotationParameter *RotationParameterList )
{
  double theta, SinTheta, CosTheta;
  double *rowp = new double[size];
  double *rowq = new double[size];
    
  // Compute the orthogonal matrix for Jacobi rotation
  theta = atan( 2 * matrix[p][q] / ( matrix[p][p] - matrix[q][q] ) ) / 2;
  RotationParameterList[index].SinTheta = SinTheta = sin(theta);
  RotationParameterList[index].CosTheta = CosTheta = cos(theta);
  RotationParameterList[index].p = p;
  RotationParameterList[index].q = q;
  RotationParameterList[index].applied = true;
    
  // Compute the off-diagonal element in the new matrix
    rowp[p] = matrix[p][p] * CosTheta * CosTheta
    + matrix[q][q] * SinTheta * SinTheta
      + 2 * SinTheta * CosTheta * matrix[p][q];
    rowq[q] = matrix[q][q] * CosTheta * CosTheta
    + matrix[p][p] * SinTheta * SinTheta
      - 2 * SinTheta * CosTheta * matrix[p][q];
    rowp[q] = rowq[p] = 0;
    
    for ( int i = 0; i < size; i++ )
      {
        if ( ( i != p ) && ( i != q ) )
	  {
            rowp[i] = CosTheta * matrix[p][i] + SinTheta * matrix[q][i];
            rowq[i] = CosTheta * matrix[q][i] - SinTheta * matrix[p][i];
	  }
      }
    
    // Update the matrix
    for ( int i = 0; i < size; i++ )
      {
        matrix[p][i] = rowp[i];
        matrix[q][i] = rowq[i];
      }
    
    delete [] rowp;
    delete [] rowq;
}

int main(int argc,char** argv)
{
  int size, ListSize, NumberOfStep, NumberOfCombination, NumberOfPair;
  double epsilon;
  bool finished = false;
  double t_begin, t_end;
    
  cin >> size >> epsilon;
    
  double** matrix = new double*[size];
    
  for ( int i = 0; i < size; i++ )
    {
      matrix[i] = new double[size];
      for ( int j = 0; j < size; j++ )
        {
	  cin >> matrix[i][j];
        }
    }
    
  if ( size % 2 == 0 )
    {
      ListSize = size;
    }
  else
    {
      ListSize = size + 1;
    }
  NumberOfPair = ListSize / 2;
    
  NumberOfCombination = 1;
  for ( int i = ListSize - 1; i <= ListSize; i++ )
    {
      NumberOfCombination *= i;
    }
  NumberOfCombination /= 2;
  NumberOfStep = ListSize - 1;

  int** ListOfIndices = new int*[NumberOfStep];
    
  ListOfIndices[0] = new int[ListSize];
  for ( int j = 0; j < size; j++ )
    {
      ListOfIndices[0][j] = j;
    }
  if ( size % 2 == 1 )
    {
      ListOfIndices[0][ListSize - 1] = -1;
    }
    
  for ( int i = 1; i < NumberOfStep; i++ )
    {
      ListOfIndices[i] = new int[ListSize];
      ListOfIndices[i][0] = 0;
      for ( int j = 1; j < ListSize - 2; j = j + 2 )
        {
	  ListOfIndices[i][j] = ListOfIndices[i - 1][j + 2];
        }
      ListOfIndices[i][ListSize - 1] = ListOfIndices[i - 1][ListSize - 2];
        
      for ( int j = ListSize - 2; j > 2; j = j - 2 )
        {
	  ListOfIndices[i][j] = ListOfIndices[i - 1][j - 2];
        }
      ListOfIndices[i][2] = ListOfIndices[i - 1][1];
    }
    
  RotationParameter* RotationParameterList = new RotationParameter[NumberOfPair];
    
  omp_lock_t mylock;
  omp_init_lock(&mylock);

#pragma omp parallel default(shared)
  {
#pragma omp single
    {
      t_begin = omp_get_wtime();
    }
    while ( !finished )
      {
#pragma omp barrier
	finished = true;
	for ( int i = 0; i < NumberOfStep; i++ )
	  {
                
#pragma omp for schedule(guided)
	    for ( int j = 0; j < NumberOfPair; j++ )
	      {
		int p = ListOfIndices[i][j * 2];
		int q = ListOfIndices[i][j * 2 + 1];
		if ( ( p != -1 ) && ( q != -1 ) )
		  {
		    if ( abs( matrix[p][q] ) > epsilon )
		      {
			finished = false;
			// This rotate fxn only applies to row p and row q
			// but not to column p and column q
			rotate( p, q, matrix, size, epsilon, j, RotationParameterList );
		      }
		    else
		      {
			RotationParameterList[j].applied = false;
		      }
		  }
		else
		  {
		    RotationParameterList[j].applied = false;
		  }
	      }
                
#pragma omp for schedule(guided)
	    for ( int k = 0; k < NumberOfPair; k++ )
	      {
		if ( RotationParameterList[k].applied )
		  {
		    RotateMore( matrix, size, k, RotationParameterList );
		  }
	      }
                
	  }
      }
#pragma omp single
    {
      t_end = omp_get_wtime();
      printf("Elapsed wall time: %f\n", t_end - t_begin);
    }
  }
    
  // Print the matrix
  for ( int i = 0; i < size; i++ )
    {
      for ( int j = 0; j < size; j++ )
        {
	  cout << matrix[i][j] << "   ";
        }
      cout << endl;
    }

  // Free memory
  omp_destroy_lock(&mylock);

  for ( int i = 0; i < size; i++ )
    {
      delete [] matrix[i];
    }
  delete [] matrix;
    
  for ( int i = 0; i < NumberOfStep; i++ )
    {
      delete [] ListOfIndices[i];
    }
  delete [] ListOfIndices;
    
  delete [] RotationParameterList;
    
  return 0;
}
