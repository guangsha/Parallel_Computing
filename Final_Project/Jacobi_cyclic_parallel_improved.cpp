////////////////////////////////////////////////////
// Improved cyclic Jacobi Algorithm for 
// diagonalization of symmetric matrices (parallel)                         
////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "omp.h"
#include <vector>

using namespace std;

struct RotationParameter
{
  bool applied;
  int p;
  int q;
  double SinTheta;
  double CosTheta;
};

void RotateMore( double **matrix, int size, int m, RotationParameter *rp )
{
  int p, q;
  double SinTheta, CosTheta;
  double newmp;
    
  p = rp->p;
  q = rp->q;
  SinTheta = rp->SinTheta;
  CosTheta = rp->CosTheta;
    
  newmp = CosTheta * matrix[m][p] + SinTheta * matrix[m][q];
  matrix[m][q] = CosTheta * matrix[m][q] - SinTheta * matrix[m][p];
  matrix[m][p] = newmp;
}

void rotate( int p, int q, double **matrix, int size, RotationParameter *rp )
{
  double theta, SinTheta, CosTheta;
  double *rowp = new double[size];
    
  // Compute the orthogonal matrix for Jacobi rotation
  theta = atan( 2 * matrix[p][q] / ( matrix[p][p] - matrix[q][q] ) ) / 2;
  rp->SinTheta = SinTheta = sin(theta);
  rp->CosTheta = CosTheta = cos(theta);
  rp->p = p;
  rp->q = q;
  rp->applied = true;
    
  // Compute the off-diagonal element in the new matrix
    rowp[p] = matrix[p][p] * CosTheta * CosTheta
    + matrix[q][q] * SinTheta * SinTheta
      + 2 * SinTheta * CosTheta * matrix[p][q];
    matrix[q][q] = matrix[q][q] * CosTheta * CosTheta
    + matrix[p][p] * SinTheta * SinTheta
      - 2 * SinTheta * CosTheta * matrix[p][q];
    rowp[q] = matrix[q][p] = 0;
    
    for ( int i = 0; i < size; i++ )
      {
        if ( ( i != p ) && ( i != q ) )
	  {
            rowp[i] = CosTheta * matrix[p][i] + SinTheta * matrix[q][i];
            matrix[q][i] = CosTheta * matrix[q][i] - SinTheta * matrix[p][i];
	  }
      }
    
    // Update the matrix
    for ( int i = 0; i < size; i++ )
      {
        matrix[p][i] = rowp[i];
      }
    
    delete [] rowp;
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
    
  vector <RotationParameter> RotationParameterList;
    
  omp_lock_t mylock;
  omp_init_lock(&mylock);
    
#pragma omp parallel default(shared)
  {
    int NumberOfThread = omp_get_num_threads();
    int tid = omp_get_thread_num();
        
    int StartPairId = ceil( tid * double(NumberOfPair) / NumberOfThread );
    int EndPairId = ceil( ( tid + 1 ) * double(NumberOfPair) / NumberOfThread ) - 1;
        
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
#pragma omp single
	    {
	      RotationParameterList.clear();
	    }
	    for ( int j = StartPairId; j <= EndPairId; j++ )
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
			RotationParameter rp;
			rotate( p, q, matrix, size, &rp );
			omp_set_lock(&mylock);
			RotationParameterList.push_back(rp);
			omp_unset_lock(&mylock);
		      }
		    else
		      {
			RotationParameter rp;
			rp.applied = false;
			omp_set_lock(&mylock);
			RotationParameterList.push_back(rp);
			omp_unset_lock(&mylock);
		      }
		  }
		else
		  {
		    RotationParameter rp;
		    rp.applied = false;
		    omp_set_lock(&mylock);
		    RotationParameterList.push_back(rp);
		    omp_unset_lock(&mylock);
		  }
	      }
	    int count = 0;
	    int parametersize = 0;
	    while ( count < NumberOfPair )
	      {
		parametersize = RotationParameterList.size();
		if ( count < parametersize )
		  {
		    if ( RotationParameterList[count].applied )
		      {
			for ( int j = StartPairId; j <= EndPairId; j++ )
			  {
			    int p = ListOfIndices[i][j * 2];
			    int q = ListOfIndices[i][j * 2 + 1];
			    if ( ( p == -1 ) )
			      {
				RotateMore( matrix, size, q, &RotationParameterList[count] );
			      }
			    else if ( ( q == -1 ) )
			      {
				RotateMore( matrix, size, p, &RotationParameterList[count] );
			      }
			    else
			      {
				if ( ( q != RotationParameterList[count].p )
				     && ( q != RotationParameterList[count].q ) )
				  {
				    RotateMore( matrix, size, p, &RotationParameterList[count] );
				    RotateMore( matrix, size, q, &RotationParameterList[count] );
				  }
			      }
			  }
		      }
		    count++;
		  }
	      }
#pragma omp barrier
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
    
  return 0;
}
