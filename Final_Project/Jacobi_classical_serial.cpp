////////////////////////////////////////////////////
// Classical Jacobi Algorithm for diagonalization 
// of symmetric matrices (serial)                         
////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

using namespace std;

void rotate( int p, int q, double **matrix, int size, double epsilon )
{
  double theta, SinTheta, CosTheta;
  double *rowp = new double[size];
    
  // Compute the orthogonal matrix for Jacobi rotation
  theta = atan( 2 * matrix[p][q] / ( matrix[p][p] - matrix[q][q] ) ) / 2;
  SinTheta = sin(theta);
  CosTheta = cos(theta);
    
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
        matrix[i][p] = matrix[p][i] = rowp[i];
        matrix[i][q] = matrix[q][i];
      }
    
    delete [] rowp;
}

void jacobi( double **matrix, int size, double epsilon )
{
  double** NewMatrix = new double*[size];
  double max, theta, SinTheta, CosTheta;
  int p, q;
  // First find the largest off-diagonal element
  max = matrix[0][1]; // row = 0, column = 1
  p = 0;
  q = 1;
  for ( int i = 0; i < size; i++ )
    {
      for ( int j = i + 1; j < size; j++ )
        {
	  if ( abs( matrix[i][j] ) > max )
            {
	      max = abs( matrix[i][j] );
	      p = i;
	      q = j;
            }
        }
    }
    
  while ( max > epsilon )
    {
      rotate( p, q, matrix, size, epsilon );
        
      // Find the largest off-diagonal element
      max = matrix[0][1]; // row = 0, column = 1
      p = 0;
      q = 1;
      for ( int i = 0; i < size; i++ )
        {
	  for ( int j = i + 1; j < size; j++ )
            {
	      if ( abs( matrix[i][j] ) > max )
                {
		  max = abs( matrix[i][j] );
		  p = i;
		  q = j;
                }
            }
        }
    }
  //  cout << sum << endl;
}


int main(int argc,char** argv)
{
  int size;
  double epsilon;
    
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
    
  jacobi( matrix, size, epsilon );
    
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
  for ( int i = 0; i < size; i++ )
    {
      delete [] matrix[i];
    }
  delete [] matrix;
    
  return 0;
}
