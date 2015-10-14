////////////////////////////////////////////////////
// Cyclic Jacobi Algorithm for diagonalization 
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
    rowp[q] = matrix[q][p] = 0;
    matrix[q][q] = matrix[q][q] * CosTheta * CosTheta
    + matrix[p][p] * SinTheta * SinTheta
      - 2 * SinTheta * CosTheta * matrix[p][q];
    
    for ( int i = 0; i < size; i++ )
      {
        if ( ( i != p ) && ( i != q ) )
	  {
            rowp[i] = CosTheta * matrix[p][i] + SinTheta * matrix[q][i];
            matrix[i][q] = matrix[q][i] = CosTheta * matrix[q][i] - SinTheta * matrix[p][i];
	  }
      }    
    // Update the matrix
    for ( int i = 0; i < size; i++ )
      {
        matrix[i][p] = matrix[p][i] = rowp[i];
      }
    
    delete [] rowp;
}

void jacobi( double **matrix, int size, double epsilon )
{
  double max, theta, SinTheta, CosTheta;
  int p, q;
    
  bool finished = false;
  int sum = 0;
  while ( !finished )
    {
      finished = true;
      for ( int p = 0; p < size; p++ )
	for ( int q = p + 1; q < size; q++ )
	  {
	    if ( abs( matrix[p][q] ) > epsilon )
	      {
		finished = false;
		rotate( p, q, matrix, size, epsilon );
	      }
	  }
    }
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
