////////////////////////////////////////////////////
// Generate N*N matrix with random numbers from 
// 0 to 10 for diagonal elements and random numbers 
// from 0 to 1 for off-diagonal elements
////////////////////////////////////////////////////

#include <iostream>
#include <random>

using namespace std;

int main(int argc,char** argv)
{
  int size = atoi(argv[1]);
  int seed = atoi(argv[2]);
    
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> dis(0, 10);
    
  double** matrix = new double*[size];
  for ( int i = 0; i < size; i++ )
    {
      matrix[i] = new double[size];
    }
    
  for ( int i = 0; i < size; i++ )
    {
      matrix[i][i] = dis(gen);
      for ( int j = i + 1; j < size; j++ )
        {
	  matrix[j][i] = matrix[i][j] = dis(gen) / 10;
        }
    }
    
  printf("%d\t%20.15f\n", size, 0.0000000001);
  for ( int i = 0; i < size; i++ )
    {
      for ( int j = 0; j < size; j++ )
        {
	  printf("%4.2f\t", matrix[i][j]);
        }
      printf("\n");
    }
    
  for ( int i = 0; i < size; i++ )
    {
      delete [] matrix[i];
    }
  delete [] matrix;
    
  return 0;
}
