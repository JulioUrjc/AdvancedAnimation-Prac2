#include "LinSys.h"
#include "Vector3.h"
#include "Matrix33.h"
#include <stdlib.h>

Vector::Vector( const float *x, int n )
{
   Create( n );
   for( register int i = 0; i < size; i++ ) elem[i] = x[i];
}

Vector::Vector( const Vector &A )
{
   Create( A.Size() );
   for( register int i = 0; i < A.Size(); i++ ) elem[i] = A(i);
}

Vector::Vector( int n )
{
   Create( n );
   for( register int i = 0; i < n; i++ ) elem[i] = 0.0;
}

void Vector::SetSize( int new_size )
{
   if( size != new_size )
   {
      delete[] elem;
      Create( new_size );
      for( register int i = 0; i < new_size; i++ ) elem[i] = 0.0;
   }
}

void Vector::AddBlock3(int i, const Vector3& b)
{
   for(int row=0; row<3; row++, i++)
   {
      elem[i] += b[row];
   }
}

Vector3 Vector::GetBlock3(int i) const
{
   return Vector3(elem[i], elem[i+1], elem[i+2]);
}

// Create a new matrix of the given size.  If n_cols is zero (the default), 
// it is assumed that the matrix is to be square; that is, n_rows x n_rows.  
// The matrix is filled with "value", which defaults to zero.
MatrixMN::MatrixMN( int n_rows, int n_cols, float value ) 
{
   rows = 0;
   cols = 0;
   elem = NULL;
   SetSize( n_rows, n_cols );
   float *e = elem;
   for( register int i = 0; i < rows * cols; i++ ) *e++ = value;
}

// Copy constructor.
MatrixMN::MatrixMN( const MatrixMN &M ) 
{
   rows = 0;
   cols = 0;
   elem = NULL;
   SetSize( M.Rows(), M.Cols() );
   register float *e = elem;
   register float *m = M.Array();
   for( register int i = 0; i < rows * cols; i++ ) *e++ = *m++;
}

MatrixMN::~MatrixMN() 
{
   SetSize( 0, 0 );
}

MatrixMN& MatrixMN::operator=(const MatrixMN& M)
{
   assert(Rows() == M.Rows() && Cols() == M.Cols());

   for (int i=0; i<Rows(); i++)
   {
      for (int j=0; j<Cols(); j++)
      {
         (*this)(i,j) = M(i,j);
      }
   }

   return *this;
}

// Re-shape the matrix.  If the number of elements in the new matrix is
// different from the original matrix, the original data is deleted and
// replaced with a new array.  If new_cols is zero (the default), it is
// assumed to be the same as new_rows -- i.e. a square matrix.
void MatrixMN::SetSize( int new_rows, int new_cols )
{
   if( new_cols == 0 ) new_cols = new_rows;
   int n = new_rows * new_cols;
   if( rows * cols != n )
   {
      if( elem != NULL ) delete[] elem;
      elem = ( n == 0 ) ? NULL : new float[ n ];
   }
   rows = new_rows;
   cols = new_cols;
}

void MatrixMN::AddBlock13(int i, int j, const Vector3& v)
{
   for(int col=0; col<3; col++, j++)
   {
      (*this)(i,j) += v[col];
   }
}

void MatrixMN::AddBlock33(int i, int j, const Matrix33& m)
{
   for(int row=0, ri=i; row<3; row++, ri++)
   {
      for(int col=0, cj=j; col<3; col++, cj++)
      {
         (*this)(ri,cj) += m[row][col];
      }
   }
}

MatrixMN MatrixMN::operator*(const MatrixMN& M) const
{
   assert(Cols() == M.Rows());

   MatrixMN R(Rows(), M.Cols());

   for (int i=0; i<Rows(); i++)
   {
      for (int j=0; j<M.Cols(); j++)
      {
         float val = 0;
         for (int k=0; k<Cols(); k++)
         {
            val += (*this)(i,k) * M(k, j);
         }
         R(i,j) = val;
      }
   }

   return R;
}

MatrixMN MatrixMN::getTranspose(void) const
{
   MatrixMN R(Cols(), Rows());

   for (int i=0; i<Cols(); i++)
   {
      for (int j=0; j<Rows(); j++)
      {
         R(i,j) = (*this)(j,i);
      }
   }

   return R;
}

int GaussElimination( const MatrixMN &A, const Vector &b, Vector &x)
{
   MatrixMN B( A );
   Vector c( b );
   x.SetSize( A.Cols() );
   int m = B.Rows();
   register int i, j, k;

   // Perform Gaussian elimination on the copies, B and c.

   for( i = 0; i < m; i++ )
   {
      for( j = i + 1; j < m; j++ )
      {
         float scale = -B(j,i) / B(i,i);
         for( k = i; k < m; k++ )
            B(j,k) += scale * B(i,k);
         B(j,i) = 0.0;
         c(j) += scale * c(i);
      }
   }

   // Now solve by back substitution.

   for( i = m - 1; i >= 0; i-- )
   {
      float a = 0.0;
      for( j = i + 1; j < m; j++ ) a += B(i,j) * x(j);
      x(i) = ( c(i) - a ) / B(i,i);
   }

   return 1;
}

void ProjectedGaussSeidel( const MatrixMN &A, const Vector &b, Vector &x, bool *flags)
{
   assert(A.Rows() == A.Cols());

   int niters = 100;
   for (int k=0; k<niters; k++)
   {
      for (int i=0; i<A.Rows(); i++)
      {
         float val = b(i);
         for (int j=0; j<i; j++)
         {
            val -= A(i,j)*x(j);
         }
         for (int j=i+1; j<A.Cols(); j++)
         {
            val -= A(i,j)*x(j);
         }
         x(i) = val / A(i,i);
         if (flags[i])
         {
            if (x(i) < 0.0f)
            {
               x(i) = 0.0f;
            }
         }
      }
   }
}

