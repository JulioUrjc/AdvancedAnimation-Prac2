#pragma once

class Vector3;
class Matrix33;

class Vector {
public:
   Vector( int size = 0   );
   Vector( const Vector & );
   Vector( const float *x, int n );
   Vector &operator=( const Vector & );
   void    operator=( float );
   void    SetSize( int );
   static  const Vector Null;

   void AddBlock3(int i, const Vector3& b);
   Vector3 GetBlock3(int i) const;

public: // Inlined functions.
   inline float  operator()( int i ) const { return elem[i]; }
   inline float& operator()( int i )       { return elem[i]; }
   inline float* Array() const { return elem; }
   inline int    Size () const { return size; }
   inline ~Vector() { delete[] elem; }

private:
   void   Create( int n = 0 ) { size = n; elem = new float[n]; }
   int    size;
   float* elem;
};

class MatrixMN {
public:
   MatrixMN( const MatrixMN & );
   MatrixMN( int num_rows = 0, int num_cols = 0, float value = 0.0 );
   ~MatrixMN();
   MatrixMN &operator=( const MatrixMN &M );
   MatrixMN &operator=( float s );
   void    SetSize( int rows, int cols = 0 );
   static  const MatrixMN Null;

   void AddBlock13(int i, int j, const Vector3& v);
   void AddBlock33(int i, int j, const Matrix33& m);
   MatrixMN operator*(const MatrixMN& M) const;
   MatrixMN getTranspose(void) const;

public: // Inlined functions.
   inline float  operator()( int i, int j ) const { return elem[ i * cols + j ]; }
   inline float &operator()( int i, int j )       { return elem[ i * cols + j ]; }
   inline int    Rows  () const { return rows; }
   inline int    Cols  () const { return cols; }
   inline float *Array () const { return elem; }

private:
   int    rows; // Number of rows in the matrix.
   int    cols; // Number of columns in the matrix.
   float *elem; // Pointer to the actual data.
};

int GaussElimination(const MatrixMN &A, 
                     const Vector &b, // This is the right-hand side.
                     Vector &x // This is the vector we are solving for.
                     );

void ProjectedGaussSeidel(const MatrixMN &A, 
                         const Vector &b, // This is the right-hand side.
                         Vector &x, // This is the vector we are solving for
                         bool *flags // Flags to indicate inequalities
                         );

