#pragma once

#include <memory.h>

#include "Vector3.h"

class Matrix33 
{
protected:
   float val[3][3];

public:
   static const Matrix33 ZERO;
   static const Matrix33 IDENTITY;

   Matrix33(void){}

   Matrix33(float v00, float v01, float v02, float v10, float v11, float v12, float v20, float v21, float v22)
   {
      set(v00, v01, v02, v10, v11, v12, v20, v21, v22);
   }

   Matrix33(const Matrix33& m)
   {
      (*this) = m;
   }

   Matrix33(const Vector3& v0, const Vector3& v1, const Vector3& v2)
   {
      setColumns(v0, v1, v2);
   }

   Matrix33(bool zero){
      if (zero)
         (*this) = ZERO;	   
      else{
         (*this) = IDENTITY;
      }			
   }

   ~Matrix33(void){}

   static inline Matrix33 ToCrossProd(const Vector3& v);

   static inline Matrix33 FromOuterProd(const Vector3& a, const Vector3& b);

   static Matrix33 interpolate(const Matrix33& m0, const Matrix33& m1, float f);

   //Get methods
   const float* operator[](int i) const
   {
      return val[i];
   }

   Vector3 getColumn(int i) const
   {
      return Vector3(val[0][i], val[1][i], val[2][i]);
   }

   Vector3 getRow(int i) const
   {
      return Vector3(val[i][0], val[i][1], val[i][2]);
   }

   //Set methods
   float* operator[](int i)
   {
      return val[i];
   }

   Matrix33& operator=(const Matrix33& m)
   {
      memcpy(val,m.val,9*sizeof(float));
      return *this;
   }

   void set(float v00, float  v01, float v02, float v10, float v11, float v12, float v20, float v21, float v22)
   {
      val[0][0] = v00; val[0][1] = v01; val[0][2] = v02;
      val[1][0] = v10; val[1][1] = v11; val[1][2] = v12;
      val[2][0] = v20; val[2][1] = v21; val[2][2] = v22;
   }

   void setColumns(const Vector3& v0, const Vector3& v1, const Vector3& v2)
   {
      set(v0[0], v1[0], v2[0], v0[1], v1[1], v2[1], v0[2], v1[2], v2[2]);
   }

   void setColumn(int i, const Vector3& v)
   {
      val[0][i] = v[0]; val[1][i] = v[1]; val[2][i] = v[2];
   }

   void setRow(int i, const Vector3& v)
   {
      val[i][0] = v[0]; val[i][1] = v[1]; val[i][2] = v[2];
   }

   //Arithmetic operators
   inline Matrix33 operator-(void) const;

   inline Matrix33 operator+(const Matrix33& m) const;

   inline Matrix33 operator-(const Matrix33& m) const;

   inline Matrix33 operator*(float f) const;

   Matrix33 operator*(const Matrix33& m) const;

   Matrix33 transposeMult(const Matrix33& m) const;

   inline Matrix33& operator+=(const Matrix33& m);

   inline Matrix33& operator-=(const Matrix33& m);

   inline Matrix33& operator*=(float f);

   inline Matrix33& operator*=(const Matrix33& m);

   inline Vector3 operator*(const Vector3& v) const;

   inline Vector3 transposeMult(const Vector3& v) const;

   inline void addDiagonal(float v);

   inline void subDiagonal(float v);

   //Friend arithmetic operators
   friend inline Matrix33 operator*(float f, const Matrix33& m);

   //Other operations
   Matrix33 getInverse(void) const;

   inline Matrix33 getTranspose(void) const;

   inline Matrix33 getAdjointTranspose(void) const;

   inline float getDet(void) const;

   inline float getFrobeniusNorm(void) const;

   inline float getNormInf(void) const;

   inline float getNormOne(void) const;

   static inline void rotate(float& a, float& b, float sinAngle, float cosAngle);
   void rotateX(float angle);
   void rotateY(float angle);
   void rotateZ(float angle);

};

inline std::ostream & operator<<(std::ostream & out, const Matrix33 & m)
{
   out << m[0][0] << " " << m[0][1] << " " << m[0][2] << " ";
   out << m[1][0] << " " << m[1][1] << " " << m[1][2] << " ";
   out << m[2][0] << " " << m[2][1] << " " << m[2][2] << " ";
   return out;
}

inline std::istream & operator>>(std::istream & in, Matrix33 & m)
{
   in >> m[0][0] >> m[0][1] >> m[0][2];
   in >> m[1][0] >> m[1][1] >> m[1][2];
   in >> m[2][0] >> m[2][1] >> m[2][2];
   return in;
}

inline Matrix33 Matrix33::ToCrossProd(const Vector3& v)
{
   return Matrix33(0.0f, -v[2], v[1], v[2], 0.0f, -v[0], -v[1], v[0], 0.0f);
}

inline Matrix33 Matrix33::FromOuterProd(const Vector3& a, const Vector3& b)
{
   return Matrix33(
      a[0] * b[0], a[0] * b[1], a[0] * b[2],
      a[1] * b[0], a[1] * b[1], a[1] * b[2],
      a[2] * b[0], a[2] * b[1], a[2] * b[2]);
}

inline Matrix33 Matrix33::operator-(void) const
{
   return Matrix33(
      -val[0][0], -val[0][1], -val[0][2],
      -val[1][0], -val[1][1], -val[1][2],
      -val[2][0], -val[2][1], -val[2][2]);
}

inline Matrix33 Matrix33::operator+(const Matrix33& m) const
{
   return Matrix33(
      val[0][0] + m[0][0], val[0][1] + m[0][1], val[0][2] + m[0][2],
      val[1][0] + m[1][0], val[1][1] + m[1][1], val[1][2] + m[1][2],
      val[2][0] + m[2][0], val[2][1] + m[2][1], val[2][2] + m[2][2]);
}

inline Matrix33 Matrix33::operator-(const Matrix33& m) const
{
   return Matrix33(
      val[0][0] - m[0][0], val[0][1] - m[0][1], val[0][2] - m[0][2],
      val[1][0] - m[1][0], val[1][1] - m[1][1], val[1][2] - m[1][2],
      val[2][0] - m[2][0], val[2][1] - m[2][1], val[2][2] - m[2][2]);
}

inline Matrix33 Matrix33::operator*(float f) const
{
   return Matrix33(
      val[0][0]*f, val[0][1]*f, val[0][2]*f,
      val[1][0]*f, val[1][1]*f, val[1][2]*f,
      val[2][0]*f, val[2][1]*f, val[2][2]*f);
}

inline Matrix33& Matrix33::operator+=(const Matrix33& m)
{
   val[0][0] += m[0][0];   val[0][1] += m[0][1];   val[0][2] += m[0][2];
   val[1][0] += m[1][0];   val[1][1] += m[1][1];   val[1][2] += m[1][2];
   val[2][0] += m[2][0];   val[2][1] += m[2][1];   val[2][2] += m[2][2];
   return *this;
}

inline Matrix33& Matrix33::operator-=(const Matrix33& m)
{
   val[0][0] -= m[0][0];   val[0][1] -= m[0][1];   val[0][2] -= m[0][2];
   val[1][0] -= m[1][0];   val[1][1] -= m[1][1];   val[1][2] -= m[1][2];
   val[2][0] -= m[2][0];   val[2][1] -= m[2][1];   val[2][2] -= m[2][2];
   return *this;
}

inline Matrix33& Matrix33::operator*=(float f)
{
   val[0][0] *= f;   val[0][1] *= f;   val[0][2] *= f;
   val[1][0] *= f;   val[1][1] *= f;   val[1][2] *= f;
   val[2][0] *= f;   val[2][1] *= f;   val[2][2] *= f;
   return *this;
}

inline Matrix33& Matrix33::operator*=(const Matrix33& m)
{
   *this = *this * m;
   return *this;
}

inline Vector3 Matrix33::operator*(const Vector3& v) const
{
   return Vector3(
      val[0][0]*v[0] + val[0][1]*v[1] + val[0][2]*v[2],
      val[1][0]*v[0] + val[1][1]*v[1] + val[1][2]*v[2],
      val[2][0]*v[0] + val[2][1]*v[1] + val[2][2]*v[2]);
}

inline Vector3 Matrix33::transposeMult(const Vector3& v) const
{
   return Vector3(
      val[0][0]*v[0] + val[1][0]*v[1] + val[2][0]*v[2],
      val[0][1]*v[0] + val[1][1]*v[1] + val[2][1]*v[2],
      val[0][2]*v[0] + val[1][2]*v[1] + val[2][2]*v[2]);
}

inline Matrix33 operator*(float f, const Matrix33& m)
{
   return Matrix33(m[0][0]*f, m[0][1]*f, m[0][2]*f,
      m[1][0]*f, m[1][1]*f, m[1][2]*f,
      m[2][0]*f, m[2][1]*f, m[2][2]*f);
}

inline void Matrix33::addDiagonal(float v) 
{
   val[0][0] += v;   val[1][1] += v;   val[2][2] += v;
}

inline void Matrix33::subDiagonal(float v) 
{
   val[0][0] -= v;   val[1][1] -= v;   val[2][2] -= v;
}

inline Matrix33 Matrix33::getTranspose(void) const
{
   return Matrix33(
      val[0][0], val[1][0], val[2][0],
      val[0][1], val[1][1], val[2][1],
      val[0][2], val[1][2], val[2][2]);
}

inline Matrix33 Matrix33::getAdjointTranspose(void) const
{
   return Matrix33(
      Vector3::crossProd(getColumn(1), getColumn(2)),
      Vector3::crossProd(getColumn(2), getColumn(0)),
      Vector3::crossProd(getColumn(0), getColumn(1)));
}

inline float Matrix33::getDet(void) const
{
   return val[0][0]*val[1][1]*val[2][2] + val[0][1]*val[1][2]*val[2][0] +
      val[0][2]*val[1][0]*val[2][1] - val[0][2]*val[1][1]*val[2][0] -
      val[0][1]*val[1][0]*val[2][2] - val[0][0]*val[1][2]*val[2][1];
}

inline float Matrix33::getFrobeniusNorm(void) const
{
   return sqrt(val[0][0]*val[0][0] + val[0][1]*val[0][1] + val[0][2]*val[0][2] +
      val[1][0]*val[1][0] + val[1][1]*val[1][1] + val[1][2]*val[1][2] +
      val[2][0]*val[2][0] + val[2][1]*val[2][1] + val[2][2]*val[2][2]);
}

inline float Matrix33::getNormInf(void) const
{
   float max = 0;
   for (int i=0;i<3;i++)
   {
      float sum = fabs(val[i][0])+fabs(val[i][1])+fabs(val[i][2]);
      if (max<sum) max = sum;
   }
   return max;
}

inline float Matrix33::getNormOne(void) const
{
   float max = 0;
   for (int i=0;i<3;i++)
   {
      float sum = fabs(val[0][i])+fabs(val[1][i])+fabs(val[2][i]);
      if (max<sum) max = sum;
   }
   return max;
}

inline void Matrix33::rotate(float& a, float& b, float sinAngle, float cosAngle)
{
   float t = a;
   a = cosAngle * a - sinAngle * b;
   b = sinAngle * t + cosAngle * b;
}

