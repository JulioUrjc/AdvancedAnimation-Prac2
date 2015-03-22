#ifndef _VECTOR_3_
#define _VECTOR_3_

#include "math.h"

#include <iostream>
using namespace std;

#include <Inventor/fields/SoSFVec3f.h>

class Vector3 
{
private:
   float val[3];

public:
   static const Vector3 ZERO;
   static const Vector3 UNIT_X;
   static const Vector3 UNIT_Y;
   static const Vector3 UNIT_Z;

   Vector3(void){}

   Vector3(float x, float y, float z)
   {
      set(x, y, z);
   }

   Vector3(const Vector3& v)
   {
      *this = v;
   }

   Vector3(const SbVec3f& v)
   {
      set(v[0], v[1], v[2]);
   }

   ~Vector3(void){}

   //Get methods
   const float * getReference(void) const
   {
      return val;
   }

   float* getReference()
   {
      return val;
   }

   float operator[](int i) const
   {
      return val[i];
   }

   bool isZero() const;

   SbVec3f ToSbVec3f(void) const
   {
      return SbVec3f(val[0], val[1], val[2]);
   }

   //Set methods
   float& operator[](int i)
   {
      return val[i];
   }

   Vector3& operator=(const Vector3& v)
   {
      set(v[0], v[1], v[2]);
      return *this;
   }

   void set(float x, float y, float z)
   {
      val[0]=x; val[1]=y; val[2]=z;
   }

   //Arithmetic operators
   inline Vector3 operator-(void) const;
   inline Vector3 operator+(const Vector3& v) const;
   inline Vector3 operator-(const Vector3& v) const;
   inline Vector3 operator*(const Vector3& v) const;
   inline Vector3 operator*(float f) const;

   friend Vector3 operator*(float f, const Vector3& v1);

   inline Vector3& operator+=(const Vector3& v);
   inline Vector3& operator-=(const Vector3& v);
   inline Vector3& operator*=(const Vector3& v);
   inline Vector3& operator*=(float f);

   inline float getSquaredLength() const;
   inline float getLength() const;

   void normalize();
   void rotate(float& x, float& y, float angle);
   void rotateX(float angle);
   void rotateY(float angle);
   void rotateZ(float angle);
   inline int maxDimension();

   static float angle(const Vector3& v0, const Vector3& v1);
   static Vector3 minVector(const Vector3& v0, const Vector3& v1);
   static Vector3 maxVector(const Vector3& v0, const Vector3& v1);
   static inline float distance(const Vector3& v0, const Vector3& v1);
   static inline float squaredDistance(const Vector3& v0, const Vector3& v1);
   static inline float dotProd(const Vector3& v0, const Vector3& v1);
   static inline Vector3 crossProd(const Vector3& v0, const Vector3& v1);
   static Vector3 interpolate(Vector3 v0, Vector3 v1, float f);

};

inline std::ostream & operator<<(std::ostream & out, const Vector3 & v)
{
   out << v[0] << " " << v[1] << " " << v[2];
   return out;
}

inline std::istream & operator>>(std::istream & in, Vector3 & v)
{
   in >> v[0] >> v[1] >> v[2];
   return in;
}

inline Vector3 Vector3::operator-(void) const
{
   return Vector3(-val[0], -val[1], -val[2]);
}

inline Vector3 Vector3::operator+(const Vector3& v) const
{
   return Vector3(val[0]+v[0], val[1]+v[1], val[2]+v[2]);
}

inline Vector3 Vector3::operator-(const Vector3& v) const
{
   return Vector3(val[0]-v[0], val[1]-v[1], val[2]-v[2]);
}

inline Vector3 Vector3::operator*(const Vector3& v) const
{
   return Vector3(val[0]*v[0], val[1]*v[1], val[2]*v[2]);
}

inline Vector3 Vector3::operator*(float f) const
{
   return Vector3(val[0]*f, val[1]*f, val[2]*f);
}

inline Vector3 operator*(float f, const Vector3& v)
{
   return Vector3(v[0]*f, v[1]*f, v[2]*f);
}

inline Vector3& Vector3::operator+=(const Vector3& v)
{
   val[0] += v[0];
   val[1] += v[1];
   val[2] += v[2];
   return *this;
}

inline Vector3& Vector3::operator-=(const Vector3& v)
{
   val[0] -= v[0];
   val[1] -= v[1];
   val[2] -= v[2];
   return *this;
}

inline Vector3& Vector3::operator*=(const Vector3& v)
{
   val[0] *= v[0];
   val[1] *= v[1];
   val[2] *= v[2];
   return *this;
}

inline Vector3& Vector3::operator*=(float f)
{
   val[0] *= f;
   val[1] *= f;
   val[2] *= f;
   return *this;
}

inline bool Vector3::isZero() const
{
   return (val[0] == 0.0f && val[1] == 0.0f && val[2] == 0.0f);
}

inline float Vector3::getLength() const
{
   return sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]);
}

inline float Vector3::getSquaredLength() const
{
   return val[0]*val[0] + val[1]*val[1] + val[2]*val[2];
}

inline void Vector3::rotate(float& x, float& y, float angle)
{
  float sinAngle = sin(angle);
	float cosAngle = cos(angle);
  float t = x;
  x = cosAngle*x - sinAngle*y;
  y = sinAngle*t + cosAngle*y;
}

inline void Vector3::rotateX(float angle)
{
   rotate(val[1], val[2], angle);
}

inline void Vector3::rotateY(float angle)
{
   rotate(val[2], val[0], angle);
}

inline void Vector3::rotateZ(float angle)
{
   rotate(val[0], val[1], angle);
}

inline float Vector3::distance(const Vector3& v0, const Vector3& v1)
{
   return (v0 - v1).getLength();
}

inline float Vector3::squaredDistance(const Vector3& v0, const Vector3& v1)
{
   return (v0 - v1).getSquaredLength();
}

inline float Vector3::dotProd(const Vector3& v0, const Vector3& v1)
{ 
   return v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]; 
}

inline Vector3 Vector3::crossProd(const Vector3& v0, const Vector3& v1)
{
   return Vector3(v0[1]*v1[2] - v1[1]*v0[2],
      v0[2]*v1[0] - v1[2]*v0[0],
      v0[0]*v1[1] - v1[0]*v0[1]);
}

inline void Vector3::normalize()
{
   float n = getLength();

   if (n == 0.0f) set(0.0f, -1.0f, 0.0f);
   else *this *= 1.0f / n;
}

inline Vector3 Vector3::minVector(const Vector3& v0, const Vector3& v1)
{
   Vector3 r(v0);

   if (v1[0] < v0[0]) r[0] = v1[0];
   if (v1[1] < v0[1]) r[1] = v1[1];
   if (v1[2] < v0[2]) r[2] = v1[2];

   return r;
}

inline Vector3 Vector3::maxVector(const Vector3& v0, const Vector3& v1)
{
   Vector3 r(v0);

   if (v1[0] > v0[0]) r[0] = v1[0];
   if (v1[1] > v0[1]) r[1] = v1[1];
   if (v1[2] > v0[2]) r[2] = v1[2];

   return r;
}

inline int Vector3::maxDimension()
{
   int idx=0;
   float max=val[0];
   for (int i=0; i<3; i++) {
      if (val[i] > max) {
         max = val[i];
         idx = i;
      }
   }

   return idx;
}

inline Vector3 Vector3::interpolate(Vector3 v0, Vector3 v1, float f)
{
   if (f < 0.0f) f = 0.0f;
   if (f > 1.0f) f = 1.0f;
   v0 *= 1.0f - f;
   v1 *= f;
   v0 += v1;
   return v0;
}

inline float Vector3::angle(const Vector3& v0, const Vector3& v1)
{
   float d = dotProd(v0, v1);
   if (d > 1.0f) return acos(1.0f);
   if (d < -1.0f) return acos(-1.0f);
   return acos(d);
}


#endif // _VECTOR_3_ 

