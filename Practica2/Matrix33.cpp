#include "Matrix33.h"

const Matrix33 Matrix33::ZERO(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f);
const Matrix33 Matrix33::IDENTITY(1.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,1.0f);


Matrix33 Matrix33::operator*(const Matrix33& m) const
{
	return Matrix33(
      val[0][0] * m[0][0] + val[0][1] * m[1][0] + val[0][2] * m[2][0],
		val[0][0] * m[0][1] + val[0][1] * m[1][1] + val[0][2] * m[2][1],
		val[0][0] * m[0][2] + val[0][1] * m[1][2] + val[0][2] * m[2][2],
      val[1][0] * m[0][0] + val[1][1] * m[1][0] + val[1][2] * m[2][0],
      val[1][0] * m[0][1] + val[1][1] * m[1][1] + val[1][2] * m[2][1],
      val[1][0] * m[0][2] + val[1][1] * m[1][2] + val[1][2] * m[2][2],
      val[2][0] * m[0][0] + val[2][1] * m[1][0] + val[2][2] * m[2][0],
      val[2][0] * m[0][1] + val[2][1] * m[1][1] + val[2][2] * m[2][1],
      val[2][0] * m[0][2] + val[2][1] * m[1][2] + val[2][2] * m[2][2]);
}

Matrix33 Matrix33::transposeMult(const Matrix33& m) const
{
  return Matrix33(val[0][0] * m[0][0] + val[0][1] * m[0][1] + val[0][2] * m[0][2],
                  val[0][0] * m[1][0] + val[0][1] * m[1][1] + val[0][2] * m[1][2],
                  val[0][0] * m[2][0] + val[0][1] * m[2][1] + val[0][2] * m[2][2],
                  val[1][0] * m[0][0] + val[1][1] * m[0][1] + val[1][2] * m[0][2],
                  val[1][0] * m[1][0] + val[1][1] * m[1][1] + val[1][2] * m[1][2],
                  val[1][0] * m[2][0] + val[1][1] * m[2][1] + val[1][2] * m[2][2],
                  val[2][0] * m[0][0] + val[2][1] * m[0][1] + val[2][2] * m[0][2],
                  val[2][0] * m[1][0] + val[2][1] * m[1][1] + val[2][2] * m[1][2],
                  val[2][0] * m[2][0] + val[2][1] * m[2][1] + val[2][2] * m[2][2]);
}

Matrix33 Matrix33::getInverse(void) const
{
	float d = getDet();

   // assert(fabs(d) > 0.0);
  
  if (d == 0.0f) 
  {	
      // hack
		return ZERO;
	}

	d = 1.0f / d;

	Matrix33 t(val[1][1]*val[2][2] - val[1][2]*val[2][1],
             val[0][2]*val[2][1] - val[0][1]*val[2][2],
             val[0][1]*val[1][2] - val[0][2]*val[1][1],
             val[1][2]*val[2][0] - val[1][0]*val[2][2],
             val[0][0]*val[2][2] - val[0][2]*val[2][0],
             val[0][2]*val[1][0] - val[0][0]*val[1][2],
             val[1][0]*val[2][1] - val[1][1]*val[2][0],
             val[0][1]*val[2][0] - val[0][0]*val[2][1],
             val[0][0]*val[1][1] - val[0][1]*val[1][0]);
   
	return t * d;
}

Matrix33 Matrix33::interpolate(const Matrix33& m0, const Matrix33& m1, float f)
{
  if (f < 0.0f) f = 0.0f;
  else if (f > 1.0f) f = 1.0f;
  
  return m0*(1.0f - f) + m1*f;
}

void Matrix33::rotateX(float angle)
{
  float sinAngle = sin(angle);
  float cosAngle = cos(angle);
  rotate(val[1][0], val[2][0], sinAngle, cosAngle);
  rotate(val[1][1], val[2][1], sinAngle, cosAngle);
  rotate(val[1][2], val[2][2], sinAngle, cosAngle);
}

void Matrix33::rotateY(float angle)
{
  float sinAngle = sin(angle);
  float cosAngle = cos(angle);
  rotate(val[2][0], val[0][0], sinAngle, cosAngle);
  rotate(val[2][1], val[0][1], sinAngle, cosAngle);
  rotate(val[2][2], val[0][2], sinAngle, cosAngle);
}

void Matrix33::rotateZ(float angle)
{
  float sinAngle = sin(angle);
  float cosAngle = cos(angle);
  rotate(val[0][0], val[1][0], sinAngle, cosAngle);
  rotate(val[0][1], val[1][1], sinAngle, cosAngle);
  rotate(val[0][2], val[1][2], sinAngle, cosAngle);
}

