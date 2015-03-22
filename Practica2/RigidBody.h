#pragma once;

#include "Vector3.h"
#include "Matrix33.h"

class RigidBody
{
protected:
   Vector3 pos, vel, omega;
   Matrix33 rot;
   Matrix33 inertia0, inertia, inertia0inv, inertiainv;
   float mass;

public:
   static const RigidBody RIGIDBODYVOID;

   RigidBody(void){}
   ~RigidBody(void){}
   void Init(const Vector3& pos, const Matrix33& rot, const Vector3& vel, const Vector3& omega, const Matrix33& inertia0, float mass)
   {
      this->pos = pos;
      this->rot = rot;
      this->vel = vel;
      this->omega = omega;
      this->mass = mass;
      this->inertia0 = inertia0;
      inertia0inv = inertia0.getInverse();
      Matrix33 rotT = rot.getTranspose();
      inertia = rot * inertia0 * rotT;
      inertiainv = rot * inertia0inv * rotT;
   }

   void UpdateInertia(void)
   {
      Matrix33 rotT = rot.getTranspose();
      inertia = rot * inertia0 * rotT;
      inertiainv = rot * inertia0inv * rotT;
   }

   const Vector3& Position(void) const { return pos; }
   const Matrix33& Rotation(void) const { return rot; }
   const Vector3& Velocity(void) const { return vel; }
   const Vector3& Omega(void) const { return omega; }
   float Mass(void) const { return mass;}
   const Matrix33& Inertia(void) const { return inertia; }
   const Matrix33& InertiaInv(void) const { return inertiainv; }

   void SetPosition(const Vector3& pos)
   {
      this->pos = pos;
   }
   void SetRotation(const Matrix33& rot)
   {
      this->rot = rot;
   }
   void SetVelocity(const Vector3& vel)
   {
      this->vel = vel;
   }
   void SetOmega(const Vector3& omega)
   {
      this->omega = omega;
   }

   Vector3 posLocalToGlobal(const Vector3& p) const //p is a point expressed in local coordinates of the body
   {
      return pos + rot*p;
   }
   Vector3 vecLocalToGlobal(const Vector3& v) const //v is a vector expressed in local coordinates of the body
   {
      return rot*v;
   }
   Vector3 pointVelocity(const Vector3& p) const //p is a point expressed in local coordinates of the body
   {
      return vel + Vector3::crossProd(omega, rot*p);
   }

};
