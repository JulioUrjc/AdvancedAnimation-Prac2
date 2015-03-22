#pragma once

#include <vector>
using namespace std;

#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/projectors/SbPlaneProjector.h>

#include "Vector3.h"
#include "Matrix33.h"

#include "RigidBody.h"

struct PickObject
{
   //Ideally, this should have a type and an id
   static enum objType {SPHERE};
   objType type;
   int id;

   PickObject(objType type, int id) : type(type), id(id) {}
};

class Scene
{
public:
   //General data
   bool constraint;
   bool collisions;

protected:
   //General data
   float step;
   int nbodies;
   float floor;

   //Objects and properties
   vector<RigidBody> rbodies;

   //Grabbing
   bool grabbed;
   SbPlaneProjector *grabProjector;
   Vector3 grabPos;

   //Rendering members
   vector<SoSphere*> point1r, point2r;
   vector<SoCylinder*> bodyar;
   vector<SoCube*> bodybr;
   
   //Transformations
   vector<SoTranslation*> point1XForm, point2XForm;
   vector<SoTransform*> bodyXForm;

   //Animation
   bool pause;

public:
   Scene(void);
   Scene(int argc, char* argv[]);
   ~Scene(void);

   //Initialization
   void Init(void);
   void PrintSettings(void);
   void Add(SoSeparator *root);

   //Update scene
   void Update(void);
   void Animate(void);
   void Pause(void);

   //Grabbing
   bool Grabbed(void) const
   {
      return grabbed;
   }
   void Grab(void *pickedObj, SbVec2f& screenPos, SbPlaneProjector *proj);
   void UnGrab(void);
   void UpdateGrab(SbVec2f& screenPos);

};
