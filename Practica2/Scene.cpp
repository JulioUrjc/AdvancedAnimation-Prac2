#include ".\Scene.h"
#include ".\Exercise.h"

#include <iostream>
using namespace std;

#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/SbTime.h>

Scene::Scene(void)
{
   Init();
   PrintSettings();
}

Scene::Scene(int argc, char* argv[])
{
	// defaults:
   step = 0.01f;
   nbodies = 1;
   floor = -0.5f * nbodies;

   int arg=1;
   while(arg<argc)
   {
      //Step size
      if(!strcmp(argv[arg], "-step"))
      {
         step=(float)atof(argv[++arg]);
         arg++;
      }
      else if(!strcmp(argv[arg], "-bodies"))
      {
         nbodies=(int)atoi(argv[++arg]);
         arg++;
      }
      //Others
      else
      {
         //Print instructions
         cerr << endl << "Unrecognized option " << argv[arg] << endl;
         cerr << "Usage: Test.exe -[option1] [settings] -[option2] [settings] ..." << endl;
         cerr << "Options:" << endl;
         cerr << "\t-step [step size in secs]" << endl;
         cerr << "\t-bodes [number of bodies]" << endl;
         break;
      }
   }

   PrintSettings();
   Init();
}

Scene::~Scene(void)
{
   rbodies.clear();
}

void Scene::PrintSettings(void)
{
   cerr << endl << "Current Settings:" << endl;
   cerr << "\t-step " << step << endl;
}

void Scene::Init(void)
{
   //Animation settings
   constraint=true;
   collisions=true;
   pause=true;
   grabbed=false;
   grabProjector=NULL;

   //Initialize simulation state
   rbodies.reserve(nbodies); rbodies.resize(nbodies);
   float mass = 0.1f;
   float length = 1.0f;
   float radius = length/8.0f;
   Matrix33 inertia0 = mass * Matrix33(1.0f/2.0f*radius*radius, 0, 0, 0, 1.0f/12.0f*(3.0f*radius*radius+length*length), 0, 0, 0, 1.0f/12.0f*(3.0f*radius*radius+length*length));
   Matrix33 rot = Matrix33::IDENTITY;
   Vector3 vel = Vector3::ZERO;
   Vector3 omega = Vector3::ZERO;
   for (int i=0; i<nbodies; i++)
   {
      rbodies[i].Init((i+0.5f)*length*Vector3::UNIT_X, rot, vel, omega, inertia0, mass);
   }
   Vector3 axis = 0.5f*length*Vector3::UNIT_X;

   //Create spheres and transformations
   point1r.reserve(nbodies); point1r.resize(nbodies);
   point2r.reserve(nbodies); point2r.resize(nbodies);
   point1XForm.reserve(nbodies); point1XForm.resize(nbodies);
   point2XForm.reserve(nbodies); point2XForm.resize(nbodies);
   bodyar.reserve(nbodies); bodyar.resize(nbodies);
   bodybr.reserve(nbodies); bodybr.resize(nbodies);
   bodyXForm.reserve(nbodies); bodyXForm.resize(nbodies);
   for (int i=0; i<nbodies; i++)
   {
      point1r[i] = new SoSphere();
      point1r[i]->radius.setValue(0.03f);
      point1XForm[i] = new SoTranslation();
      point1XForm[i]->translation.setValue(rbodies[i].posLocalToGlobal(-axis).ToSbVec3f());
      point2r[i] = new SoSphere();
      point2r[i]->radius.setValue(0.03f);
      point2XForm[i] = new SoTranslation();
      point2XForm[i]->translation.setValue(rbodies[i].posLocalToGlobal(axis).ToSbVec3f());

      //Create cylinders, cubes and transformations
      bodyar[i] = new SoCylinder();
      bodyar[i]->radius.setValue(0.01f);
      bodyar[i]->height.setValue(2*axis.getLength());
      bodybr[i] = new SoCube();
      bodybr[i]->width.setValue(0.1f);
      bodybr[i]->depth.setValue(0.1f);
      bodybr[i]->height.setValue(1.2f*axis.getLength());
      bodyXForm[i] = new SoTransform();
      bodyXForm[i]->pointAt(rbodies[i].Position().ToSbVec3f(), rbodies[i].posLocalToGlobal(-axis).ToSbVec3f());
   }
}


void Scene::Update(void)
{
   if(pause)
   {
      return;
   }

   //Perform animation
   AdvanceTimeStep(rbodies, step, constraint, collisions, floor);

   Vector3 axis(-0.5f, 0, 0);
   for (int i=0; i<nbodies; i++)
   {
      point1XForm[i]->translation.setValue(rbodies[i].posLocalToGlobal(axis).ToSbVec3f());
      point2XForm[i]->translation.setValue(rbodies[i].posLocalToGlobal(-axis).ToSbVec3f());

      bodyXForm[i]->pointAt(rbodies[i].Position().ToSbVec3f(), rbodies[i].posLocalToGlobal(axis).ToSbVec3f());
      bodyar[i]->height=2.0f*rbodies[i].vecLocalToGlobal(axis).getLength();
   }

   // reduce frame rate a bit
   //SbTime t = SbTime::getTimeOfDay();
   //while(SbTime::getTimeOfDay().getValue()-t.getValue()<0.00125) {
	  // // delay
   //}
}


void Scene::Add(SoSeparator *root)
{
   //Add separator for bottom points
   SoSeparator *pointSeparator=new SoSeparator;
   root->addChild(pointSeparator);

   //Bottom point properties
   SoMaterial *pointMaterial=new SoMaterial;
   pointMaterial->diffuseColor.setValue(0.9f, 0.0f, 0.0f);
   pointSeparator->addChild(pointMaterial);

   //Add points with xforms
   for (int i=0; i<nbodies; i++)
   {
      SoSeparator *pSeparator=new SoSeparator();
      pointSeparator->addChild(pSeparator);
      pSeparator->addChild(point2XForm[i]);
      pSeparator->addChild(point2r[i]);
   }

   //Add separator for top points
   pointSeparator=new SoSeparator;
   root->addChild(pointSeparator);

   //Top point properties
   pointMaterial=new SoMaterial;
   pointMaterial->diffuseColor.setValue(0.0f, 0.0f, 0.9f);
   pointSeparator->addChild(pointMaterial);

   //Add points with xforms
   for (int i=0; i<nbodies; i++)
   {
      SoSeparator *pSeparator=new SoSeparator;
      pointSeparator->addChild(pSeparator);
      pSeparator->addChild(point1XForm[i]);
      pSeparator->addChild(point1r[i]);
   }

   //Add separator for rigid bodies
   SoSeparator *bodySeparator=new SoSeparator;
   root->addChild(bodySeparator);

   //Spring properties
   SoMaterial *bodyMaterial=new SoMaterial;
   bodyMaterial->diffuseColor.setValue(0.9f, 0.9f, 0.9f);
   bodySeparator->addChild(bodyMaterial);

   //Create common rotation
   SoTransform *cylOrient=new SoTransform;
   cylOrient->pointAt(SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(0.0f, -1.0f, 0.0f));

   //Add bodys with xforms
   for (int i=0; i<nbodies; i++)
   {
      SoSeparator *sSeparator=new SoSeparator;
      bodySeparator->addChild(sSeparator);
      sSeparator->addChild(bodyXForm[i]);
      sSeparator->addChild(cylOrient);
      sSeparator->addChild(bodyar[i]);
      sSeparator->addChild(bodybr[i]);
   }
}

void Scene::Pause(void)
{
   pause=!pause;
}

void Scene::Grab(void *pickedObj, SbVec2f& screenPos, SbPlaneProjector *proj)
{
   if (((PickObject*)pickedObj)->type != PickObject::SPHERE && ((PickObject*)pickedObj)->id != 2) return;

   grabbed=true;
   
   grabProjector=proj;
   grabPos = Vector3(grabProjector->project(screenPos));
}

void Scene::UnGrab(void)
{
   grabbed=false;

   delete grabProjector;
   grabProjector=NULL;
}

void Scene::UpdateGrab(SbVec2f& screenPos)
{
   grabPos = Vector3(grabProjector->project(screenPos));
}


