// coin3D includes
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>

#include <Inventor/sensors/SoIdleSensor.h> 
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/events/SoLocation2Event.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/projectors/SbPlaneProjector.h>
#include <Inventor/nodes/SoCamera.h>

#include ".\Scene.h"
#include <iostream>


struct StrAllData
{
   SoWinExaminerViewer *viewer;
   Scene *scene;
};
typedef struct StrAllData AllData ;


//Idle function
static void idlefunc(void *data, SoSensor* s)
{
   Scene *scene=((AllData*)data)->scene;

   scene->Update();

   s->schedule();
}

//Keyboard callback
static void keyboardfunc(void *data, SoEventCallback *eventCB)
{
   Scene *scene=((AllData*)data)->scene;
   const SoEvent *event = eventCB->getEvent();
   if(SO_KEY_PRESS_EVENT(event, S)){
      scene->Pause();
      eventCB->setHandled();
   }
   else if (SO_KEY_PRESS_EVENT(event, C))
   {
      scene->constraint = !scene->constraint;
   }
   else if (SO_KEY_PRESS_EVENT(event, F))
   {
      scene->collisions = !scene->collisions;
   }
}

//Mouse motion callback
static void motionfunc(void *data, SoEventCallback *eventCB)
{
   AllData *appData=(AllData*)data;
   SoWinExaminerViewer *viewer=appData->viewer;
   Scene *scene=appData->scene;

   if(!scene->Grabbed())
   {
      return;
   }

   const SoMouseButtonEvent *mbe=(SoMouseButtonEvent*)eventCB->getEvent();
   SbVec2f pos=mbe->getNormalizedPosition(viewer->getViewportRegion());
   scene->UpdateGrab(pos);

}

//Mouse callback
static void mousefunc(void *data, SoEventCallback *eventCB)
{
   AllData *appData=(AllData*)data;
   SoWinExaminerViewer *viewer=appData->viewer;
   Scene *scene=appData->scene;
   const SoMouseButtonEvent *mbe=(SoMouseButtonEvent*)eventCB->getEvent();

   //Handle point grabbing
   if(mbe->getButton() == SoMouseButtonEvent::BUTTON1 && mbe->getState() == SoButtonEvent::DOWN)
   {
      //Get viewport point and search for 3D picked point
      SoRayPickAction rp(viewer->getViewportRegion());
      rp.setPoint(mbe->getPosition());
      rp.apply(viewer->getSceneManager()->getSceneGraph());
      SoPickedPoint *point = rp.getPickedPoint();
      if(point == NULL)
      {
         return;
      }
      eventCB->setHandled();

      //Get the first picked node
      SoNode *node=point->getPath()->getTail();
      if(!node->getUserData()) return;

      //Set the projection plane
      SbPlane plane(-viewer->getCamera()->getViewVolume().getProjectionDirection(), point->getPoint());
      SbPlaneProjector *proj=new SbPlaneProjector(plane);
      proj->setViewVolume(viewer->getCamera()->getViewVolume());

      //Grab the corresponding point
      SbVec2f pos=mbe->getNormalizedPosition(viewer->getViewportRegion());
      scene->Grab(node->getUserData(), pos, proj);
   }

   //Handle ungrabbing
   if(mbe->getButton() == SoMouseButtonEvent::BUTTON1 && mbe->getState() == SoButtonEvent::UP)
   {
      if(scene->Grabbed())
      {
         scene->UnGrab();
      }
      eventCB->setHandled();
   }

}

int main(int argc, char* argv[])
{
   HWND window = SoWin::init(argv[0]);
   if (window==NULL) exit(1);

   SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);

   //Create root
   SoSeparator * root = new SoSeparator;
   root->ref();

   //Initialize scene and add to root
   Scene *scene=new Scene(argc, argv);
   scene->Add(root);

   //Create application data structure
   AllData appData;
   appData.scene=scene;
   appData.viewer=viewer;

   // Create idle sensor and pass pointer to scene
   SoIdleSensor *idle = new SoIdleSensor(idlefunc, &appData);
   idle->schedule();

   // Create event handler for keyboard
   SoEventCallback *keyboardEvent = new SoEventCallback;
   keyboardEvent->addEventCallback(SoKeyboardEvent::getClassTypeId(), keyboardfunc, &appData);
   root->addChild(keyboardEvent);

   // Create event handler for mouse
   SoEventCallback *mouseEvent = new SoEventCallback;
   mouseEvent->addEventCallback(SoMouseButtonEvent::getClassTypeId(), mousefunc, &appData);
   root->addChild(mouseEvent);

   // Create event handler for mouse motion
   SoEventCallback *motionEvent = new SoEventCallback;
   motionEvent->addEventCallback(SoLocation2Event::getClassTypeId(), motionfunc, &appData);
   root->addChild(motionEvent);

   //Start show
   viewer->setSceneGraph(root);
   viewer->show();

   SoWin::show(window);
   SoWin::mainLoop();
   delete viewer;
   root->unref();
   return 0;
}

