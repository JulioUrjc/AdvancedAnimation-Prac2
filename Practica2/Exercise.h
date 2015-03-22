#pragma once;

#include <vector>
using namespace std;

#include "Vector3.h"
#include "Matrix33.h"

#include "RigidBody.h"

void AdvanceTimeStep(vector<RigidBody>& bodies, float step, bool constraint, bool collisions, float floor);
