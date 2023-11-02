#ifndef CURSES_GFX_PHYSICS_H
#define CURSES_GFX_PHYSICS_H

#include <bullet/btBulletDynamicsCommon.h>

#include "curses-gfx-types.h"

class CursesGfxPhysics {
private:
    
    btDefaultCollisionConfiguration* collisionConfiguration;

    ///use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
    btCollisionDispatcher* dispatcher;

    ///btDbvtBroadphase is a good general purpose broadphase. You can also try out btAxis3Sweep.
    btBroadphaseInterface* overlappingPairCache;

    ///the default constraint solver. For parallel processing you can use a different solver (see Extras/BulletMultiThreaded)
    btSequentialImpulseConstraintSolver* solver;

    btDiscreteDynamicsWorld* dynamicsWorld;
    
    //keep track of the shapes, we release memory at exit.
        //make sure to re-use collision shapes among rigid bodies whenever possible!
        btAlignedObjectArray<btCollisionShape*> collisionShapes;
    
public:
//    btCollisionShape* groundShape;
//    btCollisionShape* colShape;
    btRigidBody* ground;
    btRigidBody* ball;
//    btRigidBody* cube;
public:
    
    void initialize();
    
    void addBody();
    void addSphere();
    btRigidBody* addCube(Coordinates3D scale, double mass, Mat4D initial);
    
    void update(double dTime);
};

#endif
