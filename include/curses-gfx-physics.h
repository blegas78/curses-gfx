#ifndef CURSES_GFX_PHYSICS_H
#define CURSES_GFX_PHYSICS_H

#include <bullet/btBulletDynamicsCommon.h>

#include "curses-gfx-types.h"

class CursesGfxPhysics {
public:
    
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
    btRigidBody* addUprightCapsule(Coordinates3D scale, double mass, Mat4D initial);
    template <class T> btRigidBody* addStaticMesh(T* vertexInfo, int numTriangles);
    
    void update(double dTime);
};


template <class T> btRigidBody* CursesGfxPhysics::addStaticMesh(T* vertexInfo, int numTriangles) {
    //create a dynamic rigidbody
//    btBvhTriangleMeshShape btms;
    btTriangleMesh* triangleMeshTerrain = new btTriangleMesh();
    btVector3 vertex1, vertex2, vertex3;
    for(int i = 0; i < numTriangles; i++) {
        vertex1 = btVector3(vertexInfo[0 + i*3].vertex.x, vertexInfo[0 + i*3].vertex.y, vertexInfo[0 + i*3].vertex.z);
        vertex2 = btVector3(vertexInfo[1 + i*3].vertex.x, vertexInfo[1 + i*3].vertex.y, vertexInfo[1 + i*3].vertex.z);
        vertex3 = btVector3(vertexInfo[2 + i*3].vertex.x, vertexInfo[2 + i*3].vertex.y, vertexInfo[2 + i*3].vertex.z);
        
        triangleMeshTerrain->addTriangle(vertex1, vertex2, vertex3);
    }
    
    //btCollisionShape* colShape = new btBoxShape(btVector3(1,1,1));
//    btCollisionShape* colShape = new btBoxShape(btVector3(btScalar(scale.x), btScalar(scale.y), btScalar(scale.z)));
    btCollisionShape* collisionShapeTerrain = new btBvhTriangleMeshShape(triangleMeshTerrain, true);

    collisionShapes.push_back(collisionShapeTerrain);
    
    /// Create Dynamic Objects
//    btTransform startTransform;
//    startTransform.setIdentity();
//    startTransform.setOrigin(btVector3(0, 0, 10));
//    startTransform.setRotation(btQuaternion(0.707, 0.707, .1));
//
//    float mat[16];
////    trans.getOpenGLMatrix(mat);
//    for(int i = 0; i < 4; i++) {
//        for(int j = 0; j < 4; j++) {
//            mat[i + 4*j] = initial.d[i][j];
//        }
//    }
//    startTransform.setFromOpenGLMatrix(mat);
//
//    btScalar mass(m);
//
    //rigidbody is dynamic if and only if mass is non zero, otherwise static
//    bool isDynamic = (mass != 0.f);
//
//    btVector3 localInertia(0, 0, 1);
//    if (isDynamic)
//        colShape->calculateLocalInertia(mass, localInertia);
    
    
    
    
    //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
//    btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
//    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
//    btRigidBody* cube = new btRigidBody(rbInfo);
    
    btDefaultMotionState* motionState = new btDefaultMotionState(btTransform(btQuaternion(1.0/sqrt(2), 0, 0, 1/sqrt(2)), btVector3(0, 0, 0)));

    btRigidBody::btRigidBodyConstructionInfo rigidBodyConstructionInfo(0.0f, motionState, collisionShapeTerrain, btVector3(0, 0, 0));
    btRigidBody* rigidBodyTerrain = new btRigidBody(rigidBodyConstructionInfo);
    rigidBodyTerrain->setFriction(btScalar(0.9));
    
    dynamicsWorld->addRigidBody(rigidBodyTerrain);
    return rigidBodyTerrain;
}

#endif
