#include "curses-gfx-physics.h"


void CursesGfxPhysics::initialize() {
    ///-----initialization_start-----
    
    ///collision configuration contains default setup for memory, collision setup. Advanced users can create their own configuration.
    collisionConfiguration = new btDefaultCollisionConfiguration();
    
    ///use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
    dispatcher = new btCollisionDispatcher(collisionConfiguration);
    
    ///btDbvtBroadphase is a good general purpose broadphase. You can also try out btAxis3Sweep.
    overlappingPairCache = new btDbvtBroadphase();
    
    ///the default constraint solver. For parallel processing you can use a different solver (see Extras/BulletMultiThreaded)
    solver = new btSequentialImpulseConstraintSolver;
    
    dynamicsWorld = new btDiscreteDynamicsWorld(dispatcher, overlappingPairCache, solver, collisionConfiguration);
    
    dynamicsWorld->setGravity(btVector3(0, 0, -9.81));
    
    ///-----initialization_end-----
}


void CursesGfxPhysics::addBody() {
    btCollisionShape* groundShape = new btBoxShape(btVector3(btScalar(10.), btScalar(10.), btScalar(10.)));
    
    collisionShapes.push_back(groundShape);
    
    btTransform groundTransform;
    groundTransform.setIdentity();
    groundTransform.setOrigin(btVector3(0, 0, -10));
    
    btScalar mass(0.);
    
    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    bool isDynamic = (mass != 0.f);
    
    btVector3 localInertia(0, 0, 0);
    if (isDynamic)
        groundShape->calculateLocalInertia(mass, localInertia);
    
    //using motionstate is optional, it provides interpolation capabilities, and only synchronizes 'active' objects
    btDefaultMotionState* myMotionState = new btDefaultMotionState(groundTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, groundShape, localInertia);
    ground = new btRigidBody(rbInfo);
    
    
    //add the body to the dynamics world
    dynamicsWorld->addRigidBody(ground);
}


void CursesGfxPhysics::addSphere() {
    //create a dynamic rigidbody

    //btCollisionShape* colShape = new btBoxShape(btVector3(1,1,1));
    btCollisionShape* colShape = new btSphereShape(btScalar(1.));
    collisionShapes.push_back(colShape);
    
    /// Create Dynamic Objects
    btTransform startTransform;
    startTransform.setIdentity();
    
    btScalar mass(1.f);
    
    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    bool isDynamic = (mass != 0.f);
    
    btVector3 localInertia(0, 0, 0);
    if (isDynamic)
        colShape->calculateLocalInertia(mass, localInertia);
    
    startTransform.setOrigin(btVector3(0, 0, 30));
    
    //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
    btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
    ball = new btRigidBody(rbInfo);
    
    dynamicsWorld->addRigidBody(ball);
}

btRigidBody* CursesGfxPhysics::addCube(Coordinates3D scale, double m, Mat4D initial) {
    //create a dynamic rigidbody
//    btBvhTriangleMeshShape btms;
    
    //btCollisionShape* colShape = new btBoxShape(btVector3(1,1,1));
    btCollisionShape* colShape = new btBoxShape(btVector3(btScalar(scale.x), btScalar(scale.y), btScalar(scale.z)));
    collisionShapes.push_back(colShape);
    
    /// Create Dynamic Objects
    btTransform startTransform;
//    startTransform.setIdentity();
//    startTransform.setOrigin(btVector3(0, 0, 10));
//    startTransform.setRotation(btQuaternion(0.707, 0.707, .1));
    
    float mat[16];
//    trans.getOpenGLMatrix(mat);
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            mat[i + 4*j] = initial.d[i][j];
        }
    }
    startTransform.setFromOpenGLMatrix(mat);
    
    btScalar mass(m);
    
    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    bool isDynamic = (mass != 0.f);
    
    btVector3 localInertia(0, 0, 1);
    if (isDynamic)
        colShape->calculateLocalInertia(mass, localInertia);
    
    
    
    
    //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
    btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
    btRigidBody* cube = new btRigidBody(rbInfo);
    
    dynamicsWorld->addRigidBody(cube);
    return cube;
}


btRigidBody* CursesGfxPhysics::addUprightCapsule(Coordinates3D scale, double m, Mat4D initial) {
    //create a dynamic rigidbody
    //btCapsuleShapeZ represents a capsule around the Z axis the total height is height+2*radius, so the height is just the height between the center of each 'sphere' of the capsule caps.
    btCollisionShape* colShape = new btCapsuleShapeZ(btScalar(scale.x), btScalar(scale.z));
    collisionShapes.push_back(colShape);
    
    /// Create Dynamic Objects
    btTransform startTransform;
//    startTransform.setIdentity();
//    startTransform.setOrigin(btVector3(0, 0, 10));
//    startTransform.setRotation(btQuaternion(0.707, 0.707, .1));
    
    float mat[16];
//    trans.getOpenGLMatrix(mat);
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            mat[i + 4*j] = initial.d[i][j];
        }
    }
    startTransform.setFromOpenGLMatrix(mat);
    
    btScalar mass(m);
    
    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    bool isDynamic = (mass != 0.f);
    
    btVector3 localInertia(0, 0, 1);
    if (isDynamic)
        colShape->calculateLocalInertia(mass, localInertia);
    
    
    
    
    //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
    btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
    btRigidBody* cube = new btRigidBody(rbInfo);
    cube->setAngularFactor(btVector3(0.0f, 0.0f, 1.0f));
    
    dynamicsWorld->addRigidBody(cube);
    return cube;
    
}


void CursesGfxPhysics::update(double dTime) {
    dynamicsWorld->stepSimulation(dTime, 100);
}
