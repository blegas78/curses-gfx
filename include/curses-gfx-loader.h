#ifndef CURSES_GFX_LOADER_H
#define CURSES_GFX_LOADER_H


#include "curses-gfx-types.h"
#include "curses-gfx-texture.h"

typedef struct _MeshVertexInfo {
    Coordinates4D vertex;
    Coordinates4D location;
    Coordinates3D normal;
    Coordinates2Df textureCoord;
    ColorRGBA color;
} MeshVertexInfo;



class Mesh {
public:
    MeshVertexInfo* vi;
    int numTriangles;
    bool hasTexture;
    int materialId;
    
    Mesh();
    ~Mesh();
};

class Material {
public:
    Texture* textures;
    int numTextures;
    bool hasDiffuseColor;
    Coordinates3D colorDiffuse;
    
    Material();
    ~Material();
};

class Scene {
public:
    Mesh *meshes;
    Material *materials;
    int numMeshes;
    int numMaterials;
    
    int load(const char* file);
    Scene();
    ~Scene();
};




#endif
