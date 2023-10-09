#include "curses-gfx-loader.h"

#include <iostream>

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags


Mesh::Mesh()
:vi(NULL), numTriangles(0) {
}

Mesh::~Mesh() {
    if(vi) {
        delete [] vi;
    }
}

Material::Material()
: textures(NULL), numTextures(0) {
}

Material::~Material() {
    if(textures) {
        delete [] textures;
    }
}

Scene::Scene()
:meshes(NULL), materials(NULL), numMaterials(0), numMeshes(0) {
}

Scene::~Scene() {
    if(meshes) {
        delete [] meshes;
    }
    if(materials) {
        delete [] materials;
    }
}

int Scene::load(const char* file) {
    
    Assimp::Importer importer;
    std::string objectFilePath = file;
    std::string objectFileDirectory;
    const size_t last_slash_idx = objectFilePath.rfind('/');
    if (std::string::npos != last_slash_idx)
    {
        objectFileDirectory = objectFilePath.substr(0, last_slash_idx);
    }

    const aiScene* scene = importer.ReadFile( objectFilePath,
//        aiProcess_CalcTangentSpace       |
        aiProcess_Triangulate            |
//        aiProcess_JoinIdenticalVertices  |
                                    aiProcess_FlipUVs |
        aiProcess_GenSmoothNormals  |
//                                             aiProcess_MakeLeftHanded |
//                                             aiProcess_OptimizeMeshes |
//                                             aiProcess_FlipWindingOrder |
        aiProcess_SortByPType);
    
    if (nullptr == scene) {
        std::cerr << "Error loading file: " << importer.GetErrorString() << std::endl;
        return -1;
      }
//    struct Mesh {
//        MeshVertexInfo* vi;
//        int numTriangles;
//        bool hasTexture;
//        int materialId;
//    };
    
    printf("File has %d meshes\n", scene->mNumMeshes);
//    MeshVertexInfo* fileVertices = NULL;
    //int numFileEdges = 0;
    numMeshes = scene->mNumMeshes;
//    Mesh meshes[numMeshes];
    meshes = new Mesh[numMeshes];
    for(int m = 0; m < scene->mNumMeshes; m++) {
        
        printf(" - mesh[%d] has %d vertices\n", m, scene->mMeshes[m]->mNumVertices);
//        meshes[m].numTriangles = scene->mMeshes[m]->mNumVertices/3;
//        meshes[m].vi = new MeshVertexInfo[scene->mMeshes[m]->mNumVertices];
        
        printf(" - mesh[%d] has normals:   %d\n", m, scene->mMeshes[m]->HasNormals());
        printf(" - mesh[%d] has faces  :   %d\n", m, scene->mMeshes[m]->mNumFaces);
        printf(" - mesh[%d] has texcoords(0): %d\n", m, scene->mMeshes[m]->HasTextureCoords(0));
        printf(" - mesh[%d] has colors(0): %d\n", m, scene->mMeshes[m]->HasVertexColors(0));
        //        printf(" - texcooords[0][0] name: %s\n", scene->mMeshes[0]->mTextureCoordsNames[0][0].C_Str());
        
        printf(" - - Allocated!\n");
        int numTriangles = 0;
        for(int f = 0; f < scene->mMeshes[m]->mNumFaces; f++) {
            if(scene->mMeshes[m]->mFaces[f].mNumIndices == 3) {
                numTriangles++;
            }
        }
        printf(" - mesh[%d] has triangles  :   %d\n", m, numTriangles);
        meshes[m].numTriangles = numTriangles;
        meshes[m].vi = new MeshVertexInfo[numTriangles*3];
        
        int i = 0;
        for(int f = 0; f < scene->mMeshes[m]->mNumFaces; f++) {
            aiFace face = scene->mMeshes[m]->mFaces[f];
            if(face.mNumIndices == 3) {
                for(int fi = 0; fi < face.mNumIndices; fi++) {
                    meshes[m].vi[i].vertex.x = scene->mMeshes[m]->mVertices[face.mIndices[fi]].x;
                    meshes[m].vi[i].vertex.y = scene->mMeshes[m]->mVertices[face.mIndices[fi]].y;
                    meshes[m].vi[i].vertex.z = scene->mMeshes[m]->mVertices[face.mIndices[fi]].z;
                    meshes[m].vi[i].vertex.w = 1;
                    
                    meshes[m].vi[i].normal.x = scene->mMeshes[m]->mNormals[face.mIndices[fi]].x;
                    meshes[m].vi[i].normal.y = scene->mMeshes[m]->mNormals[face.mIndices[fi]].y;
                    meshes[m].vi[i].normal.z = scene->mMeshes[m]->mNormals[face.mIndices[fi]].z;
                    
                    if(scene->mMeshes[m]->HasTextureCoords(0)) {
                        meshes[m].vi[i].textureCoord.x = scene->mMeshes[m]->mTextureCoords[0][face.mIndices[fi]].x;
                        meshes[m].vi[i].textureCoord.y = scene->mMeshes[m]->mTextureCoords[0][face.mIndices[fi]].y;
                    }
                    i++;
                }
            }
        }
        
        
        meshes[m].materialId = scene->mMeshes[m]->mMaterialIndex;
        printf("meshes[%d].materialId = %d\n", m, meshes[m].materialId);
        
    }
//    struct Material {
//        Texture* texture;
//        int numTextures;
//        bool hasDiffuseColor;
//        Coordinates3D colorDiffuse;
//    };
    numMaterials = scene->mNumMaterials;
    materials = new Material[numMaterials];
//    Texture *mTextureDiffuse = new Texture[scene->mNumMaterials];   // Likely incorrect
    printf("Scene has textures: %d\n", scene->HasTextures());
    for(int i = 0; i < scene->mNumMaterials; i++) {
        aiMaterial* material = scene->mMaterials[i];
//        printf(" - Name : %s\n", scene->mMaterials[i]->GetName().C_Str());
        printf(" - Material : %d num props:%d\n", i, material->mNumProperties);

        aiPropertyTypeInfo pdas;
        for(int p = 0; p < material->mNumProperties; p++) {
//            printf(" - - Type: %d\n", material->mProperties[p]->mType);
//            printf(" - - Name: %s\n", material->mProperties[p]->mKey.C_Str());
//            printf(" - - Semantic: %d\n", material->mProperties[p]->mSemantic);
//            printf(" - - Data Length: %d\n", material->mProperties[p]->mDataLength);
//            printf(" - - mIndex %d\n", material->mProperties[p]->mIndex);
            
//            ai_real fArray[
//            switch (material->mProperties[p]->mType) {
//                case aiPropertyTypeInfo::aiPTI_Float:
//                    aiGetMaterialFloatArray(material, AI_MATKEY_COLOR_DIFFUSE, material->mProperties[p]->mType, )
//                    break;
//
//                default:
//                    break;
//            }
        }
        
        aiColor4D diffuse;
        float c[4];
        if(AI_SUCCESS == aiGetMaterialColor(material, AI_MATKEY_COLOR_DIFFUSE, &diffuse)) {
            printf(" - - - Diffuse color: (%f,%f,%f,%f)\n", diffuse.r, diffuse.g, diffuse.b, diffuse.a);
            materials[i].hasDiffuseColor = true;
            materials[i].colorDiffuse = {diffuse.r, diffuse.g, diffuse.b};
        }
        
        printf(" - Texture count : %d\n", material->GetTextureCount(aiTextureType_DIFFUSE));
//        for(int t = 0; t < scene->mMaterials[i]->GetTextureCount(aiTextureType_DIFFUSE); t++) {
//            scene->mMaterials[i]->
//        }
        materials[i].numTextures = material->GetTextureCount(aiTextureType_DIFFUSE);
        if(materials[i].numTextures > 0) {
            materials[i].textures = new Texture[materials[i].numTextures];
            //        if(material->GetTextureCount(aiTextureType_DIFFUSE)) {
            for(int t = 0; t < materials[i].numTextures; t++) {
                aiString textureName;
                material->Get(AI_MATKEY_TEXTURE(aiTextureType_DIFFUSE, t), textureName);
                std::string textureFullPath = objectFileDirectory + "/" + textureName.C_Str();
                printf(" = Texture name: %s\n", textureFullPath.c_str());
                if(materials[i].textures[t].loadPng(textureFullPath.c_str()))
                {
                    printf(" ---- Error loading texture\n");
                }
                
                printf(" Texture size: (%d,%d)\n", materials[i].textures[t].width, materials[i].textures[t].height);
                
            }
        }
        
    }
    
    printf("Loaded!\n");
    
    
    return 0;
}
