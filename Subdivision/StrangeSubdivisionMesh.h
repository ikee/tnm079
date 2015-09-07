#ifndef _strange_dubdivmesh_
#define _strange_dubdivmesh_

#include "AdaptiveLoopSubdivisionMesh.h"

class StrangeSubdivisionMesh : public AdaptiveLoopSubdivisionMesh
{
public:
  virtual void Subdivide() {
    // ....
    AdaptiveLoopSubdivisionMesh::Subdivide();
  }

protected:
  bool Subdividable(unsigned int fi){
    // Every 4th face is not subdividable - kinda strange!
    // Do something more interesting...
    //

    
    Vector3<float> pos1, pos2, pos3, facePosition;
    EdgeIterator eit = GetEdgeIterator( f(fi).edge );
    pos1 = v(eit.GetEdgeVertexIndex()).pos; eit.Next();
    pos2 = v(eit.GetEdgeVertexIndex()).pos; eit.Next();
    pos3 = v(eit.GetEdgeVertexIndex()).pos;

    facePosition = pos1 + pos2;
    facePosition = facePosition + pos3;
    facePosition = facePosition/3;


    Vector3<float> cam_vec = Vector3<float>(0,0,10) - facePosition;
    
    Vector3<float> f_norm = FaceNormal( fi );
    
    float result = f_norm * cam_vec.Normalize();

    //result = FaceCurvature(fi);
    //std::cout << result << std::endl;
    return (result > 0.0f);
  }

};

#endif
