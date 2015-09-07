/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode QuadricDecimationMesh::QuadricIsoSurfaces = NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize()
{
  // Allocate memory for the quadric array
  unsigned int numVerts = mVerts.size();
  mQuadrics.reserve(numVerts);
  std::streamsize width = std::cerr.precision(); // store stream precision
  for (unsigned int i = 0; i < numVerts; i++) {

    // Compute quadric for vertex i here
    mQuadrics.push_back(createQuadricForVert(i));


    // Calculate initial error, should be numerically close to 0

    Vector3<float> v0 = mVerts[i].pos;
    Vector4<float> v(v0[0],v0[1],v0[2],1);
    Matrix4x4<float> m = mQuadrics.back();

    float error = v*(m*v);
    //std::cerr << std::scientific << std::setprecision(2) << error << " ";
  }
  std::cerr << std::setprecision(width) << std::fixed; // reset stream precision

  // Run the initialize for the parent class to initialize the edge collapses
  DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute, DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(EdgeCollapse * collapse)
{
  // Compute collapse->position and collapse->cost here
  // based on the quadrics at the edge endpoints
  
  unsigned int v0 = mEdges[collapse->halfEdge].vert;
  unsigned int v1 =  mEdges[mEdges[collapse->halfEdge].pair].vert;
  Matrix4x4<float> Q0 = createQuadricForVert(v0);
  Matrix4x4<float> Q1 = createQuadricForVert(v1);

  Matrix4x4<float> Q, Qv;
  Q = Qv = Q0 + Q1;

  float* q = Qv.GetArrayPtr();
  q[3*4 + 0] = 0;    
  q[3*4 + 1] = 0;    
  q[3*4 + 2] = 0;    
  q[3*4 + 3] = 1;    

  Vector4<float> vec_base = Vector4<float>(0,0,0,1);

  Vector4<float> vec_res;

  if (Qv.IsSingular()) {
	vec_res = Qv.Inverse() * vec_base;
  } else {
	Vector3<float> temp = (mVerts[v0].pos + mVerts[v1].pos)*0.5f; 
	vec_res = Vector4<float>(temp[0], temp[1], temp[2], 1);
  }
	Vector4<float> temp = Q * vec_res;

  float cost = vec_res * temp;

  collapse->position = Vector3<float>(vec_res[0], vec_res[1], vec_res[2]);

  //Lab2.2** (VG) modify cost

  Vector3<float> f_norm = FaceNormal(mEdges[collapse->halfEdge].face); 
  
  Vector3<float> cam_vec = Vector3<float>(0,0,10) - collapse->position;

  float heuristic = f_norm * cam_vec.Normalize();

  collapse->cost = cost * (1 + heuristic);  

  //std::cerr << "computeCollapse in QuadricDecimationMesh not implemented.\n";
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(unsigned int ind)
{
  DecimationMesh::updateVertexProperties(ind);
  mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
Matrix4x4<float> QuadricDecimationMesh::createQuadricForVert(unsigned int indx) const{
  float q[4][4] = {{0,0,0,0},
                   {0,0,0,0},
                   {0,0,0,0},
                   {0,0,0,0}};
  Matrix4x4<float> Q(q);
  std::vector<unsigned int> f_neighbors = FindNeighborFaces(indx);
  
  for (unsigned int i = 0; i < f_neighbors.size(); i++) {
    Q += createQuadricForFace(f_neighbors[i]);
  }

  // The quadric for a vertex is the sum of all the quadrics for the adjacent faces
  // Tip: Matrix4x4 has an operator +=
  return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
Matrix4x4<float> QuadricDecimationMesh::createQuadricForFace(unsigned int indx) const{

  Vertex f_vert = mVerts [ mEdges[ mFaces[indx].edge ].vert ];

  Vector3<float> f_norm =  FaceNormal(indx);

  float a = f_norm[0];
  float b = f_norm[1];
  float c = f_norm[2];
  float d = -(f_norm * f_vert.pos);
  
  // Calculate the quadric (outer product of plane parameters) for a face
  // here using the formula from Garland and Heckbert

  float kp[4][4] = {{a*a,a*b,a*c,a*d},
                    {b*a,b*b,b*c,b*d},
                    {c*a,c*b,c*c,c*d},
                    {d*a,d*b,d*c,d*d}};

  Matrix4x4<float> f_Q(kp);

  return f_Q;
}


void QuadricDecimationMesh::Render()
{
  DecimationMesh::Render();

  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  if (mVisualizationMode == QuadricIsoSurfaces)
  {
    // Implement the quadric visualization here

    GLUquadric* quadratic = gluNewQuadric();
    if( !quadratic){
      std::cerr << "Cannot initialize quadartic n";
    }
    gluQuadricNormals(quadratic, GLU_SMOOTH);
    //gluQuadricDrawStyle( quadratic, GLU_FILL);
    gluQuadricTexture(quadratic, GL_TRUE);

    float err_max = 0.0f;

    for(unsigned int index = 0; index < mVerts.size(); index++) {

      Vertex vert = mVerts[index];

      Matrix4x4<float> Qv = mQuadrics[index];

	  /*
	  float* q = Qv.GetArrayPtr();
	  q[3*4 + 0] = 0;    
	  q[3*4 + 1] = 0;    
	  q[3*4 + 2] = 0;    
	  q[3*4 + 3] = 1;    
	  */

      Vector4<float> position = Vector4<float>(vert.pos[0], vert.pos[1], vert.pos[2], 1);
      Vector4<float> temp = Qv * position;

      float error = position * temp;
      err_max = std::max(error, err_max);

      if(error < 0.f){
        continue;
      }

      //std::cout << "Error is: " << error << std::endl;

      Matrix4x4<float> R;
      if(!Qv.CholeskyFactorization(R)) {
        continue;
      }
      R = R.Inverse();
      
      //std::cout << R << "\n";
      
      // Apply transform
      glPushMatrix(); // Push modelview matrix onto stack
      
      Matrix4x4<float> ogl_mat = R.ToGLMatrix();

      glMultMatrixf( ogl_mat.GetArrayPtr() );

      //glTranslatef( position[0], position[1], position[2]); 
      gluSphere( quadratic, error, 8, 8);

      // Restore modelview matrix
      glPopMatrix();

    }

    gluDeleteQuadric(quadratic); 
    std::cout << "err_max = " << err_max << "\n";

  }
}

