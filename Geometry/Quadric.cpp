/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "Quadric.h"

Quadric::Quadric(const Matrix4x4<float> & q){
  this->mQuadric = q;
}

Quadric::~Quadric(){}

/*!
 * Evaluation of world coordinates are done through either transformation
 * of the world-coordinates by mWorld2Obj, or transformation of the quadric
 * coefficient matrix by GetTransform() ONCE (see Section 2.2 in lab text).
 */
float Quadric::GetValue(float x, float y, float z) const
{
  Matrix4x4<float> obj_mat = mWorld2Obj;
  //mWorld2Obj(obj_mat);

  Vector4<float> position = Vector4<float>(x, y, z, 1);
  position = obj_mat * position; 

  Matrix4x4<float> Q = mQuadric;

  Vector4<float> temp = Q * position;

  float result = position * temp;
  return result;
}

/*!
 * Use the quadric matrix to evaluate the gradient.
 */
Vector3<float> Quadric::GetGradient(float x, float y, float z) const
{
  Matrix4x4<float> Q_sub = mQuadric;
  Q_sub(3,0) = 0; 
  Q_sub(3,1) = 0; 
  Q_sub(3,2) = 0; 
  Q_sub(3,3) = 0; 

  Matrix4x4<float> obj_mat = mWorld2Obj;
  //mWorld2Obj(obj_mat);

  Vector4<float> position = Vector4<float>(x, y, z, 1);
  position = obj_mat * position; 

  Vector4<float> result = Q_sub * position; 

  result = 2.0f * result;

  return Vector3<float>(result[0],result[1],result[2]);
}

