
#include "UniformCubicSplineSubdivisionCurve.h"

UniformCubicSplineSubdivisionCurve::UniformCubicSplineSubdivisionCurve(const std::vector<Vector3<float> > &joints,
                                                                       Vector3<float> lineColor,
                                                                       float lineWidth)
  : mCoefficients(joints), mControlPolygon(joints)
{
  this->mLineColor = lineColor;
  this->mLineWidth = lineWidth;
}


void UniformCubicSplineSubdivisionCurve::Subdivide()
{
  // Allocate space for new coefficients
  std::vector<Vector3<float> > newc;

  assert(mCoefficients.size() > 4 && "Need at least 5 points to subdivide"); 

  std::cout << "start size of mCoeff is: " << mCoefficients.size() << std::endl; 

  // Implement the subdivision scheme for a natural cubic spline here

  // for i = mCoefficients.size:
  // make 2 new vector
    // modify vector 1 by c'i
  //    modify vector 2 by c'i+1/2
  //  push back v1
  //  push back v2

  newc.push_back(mCoefficients.front());

  int res_size = (mCoefficients.size()* 2 -1);

  for( int i = 1; i < res_size -1; i++) {
    Vector3<float> v;

    if( i % 2 ){
      v = 1/8.0f * (4.0f* mCoefficients.at(i/2) + 4.0f * mCoefficients.at(i/2+1));
    } else {
      v = 1/8.0f * (1.0f* mCoefficients.at(i/2-1) + 6.0f * mCoefficients.at(i/2) + 1.0f *  mCoefficients.at(i/2+1));
    }

    newc.push_back(v);
  }
  
  newc.push_back(mCoefficients.back());
  
  std::cerr << "Size: " << newc.size() << " and it should be " << (mCoefficients.size()* 2 -1) << std::endl;  

  // If 'mCoefficients' had size N, how large should 'newc' be? Perform a check here!
  assert( (newc.size() == res_size)  && "Incorrect number of new coefficients!");

 mCoefficients = newc;
}


void UniformCubicSplineSubdivisionCurve::Render()
{
  // Apply transform
  glPushMatrix(); // Push modelview matrix onto stack

  // Convert transform-matrix to format matching GL matrix format
  // Load transform into modelview matrix
  glMultMatrixf( mTransform.ToGLMatrix().GetArrayPtr() );

  mControlPolygon.Render();

  // save line point and color states
  glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_CURRENT_BIT);

  // draw segments
  glLineWidth(mLineWidth);
  glColor3fv(mLineColor.GetArrayPtr());
  glBegin(GL_LINE_STRIP);
  // just draw the spline as a series of connected linear segments
  for(unsigned int i = 0; i < mCoefficients.size(); i++){
    glVertex3fv( mCoefficients.at(i).GetArrayPtr() );
  }
  glEnd();

  // restore attribs
  glPopAttrib();

  glPopMatrix();

  GLObject::Render();
}

