/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkSimilarity2DTransform_hxx
#define __itkSimilarity2DTransform_hxx

#include "itkSimilarity2DTransform.h"
#include "vnl/vnl_math.h"

namespace itk
{
// Constructor with default arguments
template< class TScalarType >
Similarity2DTransform< TScalarType >
::Similarity2DTransform():Superclass(OutputSpaceDimension, ParametersDimension)
{
  m_Scale = 1.0f;
}

// Constructor with arguments
template< class TScalarType >
Similarity2DTransform< TScalarType >::Similarity2DTransform(unsigned int spaceDimension,
                                                            unsigned int parametersDimension):
  Superclass(spaceDimension, parametersDimension)
{
  m_Scale = 1.0f;
}

// Set Parameters
template< class TScalarType >
void
Similarity2DTransform< TScalarType >
::SetParameters(const ParametersType & parameters)
{
  itkDebugMacro(<< "Setting parameters " << parameters);

  //Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &parameters != &(this->m_Parameters) )
    {
    this->m_Parameters = parameters;
    }

  // Set scale
  const TScalarType scale = parameters[0];
  this->SetVarScale(scale);

  // Set angle
  const TScalarType angle = parameters[1];
  this->SetVarAngle(angle);

  // Set translation
  OffsetType translation;
  for ( unsigned int i = 0; i < SpaceDimension; i++ )
    {
    translation[i] = parameters[i + 2];
    }
  this->SetVarTranslation(translation);

  this->ComputeMatrix();
  this->ComputeOffset();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();

  itkDebugMacro(<< "After setting parameters ");
}

// Get Parameters
template< class TScalarType >
const typename Similarity2DTransform< TScalarType >::ParametersType &
Similarity2DTransform< TScalarType >
::GetParameters(void) const
{
  itkDebugMacro(<< "Getting parameters ");

  this->m_Parameters[0] = this->GetScale();
  this->m_Parameters[1] = this->GetAngle();

  // Get the translation
  OffsetType translation = this->GetTranslation();
  for ( unsigned int i = 0; i < SpaceDimension; i++ )
    {
    this->m_Parameters[i + 2] = translation[i];
    }

  itkDebugMacro(<< "After getting parameters " << this->m_Parameters);

  return this->m_Parameters;
}

// Set Scale Part
template< class TScalarType >
void
Similarity2DTransform< TScalarType >
::SetScale(ScaleType scale)
{
  m_Scale = scale;
  this->ComputeMatrix();
  this->ComputeOffset();
}

// Compute the matrix
template< class TScalarType >
void
Similarity2DTransform< TScalarType >
::ComputeMatrix(void)
{
  const double angle = this->GetAngle();

  const double cc = vcl_cos(angle);
  const double ss = vcl_sin(angle);

  const MatrixValueType ca = cc * m_Scale;
  const MatrixValueType sa = ss * m_Scale;

  MatrixType matrix;

  matrix[0][0] = ca; matrix[0][1] = -sa;
  matrix[1][0] = sa; matrix[1][1] = ca;

  this->SetVarMatrix(matrix);
}

/** Compute the Angle from the Rotation Matrix */
template< class TScalarType >
void
Similarity2DTransform< TScalarType >
::ComputeMatrixParameters(void)
{
  m_Scale = vcl_sqrt( vnl_math_sqr(this->GetMatrix()[0][0])
                      + vnl_math_sqr(this->GetMatrix()[0][1]) );

  this->SetVarAngle( vcl_acos(this->GetMatrix()[0][0] / m_Scale) );

  if ( this->GetMatrix()[1][0] < 0.0 )
    {
    this->SetVarAngle( -this->GetAngle() );
    }

  if ( ( this->GetMatrix()[1][0] / m_Scale ) - vcl_sin( this->GetAngle() ) > 0.000001 )
    {
    std::cout << "Bad Rotation Matrix" << std::endl;
    }
}

// Compute the transformation Jacobian
template< class TScalarType >
const typename Similarity2DTransform< TScalarType >::JacobianType &
Similarity2DTransform< TScalarType >::
GetJacobian(const InputPointType & p) const
{
  GetJacobianWithRespectToParameters( p, this->m_Jacobian );
  return this->m_Jacobian;
}

template< class TScalarType >
void
Similarity2DTransform< TScalarType >
::GetJacobianWithRespectToParameters(const InputPointType & p, JacobianType & jacobian) const
{
  const double angle = this->GetAngle();
  const double ca = vcl_cos(angle);
  const double sa = vcl_sin(angle);

  jacobian.SetSize( 2, this->GetNumberOfLocalParameters() );
  jacobian.Fill(0.0);

  const InputPointType center = this->GetCenter();
  const double         cx = center[0];
  const double         cy = center[1];

  // derivatives with respect to the scale
  jacobian[0][0] =    ca * ( p[0] - cx ) - sa * ( p[1] - cy );
  jacobian[1][0] =    sa * ( p[0] - cx ) + ca * ( p[1] - cy );

  // derivatives with respect to the angle
  jacobian[0][1] = ( -sa * ( p[0] - cx ) - ca * ( p[1] - cy ) ) * m_Scale;
  jacobian[1][1] = ( ca * ( p[0] - cx ) - sa * ( p[1] - cy ) ) * m_Scale;

  // compute derivatives with respect to the translation part
  // first with respect to tx
  jacobian[0][2] = 1.0;
  jacobian[1][2] = 0.0;
  // first with respect to ty
  jacobian[0][3] = 0.0;
  jacobian[1][3] = 1.0;
}

// Set Identity
template< class TScalarType >
void
Similarity2DTransform< TScalarType >
::SetIdentity(void)
{
  this->Superclass::SetIdentity();
  m_Scale = static_cast< TScalarType >( 1.0f );
}

// Print self
template< class TScalarType >
void
Similarity2DTransform< TScalarType >::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Scale =" << m_Scale  << std::endl;
}

// Create and return an inverse transformation
template< class TScalarType >
void
Similarity2DTransform< TScalarType >::CloneInverseTo(Pointer & result) const
{
  result = New();
  this->GetInverse( result.GetPointer() );
}

// return an inverse transformation
template< class TScalarType >
bool
Similarity2DTransform< TScalarType >::GetInverse(Self *inverse) const
{
  if ( !inverse )
    {
    return false;
    }

  inverse->SetCenter( this->GetCenter() );  // inverse have the same center
  inverse->SetScale( NumericTraits< double >::One / this->GetScale() );
  inverse->SetAngle( -this->GetAngle() );
  inverse->SetTranslation( -( this->GetInverseMatrix() * this->GetTranslation() ) );

  return true;
}

// Return an inverse of this transform
template< class TScalarType >
typename Similarity2DTransform< TScalarType >::InverseTransformBasePointer
Similarity2DTransform< TScalarType >
::GetInverseTransform() const
{
  Pointer inv = New();

  if ( this->GetInverse(inv) )
    {
    return inv.GetPointer();
    }
  return NULL;
}

// Create and return a clone of the transformation
template< class TScalarType >
void
Similarity2DTransform< TScalarType >::CloneTo(Pointer & result) const
{
  result = New();
  result->SetCenter( this->GetCenter() );
  result->SetScale( this->GetScale() );
  result->SetAngle( this->GetAngle() );
  result->SetTranslation( this->GetTranslation() );
}

// Set the similarity matrix
template< class TScalarType >
void
Similarity2DTransform< TScalarType >::SetMatrix(const MatrixType & matrix)
{
  itkDebugMacro("setting  m_Matrix  to " << matrix);

  typename MatrixType::InternalMatrixType test =
    matrix.GetVnlMatrix() * matrix.GetTranspose();

  test /= test[0][0]; // factor out the scale

  const double tolerance = 1e-10;
  if ( !test.is_identity(tolerance) )
    {
    itk::ExceptionObject ex(__FILE__, __LINE__, "Attempt to set a Non-Orthogonal matrix", ITK_LOCATION);
    throw ex;
    }

  this->SetVarMatrix(matrix);
  this->ComputeOffset();
  this->ComputeMatrixParameters();
  this->Modified();
}
} // namespace

#endif
