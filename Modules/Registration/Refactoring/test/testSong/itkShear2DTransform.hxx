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
#ifndef __itkShear2DTransform_hxx
#define __itkShear2DTransform_hxx

#include "itkShear2DTransform.h"
#include "itkMath.h"

namespace itk
{
// Constructor with default arguments
template< class ScalarType, unsigned int NDimensions >
Shear2DTransform< ScalarType, NDimensions >::Shear2DTransform():Superclass(SpaceDimension, ParametersDimension)
{
    // m_Scale.Fill(NumericTraits< ScalarType >::One);
    m_K = NumericTraits< ScalarType >::One;
    m_Center.Fill(NumericTraits< ScalarType >::Zero);
}

// Destructor
template< class ScalarType, unsigned int NDimensions >
Shear2DTransform< ScalarType, NDimensions >::
~Shear2DTransform()
{
    return;
}

// Set the parameters
template< class ScalarType, unsigned int NDimensions >
void
Shear2DTransform< ScalarType, NDimensions >
::SetParameters(const ParametersType & parameters) {
    typedef typename ParametersType::ValueType ParameterValueType;
    //  for ( unsigned int i = 0; i < SpaceDimension; i++ )
    //    {
    //    m_Scale[i] = parameters[i];
    //    }

    m_K = parameters[0];
    this->m_Parameters = parameters;
    this->ComputeMatrix();

    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();
}

// Get Parameters
template< class TScalarType, unsigned int NDimensions >
const typename Shear2DTransform< TScalarType, NDimensions >::ParametersType &
Shear2DTransform< TScalarType, NDimensions >
::GetParameters(void) const
 {
    itkDebugMacro(<< "Getting parameters ");

    //  // Transfer the translation part
    //  for ( unsigned int i = 0; i < SpaceDimension; i++ )
    //    {
    //    this->m_Parameters[i] = m_Scale[i];
    //    }

    this->m_Parameters[0] = m_K;

    itkDebugMacro(<< "After getting parameters " << this->m_Parameters);

    return this->m_Parameters;
 }

// Print self
template< class ScalarType, unsigned int NDimensions >
void
Shear2DTransform< ScalarType, NDimensions >
::PrintSelf(std::ostream & os, Indent indent) const
 {
    Superclass::PrintSelf(os, indent);

    os << indent << "K: " << m_K << std::endl;
    os << indent << "Center: " << m_Center << std::endl;
 }

// Compose with another affine transformation
template< class ScalarType, unsigned int NDimensions >
void
Shear2DTransform< ScalarType, NDimensions >::Compose(const Self *other, bool)
{
    //  for ( unsigned int i = 0; i < SpaceDimension; i++ )
    //    {
    //    m_Scale[i] *= other->m_Scale[i];
    //    }
    m_K *= other->m_K;
    return;
}

//// Compose with a scale
//template< class ScalarType, unsigned int NDimensions >
//void
//Shear2DTransform< ScalarType, NDimensions >::Scale(const ScaleType & scale, bool)
//{
//  for ( unsigned int i = 0; i < SpaceDimension; i++ )
//    {
//    m_Scale[i] *= scale[i];
//    }
//  return;
//}

// Transform a point
template< class ScalarType, unsigned int NDimensions >
typename Shear2DTransform< ScalarType, NDimensions >::OutputPointType
Shear2DTransform< ScalarType, NDimensions >::TransformPoint(const InputPointType & point) const
{
    OutputPointType result;

    //  for ( unsigned int i = 0; i < SpaceDimension; i++ )
    //    {
    //    result[i] = ( point[i] - m_Center[i] ) * m_Scale[i] + m_Center[i];
    //    }

    result[0] = ( point[0] - m_Center[0] ) + ( point[1] - m_Center[1] ) * m_K;
    result[1] = ( point[1] - m_Center[1] );

    return result;
}

// Transform a vector
template< class ScalarType, unsigned int NDimensions >
typename Shear2DTransform< ScalarType, NDimensions >::OutputVectorType
Shear2DTransform< ScalarType, NDimensions >::TransformVector(const InputVectorType & vect) const
{
  OutputVectorType result;

  result.Fill(0);

//  for ( unsigned int i = 0; i < SpaceDimension; i++ )
//    {
//    result[i] = vect[i] * m_Scale[i];
//    }
  return result;
}

// Transform a vnl_vector_fixed
template< class ScalarType, unsigned int NDimensions >
typename Shear2DTransform< ScalarType, NDimensions >::OutputVnlVectorType
Shear2DTransform< ScalarType, NDimensions >::TransformVector(const InputVnlVectorType & vect) const
{
  OutputVnlVectorType result;
  result.fill(0);

//  for ( unsigned int i = 0; i < SpaceDimension; i++ )
//    {
//    result[i] = vect[i] * m_Scale[i];
//    }
  return result;
}

// Transform a CovariantVector
template< class ScalarType, unsigned int NDimensions >
typename Shear2DTransform< ScalarType, NDimensions >::OutputCovariantVectorType
Shear2DTransform< ScalarType, NDimensions >::TransformCovariantVector(const InputCovariantVectorType & vect) const
{
  // Covariant Vectors are scaled by the inverse
  OutputCovariantVectorType result;

//  for ( unsigned int i = 0; i < SpaceDimension; i++ )
//    {
//    result[i] = vect[i] / m_Scale[i];
//    }
  result.Fill(0);
  return result;
}

// Create and return an inverse transformation
template< class ScalarType, unsigned int NDimensions >
bool
Shear2DTransform< ScalarType, NDimensions >::GetInverse(Self *inverse) const
{
    if ( !inverse )
    {
        return false;
    }

    //  for ( unsigned int i = 0; i < SpaceDimension; i++ )
    //    {
    //    inverse->m_Scale[i] = NumericTraits< double >::One / m_Scale[i];
    //    }

    inverse->m_K = -1 * NumericTraits< double >::One / m_K;

    return true;
}

// Return an inverse of this transform
template< class ScalarType, unsigned int NDimensions >
typename Shear2DTransform< ScalarType, NDimensions >::InverseTransformBasePointer
Shear2DTransform< ScalarType, NDimensions >
::GetInverseTransform() const
 {
    Pointer inv = New();

    if ( this->GetInverse(inv) )
    {
        return inv.GetPointer();
    }
    return NULL;
 }

// Compute the Jacobian of the transformation
// It follows the same order of Parameters vector
template< class ScalarType, unsigned int NDimensions >
const typename Shear2DTransform< ScalarType, NDimensions >::JacobianType &
Shear2DTransform< ScalarType, NDimensions >
::GetJacobian(const InputPointType & p) const
 {
    GetJacobianWithRespectToParameters( p, this->m_Jacobian );
    return this->m_Jacobian;
 }

// Compute the Jacobian of the transformation
// It follows the same order of Parameters vector
template< class ScalarType, unsigned int NDimensions >
void
Shear2DTransform< ScalarType, NDimensions >
::GetJacobianWithRespectToParameters(const InputPointType & p, JacobianType &j) const
 {
    j.SetSize( 2, 1 );
    j.Fill(0);
    //  for ( unsigned int dim = 0; dim < SpaceDimension; dim++ )
    //    {
    //    j(dim, dim) = p[dim] - m_Center[dim];
    //    }
    j(0, 0) = p[1] - m_Center[1];
    j(1, 0) = 0;
 }

template< class ScalarType, unsigned int NDimensions >
const typename Shear2DTransform< ScalarType, NDimensions >::ParametersType &
Shear2DTransform< ScalarType, NDimensions >
::GetFixedParameters(void) const
 {
    m_FixedParameters.SetSize(0);
    return m_FixedParameters;
 }

template< class ScalarType, unsigned int NDimensions >
void
Shear2DTransform< ScalarType, NDimensions >
::ComputeMatrix(void)
 {

    MatrixType matrix;
    matrix.SetIdentity();
    matrix[0][1] = m_K;

    this->SetVarMatrix(matrix);

 }


} // namespace

#endif
