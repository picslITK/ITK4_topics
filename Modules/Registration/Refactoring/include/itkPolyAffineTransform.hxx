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
#ifndef __itkPolyAffineTransform_hxx
#define __itkPolyAffineTransform_hxx

#include "itkNumericTraits.h"
#include "itkPolyAffineTransform.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "itkMath.h"
#include "itkCrossHelper.h"

namespace itk
{
// Constructor with default arguments
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::PolyAffineTransform():
  Superclass()
{
}

// Destructor
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::~PolyAffineTransform()
{
  return;
}

// Print self
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
void
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  unsigned int t;

  os << indent << "PolyAffineTransform: contains "
     << m_AtomicTransformList.size() << " transforms " << std::endl;

  for ( t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    os << indent << "Transforms[" << t << "] = ";
    os << indent.GetNextIndent();
    //m_AtomicTransformList[t]->PrintSelf(os, indent.GetNextIndent());
    os << m_AtomicTransformList[t]->GetTransformTypeAsString();
    os << std::endl;
    os << indent.GetNextIndent() << "NumberOfParameters ";
    os << m_AtomicTransformList[t]->GetNumberOfParameters();
    os << std::endl;
    }
}
// Constructor with explicit arguments
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
void
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::SetIdentity(void)
{
  for ( unsigned int t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    m_AtomicTransformList[t]->SetIdentity();
    }
  this->GetParameters();
  this->GetFixedParameters();
  this->Modified();
}

// Transform a point
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
typename PolyAffineTransform< TScalarType,
                                    NInputDimensions,
                                    NOutputDimensions >::OutputPointType
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::TransformPoint(const InputPointType & point) const
{
  if (m_AtomicTransformList.size() <= 0)
    {
    return point;
    }

  unsigned int i;
  OutputPointType outPoint, sumPoint;
  double weight, sumWeight;

  sumPoint.Fill(0.0);
  sumWeight = 0.0;

  for ( i = 0; i < m_AtomicTransformList.size(); i++ )
    {
    outPoint = m_AtomicTransformList[i]->TransformPoint( point );
    weight = m_WeightFunctorList[i]->GetWeight( point );

    sumWeight = sumWeight + weight;
    for ( unsigned d = 0; d < NOutputDimensions; d++ )
      {
      sumPoint[d] += outPoint[d] * weight;
      }
    }

  for ( unsigned d = 0; d < NOutputDimensions; d++ )
    {
    sumPoint[d] /= sumWeight;
    }

  return sumPoint;
}

// Transform a vector
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
typename PolyAffineTransform< TScalarType,
                                    NInputDimensions,
                                    NOutputDimensions >::OutputVectorType
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::TransformVector(const InputVectorType & vect) const
{
  itkExceptionMacro("TransformVector not yet implemented.");
  OutputVectorType output;
  return output;
}

// Transform a vnl_vector_fixed
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
typename PolyAffineTransform< TScalarType,
                                    NInputDimensions,
                                    NOutputDimensions >::OutputVnlVectorType
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::TransformVector(const InputVnlVectorType & vect) const
{
  //itkExceptionMacro("TransformVector not yet implemented.");
  //OutputVnlVectorType output;
  //return output;
  return vect;
}

// Transform a CovariantVector
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
typename PolyAffineTransform< TScalarType,
                                    NInputDimensions,
                                    NOutputDimensions >::OutputCovariantVectorType
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::TransformCovariantVector(const InputCovariantVectorType & vect) const
{
  //itkExceptionMacro("TransformCovariantVector not yet implemented.");
  //OutputCovariantVectorType result;     // Converted vector
  //return result;
  return vect;
}

// return an inverse transformation
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
bool
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::GetInverse(Self *inverse) const
{
  if ( !inverse )
    {
    return false;
    }

  unsigned int i;

  for ( i = 0; i < m_AtomicTransformList.size(); i++ )
    {
      if ( !m_AtomicTransformList[i]->GetInverse( inverse->GetAtomicTransformList()[i].GetPointer() ) )
        {
        return false;
        }
    }

  return true;
}

// Return an inverse of this transform
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
typename PolyAffineTransform< TScalarType, NInputDimensions,
                                    NOutputDimensions >::InverseTransformBasePointer
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::GetInverseTransform() const
{
  Pointer inv = New();

  return GetInverse(inv) ? inv.GetPointer() : NULL;
}

// Get fixed parameters
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
void
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::SetFixedParameters(const ParametersType & fp)
{
  unsigned int par = 0;

  //Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &fp != &(this->m_FixedParameters) )
    {
    this->m_FixedParameters = fp;
    }

  unsigned int paramSize;

  for ( unsigned int t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    //assuming each transform object has m_FixedParameters set initially
    AtomParametersType param = m_AtomicTransformList[t]->GetFixedParameters();
    paramSize = param.GetSize();

    for ( unsigned int par1 = 0; par1 < paramSize; par1++ )
      {
      param[par1] = this->m_FixedParameters[par];
      par++;
      }

    m_AtomicTransformList[t]->SetFixedParameters(param);
    }

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();
}

/** Get the Fixed Parameters. */
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
const typename PolyAffineTransform< TScalarType,
                                          NInputDimensions,
                                          NOutputDimensions >::ParametersType &
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::GetFixedParameters(void) const
{
  unsigned int paramSize = 0, par = 0;

  for ( unsigned int t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    paramSize += m_AtomicTransformList[t]->GetFixedParameters().GetSize();
    }

  this->m_FixedParameters.SetSize (paramSize);

  for ( unsigned int t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    const AtomParametersType &param = m_AtomicTransformList[t]->GetFixedParameters();
    for ( unsigned int par1 = 0; par1 < param.GetSize(); par1++ )
      {
      this->m_FixedParameters[par] = param[par1];
      par++;
      }
    }

  return this->m_FixedParameters;
}

// Get parameters
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
const typename PolyAffineTransform< TScalarType,
                                          NInputDimensions,
                                          NOutputDimensions >::ParametersType &
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::GetParameters(void) const
{
  unsigned int paramSize = 0, par = 0;

  for ( unsigned int t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    paramSize += m_AtomicTransformList[t]->GetNumberOfParameters();
    }

  this->m_Parameters.SetSize (paramSize);

  for ( unsigned int t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    const AtomParametersType &param = m_AtomicTransformList[t]->GetParameters();
    for ( unsigned int par1 = 0; par1 < param.GetSize(); par1++ )
      {
      this->m_Parameters[par] = param[par1];
      par++;
      }
    }

  return this->m_Parameters;
}

// Set parameters
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
void
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::SetParameters(const ParametersType & parameters)
{
  unsigned int par = 0;

  //Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &parameters != &(this->m_Parameters) )
    {
    this->m_Parameters = parameters;
    }

  unsigned int paramSize;

  for ( unsigned int t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    //assuming each transform object has m_Parameters set initially
    AtomParametersType param = m_AtomicTransformList[t]->GetParameters();
    paramSize = param.GetSize();

    for ( unsigned int par1 = 0; par1 < paramSize; par1++ )
      {
      param[par1] = this->m_Parameters[par];
      par++;
      }

    m_AtomicTransformList[t]->SetParameters(param);
    }

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();
}

// Compute the Jacobian in one position, without setting values to m_Jacobian
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
void
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::ComputeJacobianWithRespectToParameters(const InputPointType & p, JacobianType &j) const
{
  //This will not reallocate memory if the dimensions are equal
  // to the matrix's current dimensions.

  j.SetSize( NOutputDimensions, this->GetNumberOfLocalParameters() );
  j.Fill(0.0);

  unsigned int t, d, par1, par = 0;
  double *weights, sumWeight;

  sumWeight = 0;
  weights = new double[ m_AtomicTransformList.size() ];

  for ( t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    weights[t] = m_WeightFunctorList[t]->GetWeight( p );
    sumWeight = sumWeight + weights[t];
    }
  for ( t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    weights[t] /= sumWeight; //normalizing weights
    }

  for ( t = 0; t < m_AtomicTransformList.size(); t++ )
    {
    typename AtomicTransformType::JacobianType atomicJacobian;
    m_AtomicTransformList[t]->ComputeJacobianWithRespectToParameters( p, atomicJacobian );
    for ( par1 = 0; par1 < m_AtomicTransformList[t]->GetNumberOfLocalParameters(); par1++ )
      {
      for ( d = 0; d < NOutputDimensions; d++ )
        {
        j(d, par) = weights[t] * atomicJacobian[d][par1];
        }
      par++;
      }
    }

  delete weights;

  return;

}

// Output the deformation field as an image
template< class TScalarType, unsigned int NInputDimensions,
          unsigned int NOutputDimensions >
template< class TDisplacementField >
void
PolyAffineTransform< TScalarType, NInputDimensions, NOutputDimensions >
::OutputDisplacementField(typename TDisplacementField::Pointer DisplacementField) const
{
  typedef typename TDisplacementField::IndexType      IndexType;

  InputPointType point, sumPoint;
  typename TDisplacementField::PixelType displacement;

  ImageRegionIteratorWithIndex< TDisplacementField >
      iter( DisplacementField, DisplacementField->GetLargestPossibleRegion() );

  /* Iterate over the region */
  iter.GoToBegin();
  while( !iter.IsAtEnd() )
    {
    const IndexType &index = iter.GetIndex();
    DisplacementField->TransformIndexToPhysicalPoint(index, point);

    sumPoint = this->TransformPoint( point );
    displacement = sumPoint - point;
    iter.Set( displacement );

    ++iter;
    }
}

} // namespace

#endif
