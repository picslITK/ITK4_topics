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
#ifndef __itkTimeVaryingVelocityFieldIntegrationImageFilter_hxx
#define __itkTimeVaryingVelocityFieldIntegrationImageFilter_hxx

#include "itkTimeVaryingVelocityFieldIntegrationImageFilter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorLinearInterpolateImageFunction.h"

namespace itk
{

/*
 * TimeVaryingVelocityFieldIntegrationImageFilter class definitions
 */
template<class TTimeVaryingVelocityField, class TDeformationField>
TimeVaryingVelocityFieldIntegrationImageFilter
  <TTimeVaryingVelocityField, TDeformationField>
::TimeVaryingVelocityFieldIntegrationImageFilter() :
  m_LowerTimeBound( 0.0 ),
  m_UpperTimeBound( 1.0 ),
  m_NumberOfIntegrationSteps( 10 )
{
  this->SetNumberOfRequiredInputs( 1 );

  if( InputImageDimension - 1 != OutputImageDimension )
    {
    itkExceptionMacro( "The time-varying velocity field (input) should have "
      << "dimensionality of 1 greater than the deformation field (output). " );
    }

  typedef VectorLinearInterpolateImageFunction<TimeVaryingVelocityFieldType,
    RealType> DefaultVelocityFieldInterpolatorType;
  typename DefaultVelocityFieldInterpolatorType::Pointer
    velocityFieldInterpolator = DefaultVelocityFieldInterpolatorType::New();
  this->m_VelocityFieldInterpolator = velocityFieldInterpolator;

  typedef VectorLinearInterpolateImageFunction<DeformationFieldType,
    RealType> DefaultDeformationFieldInterpolatorType;
  typename DefaultDeformationFieldInterpolatorType::Pointer
    deformationFieldInterpolator =
    DefaultDeformationFieldInterpolatorType::New();
  this->m_DeformationFieldInterpolator = deformationFieldInterpolator;
}

template<class TTimeVaryingVelocityField, class TDeformationField>
TimeVaryingVelocityFieldIntegrationImageFilter
  <TTimeVaryingVelocityField, TDeformationField>
::~TimeVaryingVelocityFieldIntegrationImageFilter()
{
}

template<class TTimeVaryingVelocityField, class TDeformationField>
void
TimeVaryingVelocityFieldIntegrationImageFilter
  <TTimeVaryingVelocityField, TDeformationField>
::GenerateOutputInformation()
{
  typename TimeVaryingVelocityFieldType::ConstPointer input = this->GetInput();
  typename DeformationFieldType::Pointer output = this->GetOutput();

  if( !input || !output )
    {
    return;
    }

  typename DeformationFieldType::SizeType size;
  typename DeformationFieldType::SpacingType spacing;
  typename DeformationFieldType::PointType origin;
  typename DeformationFieldType::DirectionType direction;

  for( unsigned int i = 0; i < OutputImageDimension; i++ )
    {
    size[i] = input->GetRequestedRegion().GetSize()[i];
    spacing[i] = input->GetSpacing()[i];
    origin[i] = input->GetOrigin()[i];

    for( unsigned int j = 0; j < OutputImageDimension; j++ )
      {
      direction[i][j] = input->GetDirection()[i][j];
      }
    }

  output->SetOrigin( origin );
  output->SetSpacing( spacing );
  output->SetDirection( direction );
  output->SetRegions( size );
}

template<class TTimeVaryingVelocityField, class TDeformationField>
void
TimeVaryingVelocityFieldIntegrationImageFilter
  <TTimeVaryingVelocityField, TDeformationField>
::BeforeThreadedGenerateData()
{
  VectorType zeroVector( 0.0 );

  this->AllocateOutputs();
  this->GetOutput()->FillBuffer( zeroVector );

  this->m_VelocityFieldInterpolator->SetInputImage( this->GetInput() );

  if( !this->m_InitialDiffeomorphism.IsNull() )
    {
    this->m_DeformationFieldInterpolator->SetInputImage(
      this->m_InitialDiffeomorphism );
    }
}

template<class TTimeVaryingVelocityField, class TDeformationField>
void
TimeVaryingVelocityFieldIntegrationImageFilter
  <TTimeVaryingVelocityField, TDeformationField>
::ThreadedGenerateData( const OutputRegionType &region,
  ThreadIdType itkNotUsed( threadId ) )
{
  if( this->m_LowerTimeBound == this->m_UpperTimeBound ||
    this->m_NumberOfIntegrationSteps == 0 )
    {
    return;
    }

  ImageRegionIteratorWithIndex<DeformationFieldType> It( this->GetOutput(),
    region );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    PointType point;
    this->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );
    VectorType displacement = this->IntegrateVelocityAtPoint( point );
    It.Set( displacement );
    }
}

template<class TTimeVaryingVelocityField, class TDeformationField>
typename TimeVaryingVelocityFieldIntegrationImageFilter
  <TTimeVaryingVelocityField, TDeformationField>::VectorType
TimeVaryingVelocityFieldIntegrationImageFilter
  <TTimeVaryingVelocityField, TDeformationField>
::IntegrateVelocityAtPoint( const PointType initialSpatialPoint )
{
  // Solve the initial value problem using fourth-order Runge-Kutta
  //    y' = f(t, y), y(t_0) = y_0

  VectorType zeroVector;
  zeroVector.Fill( 0.0 );

  // Initial conditions

  PointType spatialPoint = initialSpatialPoint;
  if( !this->m_InitialDiffeomorphism.IsNull() )
    {
    if( this->m_DeformationFieldInterpolator->IsInsideBuffer( spatialPoint ) )
      {
      spatialPoint +=
        this->m_DeformationFieldInterpolator->Evaluate( spatialPoint );
      }
    }

  // Perform the integration

  // Need to know how to map the time dimension of the input image to the
  // assumed domain of [0,1].

  RealType timePoint = this->GetInput()->GetOrigin()[InputImageDimension-1];
  RealType timeDomain = this->GetInput()->GetSpacing()[InputImageDimension-1] *
    static_cast<RealType>( this->GetInput()->
    GetLargestPossibleRegion().GetSize()[InputImageDimension-1] - 1 );

  timePoint += this->m_LowerTimeBound * timeDomain;

  // Calculate the delta time used for integration
  RealType deltaTime = ( this->m_UpperTimeBound -
    this->m_LowerTimeBound ) / static_cast<RealType>(
    this->m_NumberOfIntegrationSteps );

  if( timeDomain == 0.0 )
    {
    return zeroVector;
    }

  for( unsigned int n = 0; n < this->m_NumberOfIntegrationSteps; n++ )
    {
    typename TimeVaryingVelocityFieldType::PointType x1;
    for( unsigned int d = 0; d < OutputImageDimension; d++ )
      {
      x1[d] = spatialPoint[d];
      }
    x1[InputImageDimension-1] = timePoint;
    VectorType f1 = zeroVector;
    if( this->m_VelocityFieldInterpolator->IsInsideBuffer( x1 ) )
      {
      f1 = this->m_VelocityFieldInterpolator->Evaluate( x1 );
      }

    typename TimeVaryingVelocityFieldType::PointType x2;
    for( unsigned int d = 0; d < OutputImageDimension; d++ )
      {
      x2[d] = spatialPoint[d] + deltaTime * f1[d] * 0.5;
      }
    x2[InputImageDimension-1] = timePoint + 0.5 * deltaTime * timeDomain;
    VectorType f2 = zeroVector;
    if( this->m_VelocityFieldInterpolator->IsInsideBuffer( x2 ) )
      {
      f2 = this->m_VelocityFieldInterpolator->Evaluate( x2 );
      }

    typename TimeVaryingVelocityFieldType::PointType x3;
    for( unsigned int d = 0; d < OutputImageDimension; d++ )
      {
      x3[d] = spatialPoint[d] + deltaTime * f2[d] * 0.5;
      }
    x3[InputImageDimension-1] = timePoint + 0.5 * deltaTime * timeDomain;
    VectorType f3 = zeroVector;
    if( this->m_VelocityFieldInterpolator->IsInsideBuffer( x3 ) )
      {
      f3 = this->m_VelocityFieldInterpolator->Evaluate( x3 );
      }

    typename TimeVaryingVelocityFieldType::PointType x4;
    for( unsigned int d = 0; d < OutputImageDimension; d++ )
      {
      x4[d] = spatialPoint[d] + deltaTime * f3[d];
      }
    x4[InputImageDimension-1] = timePoint + deltaTime * timeDomain;
    VectorType f4 = zeroVector;
    if( this->m_VelocityFieldInterpolator->IsInsideBuffer( x4 ) )
      {
      f4 = this->m_VelocityFieldInterpolator->Evaluate( x4 );
      }

    spatialPoint += ( deltaTime * ( f1 + f2 * 2.0 + f3 * 2.0 + f4 ) / 6.0 );
    timePoint += deltaTime * timeDomain;
    }
  VectorType displacement = spatialPoint - initialSpatialPoint;

  return displacement;
}

template<class TTimeVaryingVelocityField, class TDeformationField>
void
TimeVaryingVelocityFieldIntegrationImageFilter
  <TTimeVaryingVelocityField, TDeformationField>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "VelocityFieldInterpolator: "
    << this->m_VelocityFieldInterpolator << std::endl;
  os << indent << "LowerTimeBound: " << this->m_LowerTimeBound << std::endl;
  os << indent << "UpperTimeBound: " << this->m_UpperTimeBound << std::endl;
  os << indent << "NumberOfIntegrationSteps: "
    << this->m_NumberOfIntegrationSteps << std::endl;

  if( !this->m_InitialDiffeomorphism.IsNull() )
    {
    os << indent << "InitialDiffeomorphism: " << this->m_InitialDiffeomorphism
      << std::endl;
    os << indent << "DeformationFieldInterpolator: "
      << this->m_DeformationFieldInterpolator << std::endl;
    }
}

}  //end namespace itk

#endif
