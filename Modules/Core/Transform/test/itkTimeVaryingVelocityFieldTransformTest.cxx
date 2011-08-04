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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImage.h"
#include "itkTimeVaryingVelocityFieldIntegrationImageFilter.h"
#include "itkTimeVaryingVelocityFieldTransform.h"
#include "itkVector.h"

int itkTimeVaryingVelocityFieldTransformTest( int, char* [] )
{
  typedef itk::Vector<double, 3>    VectorType;
  typedef itk::Image<VectorType, 3> DeformationFieldType;
  typedef itk::Image<VectorType, 4> TimeVaryingVelocityFieldType;

  TimeVaryingVelocityFieldType::PointType origin;
  origin.Fill( 0.0 );
  TimeVaryingVelocityFieldType::SpacingType spacing;
  spacing.Fill( 2.0 );
  TimeVaryingVelocityFieldType::SizeType size;
  size.Fill( 25 );
  VectorType velocity;
  velocity.Fill( 0.1 );

  TimeVaryingVelocityFieldType::Pointer timeVaryingVelocityField =
    TimeVaryingVelocityFieldType::New();
  timeVaryingVelocityField->SetOrigin( origin );
  timeVaryingVelocityField->SetSpacing( spacing );
  timeVaryingVelocityField->SetRegions( size );
  timeVaryingVelocityField->Allocate();
  timeVaryingVelocityField->FillBuffer( velocity );


  DeformationFieldType::PointType df_origin;
  df_origin.Fill( 0.0 );
  DeformationFieldType::SpacingType df_spacing;
  df_spacing.Fill( 2.0 );
  DeformationFieldType::SizeType df_size;
  df_size.Fill( 25 );
  VectorType zeroVector;
  zeroVector.Fill( 0.0 );

  DeformationFieldType::Pointer initialDiffeomorphism =
    DeformationFieldType::New();
  initialDiffeomorphism->SetOrigin( df_origin );
  initialDiffeomorphism->SetSpacing( df_spacing );
  initialDiffeomorphism->SetRegions( df_size );
  initialDiffeomorphism->Allocate();
  initialDiffeomorphism->FillBuffer( zeroVector );

  typedef itk::TimeVaryingVelocityFieldIntegrationImageFilter
    <TimeVaryingVelocityFieldType, DeformationFieldType> IntegratorType;

  IntegratorType::Pointer integrator = IntegratorType::New();
  integrator->SetInput( timeVaryingVelocityField );
  integrator->SetInitialDiffeomorphism( initialDiffeomorphism );
  integrator->SetLowerTimeBound( 0.3 );
  integrator->SetUpperTimeBound( 0.75 );
  integrator->SetNumberOfIntegrationSteps( 10 );
  integrator->Update();

  integrator->Print( std::cout, 3 );

  DeformationFieldType::IndexType index;
  index.Fill( 0 );
  VectorType displacement;

  // This integration should result in a constant image of value
  // 0.75 * 0.1 - 0.3 * 0.1 = 0.045 with ~epsilon deviation
  // due to numerical computations
  displacement = integrator->GetOutput()->GetPixel( index );
  if( vnl_math_abs( displacement[0] - 0.045 ) > 0.0001 )
    {
    std::cerr << "Failed to produce the correct forward integration."
      << std::endl;
    std::cerr << "  instead of 0.045, we got " << displacement[0] << std::endl;
    return EXIT_FAILURE;
    }

  IntegratorType::Pointer inverseIntegrator = IntegratorType::New();
  inverseIntegrator->SetInput( timeVaryingVelocityField );
  inverseIntegrator->SetLowerTimeBound( 1.0 );
  inverseIntegrator->SetUpperTimeBound( 0.0 );
  inverseIntegrator->SetNumberOfIntegrationSteps( 10 );
  inverseIntegrator->Update();

  // This integration should result in a constant image of value
  // -( 0.1 * 1.0 - ( 0.1 * 0.0 ) ) = -0.1 with ~epsilon deviation
  // due to numerical computations
  displacement = inverseIntegrator->GetOutput()->GetPixel( index );
  if( vnl_math_abs( displacement[0] + 0.1 ) > 0.0001 )
    {
    std::cerr << "Failed to produce the correct inverse integration."
      << std::endl;
    std::cerr << "  instead of -0.1, we got " << displacement[0] << std::endl;
    return EXIT_FAILURE;
    }

  // Now test the transform

  typedef itk::TimeVaryingVelocityFieldTransform<double, 3> TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetLowerTimeBound( 0.0 );
  transform->SetUpperTimeBound( 1.0 );
  transform->SetTimeVaryingVelocityField( timeVaryingVelocityField );

  TransformType::InputPointType point;
  point.Fill( 1.3 );

  TransformType::OutputPointType transformedPoint =
    transform->TransformPoint( point );

  point += velocity;
  if( point.EuclideanDistanceTo( transformedPoint ) > 0.001 )
    {
    std::cerr << "Failed to produce the expected transformed point."
      << std::endl;
    return EXIT_FAILURE;
    }
  point -= velocity;

  TransformType::InputPointType point2;
  point2.CastFrom( transformedPoint );
  transformedPoint = transform->GetInverseTransform()->TransformPoint( point2 );

  if( point.EuclideanDistanceTo( transformedPoint ) > 0.001 )
    {
    std::cerr << "Failed to produce the expected inverse transformed point."
      << std::endl;
    return EXIT_FAILURE;
    }

  transform->Print( std::cout, 3 );

  return EXIT_SUCCESS;
}
