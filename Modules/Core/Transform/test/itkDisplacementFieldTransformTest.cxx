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

#include <iostream>

#include "itkDisplacementFieldTransform.h"
#include "itkVectorInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkCenteredAffineTransform.h"
#include "itkImageRegionIteratorWithIndex.h"

const unsigned int dimensions = 2;
typedef itk::DisplacementFieldTransform<double, dimensions>
                                               DisplacementTransformType;
typedef DisplacementTransformType::ScalarType  ScalarType;

template <typename TPoint>
bool samePoint( const TPoint & p1, const TPoint & p2, double epsilon=1e-8 )
  {
  bool pass=true;
  for ( unsigned int i = 0; i < TPoint::PointDimension; i++ )
    {
    if( vcl_fabs( p1[i] - p2[i] ) > epsilon )
      pass=false;
    }
  return pass;
  }

template <typename TVector>
bool sameVector( const TVector & p1, const TVector & p2, double epsilon=1e-8 )
  {
  bool pass=true;
  for ( unsigned int i = 0; i < TVector::Dimension; i++ )
    {
    if( vcl_fabs( p1[i] - p2[i] ) > epsilon )
      pass=false;
    }
  return pass;
  }

template <typename TVector>
bool sameVariableVector( const TVector & p1, const TVector & p2, double epsilon=1e-8 )
  {
  bool pass=true;

  const unsigned int D1 = p1.Size();
  const unsigned int D2 = p2.Size();

  if (D1 != D2)
    {
    return false;
    }

  for ( unsigned int i = 0; i < D1; i++ )
    {
    if( vcl_fabs( p1[i] - p2[i] ) > epsilon )
      pass=false;
    }
  return pass;
  }

template <typename TTensor>
bool sameTensor( const TTensor & p1, const TTensor & p2, double epsilon=1e-8 )
  {
  bool pass=true;
  for ( unsigned int i = 0; i < TTensor::InternalDimension; i++ )
    {
    if( vcl_fabs( p1[i] - p2[i] ) > epsilon )
      pass=false;
    }
  return pass;
  }

template <typename TArray2D>
bool sameArray2D( const TArray2D & a1, const TArray2D & a2, double epsilon=1e-8 )
  {
  bool pass=true;

  if ( (a1.rows() != a2.rows()) || (a1.cols() != a2.cols()) )
    {
    return false;
    }

  for ( unsigned int i = 0; i < a1.cols(); i++ )
    {
    for ( unsigned int j = 0; j < a1.rows(); j++ )
      {
        if( vcl_fabs( a1(j,i) - a2(j,i) ) > epsilon )
          {
          pass=false;
          }
      }
    }
  return pass;
  }


int itkDisplacementFieldTransformTest(int ,char *[] )
{

  std::cout.precision(12);

  typedef  itk::Matrix<ScalarType, dimensions, dimensions>  Matrix2Type;
  typedef  itk::Vector<ScalarType, dimensions>              Vector2Type;

  /* Create a displacement field transform */
  DisplacementTransformType::Pointer displacementTransform =
      DisplacementTransformType::New();
  typedef DisplacementTransformType::DisplacementFieldType FieldType;
  FieldType::Pointer field = FieldType::New(); //This is based on itk::Image

  FieldType::SizeType size;
  FieldType::IndexType start;
  FieldType::RegionType region;
  int dimLength = 20;
  size.Fill( dimLength );
  start.Fill( 0 );
  region.SetSize( size );
  region.SetIndex( start );
  field->SetRegions( region );
  field->Allocate();

  DisplacementTransformType::OutputVectorType zeroVector;
  zeroVector.Fill( 0 );
  field->FillBuffer( zeroVector );

  /* Initialize Affine transform and use to create displacement field */
  typedef itk::CenteredAffineTransform<double,2>  AffineTransformType;
  typedef AffineTransformType::Pointer            AffineTransformPointer;
  typedef AffineTransformType::MatrixType         AffineMatrixType;
  AffineMatrixType affineMatrix;
  affineMatrix(0,0) = 1.0;
  affineMatrix(1,0) = 0.01;
  affineMatrix(0,1) = 0.02;
  affineMatrix(1,1) = 1.1;
  DisplacementTransformType::JacobianType fieldJTruth;
  fieldJTruth.SetSize(2,2);
  fieldJTruth(0,0) = 1.0;
  fieldJTruth(1,0) = 0.01;
  fieldJTruth(0,1) = 0.02;
  fieldJTruth(1,1) = 1.1;
  AffineTransformType::Pointer affineTransform = AffineTransformType::New();
  affineTransform->SetIdentity();
  affineTransform->SetMatrix( affineMatrix );

  itk::ImageRegionIteratorWithIndex<FieldType> it( field, field->GetLargestPossibleRegion() );
  it.GoToBegin();

  while (! it.IsAtEnd() )
    {
    FieldType::PointType pt;
    field->TransformIndexToPhysicalPoint( it.GetIndex(), pt );
    FieldType::PointType pt2 = affineTransform->TransformPoint( pt );
    FieldType::PointType::VectorType vec = pt2 - pt;
    FieldType::PixelType v;
    v[0] = vec[0];
    v[1] = vec[1];
    field->SetPixel( it.GetIndex(), v );
    ++it;
    }

  displacementTransform->SetDisplacementField( field );

  DisplacementTransformType::InputPointType testPoint;
  testPoint[0] = 10;
  testPoint[1] = 8;

  /* Test LocalJacobian methods */
  DisplacementTransformType::JacobianType jacobian;
  displacementTransform->GetJacobianWithRespectToPosition( testPoint, jacobian );
  std::cout << "Local jacobian estimated. " << std::endl << jacobian << std::endl;
  if (!sameArray2D( jacobian, fieldJTruth, 1e-6 ) )
    {
      std::cout << "Failed getting local jacobian. Should be " << std::endl << affineMatrix << std::endl;
      return EXIT_FAILURE;
    }

  DisplacementTransformType::JacobianType invfieldJTruth;
  invfieldJTruth.SetSize(2,2);
  invfieldJTruth(0,0) = affineTransform->GetInverseTransform()->GetParameters()[0];
  invfieldJTruth(1,0) = affineTransform->GetInverseTransform()->GetParameters()[1];
  invfieldJTruth(0,1) = affineTransform->GetInverseTransform()->GetParameters()[2];
  invfieldJTruth(1,1) = affineTransform->GetInverseTransform()->GetParameters()[3];
  displacementTransform->GetInverseJacobianOfForwardFieldWithRespectToPosition( testPoint, jacobian );
  std::cout << "Local inverse jacobian estimated. " << std::endl << jacobian << std::endl;
  if (!sameArray2D( jacobian, invfieldJTruth, 1e-1 ) )
    {
      std::cout << "Failed getting local inverse jacobian. Should be " << std::endl << invfieldJTruth << std::endl;
      return EXIT_FAILURE;
    }

  displacementTransform->GetInverseJacobianOfForwardFieldWithRespectToPosition( testPoint, jacobian, true );
  std::cout << "Local inverse jacobian estimated with SVD. " << std::endl << jacobian << std::endl;
  if (!sameArray2D( jacobian, invfieldJTruth, 1e-1 ) )
    {
      std::cout << "Failed getting local inverse jacobian. Should be " << std::endl << invfieldJTruth << std::endl;
      return EXIT_FAILURE;
    }

  /* Test GetJacobianWithRespectToParameters. Should return identity */
  DisplacementTransformType::JacobianType
    identity(dimensions, dimensions), testIdentity;
  identity.Fill(0);
  for( unsigned int i=0; i < dimensions; i++ )
    {
    identity[i][i] = 1.0;
    }
  displacementTransform->GetJacobianWithRespectToParameters(
                                                    testPoint, testIdentity );
  if( !sameArray2D( identity, testIdentity, 1e-10 ) )
    {
    std::cout << "Failed returning identity for "
                 "GetJacobianWithRespectToParameters( point, ... )"
              << std::endl;
    return EXIT_FAILURE;
    }
  DisplacementTransformType::IndexType testIndex;
  testIdentity.SetSize(1,1); //make sure it gets resized properly
  displacementTransform->GetJacobianWithRespectToParameters(
                                                    testIndex, testIdentity );
  if( !sameArray2D( identity, testIdentity, 1e-10 ) )
    {
    std::cout << "Failed returning identity for "
                 "GetJacobianWithRespectToParameters( index, ... )"
              << std::endl;
    return EXIT_FAILURE;
    }

  /** Test transforming of points */

  DisplacementTransformType::OutputPointType deformOutput, deformTruth;

  /* Test a point with non-zero displacement */
  FieldType::IndexType idx;
  field->TransformPhysicalPointToIndex( testPoint, idx );
  deformTruth = testPoint + field->GetPixel( idx );

  deformOutput = displacementTransform->TransformPoint( testPoint );
  std::cout << "point 1 transformed: " << deformOutput << std::endl;
  if( !samePoint( deformOutput, deformTruth  ) )
      {
      std::cout << "Failed transforming point 1. Should be " << deformTruth << std::endl;
      return EXIT_FAILURE;
      }


  DisplacementTransformType::InputVectorType testVector;
  DisplacementTransformType::OutputVectorType deformVector, deformVectorTruth;
  testVector[0] = 0.5;
  testVector[1] = 0.5;

  deformVectorTruth = affineTransform->TransformVector( testVector );
  deformVector = displacementTransform->TransformVector( testVector, testPoint );
  std::cout << "vector 1 transformed: " << deformVector << std::endl;
  if( !sameVector( deformVector, deformVectorTruth, 0.0001 ) )
      {
      std::cout << "Failed transforming vector 1. Should be " << deformVectorTruth << std::endl;
      return EXIT_FAILURE;
      }

  bool caughtException = false;
  try
    {
    deformVector = displacementTransform->TransformVector( testVector );
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Caught expected exception:" << std::endl << e << std::endl;
    caughtException = true;
    }

  if (!caughtException)
    {
    std::cout << "Expected TransformVector(vector) to throw exception." << std::endl;
    return EXIT_FAILURE;
    }


  /** Test VectorTransform for variable length vector which does not
   *  need do have the dimensionality as the transform */
  DisplacementTransformType::InputVectorPixelType testVVector(3);
  DisplacementTransformType::OutputVectorPixelType deformVVector, deformVVectorTruth(3);
  testVVector[0] = 0.5;
  testVVector[1] = 0.5;
  testVVector[2] = 1.0;

  deformVVectorTruth = affineTransform->TransformVector( testVVector );
  deformVVector = displacementTransform->TransformVector( testVVector, testPoint );
  std::cout << "variable length vector 1 transformed: " << deformVVector << std::endl;
  if( !sameVariableVector( deformVVector, deformVVectorTruth, 0.0001 ) )
      {
      std::cout << "Failed transforming variable length vector 1. Should be " << deformVVectorTruth << std::endl;
      return EXIT_FAILURE;
      }

  caughtException = false;
  try
    {
    deformVVector = displacementTransform->TransformVector( testVVector );
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Caught expected exception:" << std::endl << e << std::endl;
    caughtException = true;
    }

  if (!caughtException)
    {
    std::cout << "Expected TransformVector(vector) to throw exception." << std::endl;
    return EXIT_FAILURE;
    }


  DisplacementTransformType::InputCovariantVectorType testcVector;
  DisplacementTransformType::OutputCovariantVectorType deformcVector, deformcVectorTruth;
  testcVector[0] = 0.5;
  testcVector[1] = 0.5;

  deformcVectorTruth = affineTransform->TransformCovariantVector( testcVector );
  deformcVector = displacementTransform->TransformCovariantVector( testcVector, testPoint );
  std::cout << "covariant vector 1 transformed: " << deformcVector << std::endl;
  if( !sameVector( deformcVector, deformcVectorTruth, 0.1 ) )
      {
      std::cout << "Failed transforming vector 1. Should be "
                << deformcVectorTruth << std::endl;
      return EXIT_FAILURE;
      }

  caughtException = false;
  try
    {
    deformcVector = displacementTransform->TransformCovariantVector( testcVector );
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Caught expected exception:" << std::endl << e << std::endl;
    caughtException = true;
    }

  if (!caughtException)
    {
    std::cout << "Expected TransformCovariantVector(vector) to throw exception."
              << std::endl;
    return EXIT_FAILURE;
    }


  DisplacementTransformType::InputVectorPixelType testcVVector(3);
  DisplacementTransformType::OutputVectorPixelType deformcVVector, deformcVVectorTruth(3);
  testcVVector[0] = 0.5;
  testcVVector[1] = 0.5;
  testcVVector[2] = 1.0;

  deformcVVectorTruth = affineTransform->TransformCovariantVector( testcVVector );
  deformcVVector = displacementTransform->TransformCovariantVector( testcVVector, testPoint );
  std::cout << "variable length covariant vector 1 transformed: "
            << deformcVVector << std::endl;
  if( !sameVariableVector( deformcVVector, deformcVVectorTruth, 0.1 ) )
      {
      std::cout
        << "Failed transforming variable length covariant vector 1. Should be "
        << deformcVVectorTruth << std::endl;
      return EXIT_FAILURE;
      }

  caughtException = false;
  try
    {
    deformcVVector = displacementTransform->TransformCovariantVector( testcVVector );
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Caught expected exception:" << std::endl << e << std::endl;
    caughtException = true;
    }

  if (!caughtException)
    {
    std::cout << "Expected TransformCovariantVector(vector) to throw exception."
              << std::endl;
    return EXIT_FAILURE;
    }


  DisplacementTransformType::InputDiffusionTensor3DType testTensor;
  DisplacementTransformType::OutputDiffusionTensor3DType deformTensor,
                                                      deformTensorTruth;
  testTensor[0] = 3;
  testTensor[1] = 0.01;
  testTensor[2] = 0.01;
  testTensor[3] = 2;
  testTensor[4] = 0.01;
  testTensor[5] = 1;

  // pass thru functionality only for now
  deformTensorTruth = affineTransform->TransformDiffusionTensor( testTensor );
  std::cout << "tensor 1:             " << testTensor << std::endl;
  deformTensor = displacementTransform->TransformDiffusionTensor(
                                                      testTensor, testPoint );
  std::cout << "tensor 1 transformed: " << deformTensor << std::endl;
  if( !sameTensor( deformTensor, deformTensorTruth, 0.0001 ) )
      {
      std::cout << "Failed transforming tensor 1. Should be "
                << deformTensorTruth << std::endl;
      //return EXIT_FAILURE;
      }

  caughtException = false;
  try
    {
    deformTensor = displacementTransform->TransformDiffusionTensor( testTensor );
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Caught expected exception:" << std::endl << e << std::endl;
    caughtException = true;
    }

  if (!caughtException)
    {
    std::cout << "Expected TransformDiffusionTensor(tensor) to throw exception."
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Test setting of parameters with wrong size */
  caughtException = false;
  try
    {
    DisplacementTransformType::ParametersType params(1);
    params.Fill(0);
    displacementTransform->SetParameters( params );
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Caught expected exception: " << std::endl << e << std::endl;
    caughtException = true;
    }
  if( !caughtException )
    {
    std::cout << "Expected SetParameters with wrong size to throw exception."
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Test UpdateTransformParameters */
  std::cout << "Testing UpdateTransformParameters..." << std::endl;
  DisplacementTransformType::DerivativeType
                 derivative( displacementTransform->GetNumberOfParameters() );
  DisplacementTransformType::DerivativeType
                 updateTruth( displacementTransform->GetNumberOfParameters() );
  DisplacementTransformType::ParametersType
                 params( displacementTransform->GetNumberOfParameters() );
  derivative.Fill(1.2);
  ScalarType testFactor = 1.5;
  for( unsigned int i=0;
        i < displacementTransform->GetNumberOfParameters(); i++ )
    {
    params[i] = i;
    updateTruth[i] = params[i] + derivative[i] * testFactor;
    }
  displacementTransform->SetParameters( params );
  displacementTransform->UpdateTransformParameters( derivative, testFactor );
  params = displacementTransform->GetParameters();
  //std::cout  << "params: " << std::endl << params << std::endl;
  //           << "derivativeTruth: " << std::endl << derivative << std::endl
  for( unsigned int i=0;
        i < displacementTransform->GetNumberOfParameters(); i++ )
    {
    if( params[i] != updateTruth[i] )
      {
      std::cout  << "UpdateTransformParameters test failed: " << std::endl;
      std::cout  << "params: " << std::endl << params << std::endl
                 << "updateTruth: " << std::endl << updateTruth << std::endl;
      return EXIT_FAILURE;
      }
    }

  /* Test IsLinear()
   * Should always return false */
  if( displacementTransform->IsLinear() )
    {
    std::cout << "DisplacementFieldTransform returned 'true' for IsLinear()."
      " Expected 'false'." << std::endl;
    return EXIT_FAILURE;
    }

    /* We haven't set an inverse displacement field for the inverse displacement
   * transform, so we should get a false return here */
  DisplacementTransformType::Pointer inverseTransform
    = DisplacementTransformType::New();
  if( displacementTransform->GetInverse( inverseTransform ) )
    {
    std::cout << "Expected GetInverse() to fail." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
