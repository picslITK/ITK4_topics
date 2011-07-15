/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkDeformationFieldTransformTest.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>

#include "itkDeformationFieldTransform.h"
#include "itkVectorInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkCenteredAffineTransform.h"
#include "itkImageRegionIteratorWithIndex.h"

const unsigned int dimensions = 2;
typedef itk::DeformationFieldTransform<double, dimensions>
                                              DeformationTransformType;
typedef DeformationTransformType::ScalarType  ScalarType;

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


int itkDeformationFieldTransformTest(int ,char *[] )
{

  /* NOTE
   * Requires updating when DeformationFieldTransform is
   * fully implemented with operations for transforming
   * vectors as well. Currently this just tests transforming
   * points. */

  std::cout.precision(12);

  typedef  itk::Matrix<ScalarType, dimensions, dimensions>  Matrix2Type;
  typedef  itk::Vector<ScalarType, dimensions>               Vector2Type;

  /* Create a deformation field transform */
  DeformationTransformType::Pointer deformationTransform =
      DeformationTransformType::New();
  typedef DeformationTransformType::DeformationFieldType FieldType;
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

  DeformationTransformType::OutputVectorType zeroVector;
  zeroVector.Fill( 0 );
  field->FillBuffer( zeroVector );

  /* Initialize Affine transform and use to create deformation field */
  typedef itk::CenteredAffineTransform<double,2> AffineTransformType;
  typedef AffineTransformType::Pointer AffineTransformPointer;
  typedef AffineTransformType::MatrixType AffineMatrixType;
  AffineMatrixType affineMatrix;
  affineMatrix(0,0) = 1.0;
  affineMatrix(1,0) = 0.01;
  affineMatrix(0,1) = 0.02;
  affineMatrix(1,1) = 1.1;
  DeformationTransformType::JacobianType fieldJTruth;
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

  deformationTransform->SetDeformationField( field );

  DeformationTransformType::InputPointType testPoint;
  testPoint[0] = 10;
  testPoint[1] = 8;

  /* Test LocalJacobian methods */
  DeformationTransformType::JacobianType jacobian;
  deformationTransform->GetJacobianWithRespectToPosition( testPoint, jacobian );
  std::cout << "Local jacobian estimated. " << std::endl << jacobian << std::endl;
  if (!sameArray2D( jacobian, fieldJTruth, 1e-6 ) )
    {
      std::cout << "Failed getting local jacobian. Should be " << std::endl << affineMatrix << std::endl;
      return EXIT_FAILURE;
    }

  DeformationTransformType::JacobianType invfieldJTruth;
  invfieldJTruth.SetSize(2,2);
  invfieldJTruth(0,0) = affineTransform->GetInverseTransform()->GetParameters()[0];
  invfieldJTruth(1,0) = affineTransform->GetInverseTransform()->GetParameters()[1];
  invfieldJTruth(0,1) = affineTransform->GetInverseTransform()->GetParameters()[2];
  invfieldJTruth(1,1) = affineTransform->GetInverseTransform()->GetParameters()[3];
  deformationTransform->GetInverseJacobianOfForwardFieldWithRespectToPosition( testPoint, jacobian );
  std::cout << "Local inverse jacobian estimated. " << std::endl << jacobian << std::endl;
  if (!sameArray2D( jacobian, invfieldJTruth, 1e-1 ) )
    {
      std::cout << "Failed getting local inverse jacobian. Should be " << std::endl << invfieldJTruth << std::endl;
      return EXIT_FAILURE;
    }

  deformationTransform->GetInverseJacobianOfForwardFieldWithRespectToPosition( testPoint, jacobian, true );
  std::cout << "Local inverse jacobian estimated with SVD. " << std::endl << jacobian << std::endl;
  if (!sameArray2D( jacobian, invfieldJTruth, 1e-1 ) )
    {
      std::cout << "Failed getting local inverse jacobian. Should be " << std::endl << invfieldJTruth << std::endl;
      return EXIT_FAILURE;
    }

  /* Test GetJacobianWithRespectToParameters. Should return identity */
  DeformationTransformType::JacobianType
    identity(dimensions, dimensions), testIdentity;
  identity.Fill(0);
  for( unsigned int i=0; i < dimensions; i++ )
    {
    identity[i][i] = 1.0;
    }
  deformationTransform->GetJacobianWithRespectToParameters(
                                                    testPoint, testIdentity );
  if( !sameArray2D( identity, testIdentity, 1e-10 ) )
    {
    std::cout << "Failed returning identity for "
                 "GetJacobianWithRespectToParameters( point, ... )"
              << std::endl;
    return EXIT_FAILURE;
    }
  DeformationTransformType::IndexType testIndex;
  testIdentity.SetSize(1,1); //make sure it gets resized properly
  deformationTransform->GetJacobianWithRespectToParameters(
                                                    testIndex, testIdentity );
  if( !sameArray2D( identity, testIdentity, 1e-10 ) )
    {
    std::cout << "Failed returning identity for "
                 "GetJacobianWithRespectToParameters( index, ... )"
              << std::endl;
    return EXIT_FAILURE;
    }

  /** Test transforming of points */

  DeformationTransformType::OutputPointType deformOutput, deformTruth;

  /* Test a point with non-zero deformation */
  FieldType::IndexType idx;
  field->TransformPhysicalPointToIndex( testPoint, idx );
  deformTruth = testPoint + field->GetPixel( idx );

  deformOutput = deformationTransform->TransformPoint( testPoint );
  std::cout << "point 1 transformed: " << deformOutput << std::endl;
  if( !samePoint( deformOutput, deformTruth  ) )
      {
      std::cout << "Failed transforming point 1. Should be " << deformTruth << std::endl;
      return EXIT_FAILURE;
      }


  DeformationTransformType::InputVectorType testVector;
  DeformationTransformType::OutputVectorType deformVector, deformVectorTruth;
  testVector[0] = 0.5;
  testVector[1] = 0.5;

  deformVectorTruth = affineTransform->TransformVector( testVector );
  deformVector = deformationTransform->TransformVector( testVector, testPoint );
  std::cout << "vector 1 transformed: " << deformVector << std::endl;
  if( !sameVector( deformVector, deformVectorTruth, 0.0001 ) )
      {
      std::cout << "Failed transforming vector 1. Should be " << deformVectorTruth << std::endl;
      return EXIT_FAILURE;
      }

  bool caughtException = false;
  try
    {
    deformVector = deformationTransform->TransformVector( testVector );
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
  DeformationTransformType::InputVectorPixelType testVVector(3);
  DeformationTransformType::OutputVectorPixelType deformVVector, deformVVectorTruth(3);
  testVVector[0] = 0.5;
  testVVector[1] = 0.5;
  testVVector[2] = 1.0;

  deformVVectorTruth = affineTransform->TransformVector( testVVector );
  deformVVector = deformationTransform->TransformVector( testVVector, testPoint );
  std::cout << "variable length vector 1 transformed: " << deformVVector << std::endl;
  if( !sameVariableVector( deformVVector, deformVVectorTruth, 0.0001 ) )
      {
      std::cout << "Failed transforming variable length vector 1. Should be " << deformVVectorTruth << std::endl;
      return EXIT_FAILURE;
      }

  caughtException = false;
  try
    {
    deformVVector = deformationTransform->TransformVector( testVVector );
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


  DeformationTransformType::InputCovariantVectorType testcVector;
  DeformationTransformType::OutputCovariantVectorType deformcVector, deformcVectorTruth;
  testcVector[0] = 0.5;
  testcVector[1] = 0.5;

  deformcVectorTruth = affineTransform->TransformCovariantVector( testcVector );
  deformcVector = deformationTransform->TransformCovariantVector( testcVector, testPoint );
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
    deformcVector = deformationTransform->TransformCovariantVector( testcVector );
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


  DeformationTransformType::InputVectorPixelType testcVVector(3);
  DeformationTransformType::OutputVectorPixelType deformcVVector, deformcVVectorTruth(3);
  testcVVector[0] = 0.5;
  testcVVector[1] = 0.5;
  testcVVector[2] = 1.0;

  deformcVVectorTruth = affineTransform->TransformCovariantVector( testcVVector );
  deformcVVector = deformationTransform->TransformCovariantVector( testcVVector, testPoint );
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
    deformcVVector = deformationTransform->TransformCovariantVector( testcVVector );
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


  DeformationTransformType::InputTensorType testTensor;
  DeformationTransformType::OutputTensorType deformTensor, deformTensorTruth;
  testTensor[0] = 3;
  testTensor[1] = 0.01;
  testTensor[2] = 0.01;
  testTensor[3] = 2;
  testTensor[4] = 0.01;
  testTensor[5] = 1;

  // pass thru functionality only for now
  deformTensorTruth = affineTransform->TransformTensor( testTensor );
  std::cout << "tensor 1:             " << testTensor << std::endl;
  deformTensor = deformationTransform->TransformTensor( testTensor, testPoint );
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
    deformTensor = deformationTransform->TransformTensor( testTensor );
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Caught expected exception:" << std::endl << e << std::endl;
    caughtException = true;
    }

  if (!caughtException)
    {
    std::cout << "Expected TransformTensor(tensor) to throw exception."
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Test setting of parameters with wrong size */
  caughtException = false;
  try
    {
    DeformationTransformType::ParametersType params(1);
    params.Fill(0);
    deformationTransform->SetParameters( params );
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

  /* Test SmoothDeformationFieldGauss */
  std::cout << "Test SmoothDeformationFieldGauss" << std::endl;
  DeformationTransformType::ParametersType params;
  DeformationTransformType::ParametersType
                  paramsFill( deformationTransform->GetNumberOfParameters() );
  DeformationTransformType::ParametersValueType paramsFillValue = 1.0;
  paramsFill.Fill( paramsFillValue );
  // Add an outlier to see that some smoothing is taking place.
  unsigned int outlier = dimLength*dimensions*4 + dimLength*dimensions / 2;
  paramsFill( outlier ) = 9999.0;
  paramsFill( outlier + 1 ) = 9999.0;
  deformationTransform->SetGaussianSmoothSigma(3);
  deformationTransform->SetParameters( paramsFill );
  params = deformationTransform->GetParameters();
  //std::cout << "field->GetPixelContainter *before* Smooth: "
  //          << field->GetPixelContainer() << std::endl;
  //std::cout << "params *before* SmoothDeformationFieldGauss: " << std::endl
  //          << params << std::endl;
  deformationTransform->SmoothDeformationFieldGauss();
  params = deformationTransform->GetParameters();
  //std::cout << "field->GetPixelContainter *after* Smooth: "
  //          << field->GetPixelContainer() << std::endl;
  /* print out result. We should see 0's on all boundaries from the smoothing
   * routine, and then some smoothing around the outlier we set above. */
  std::cout << "Parameters *after* SmoothDeformationFieldGauss, around "
            << "outlier: " << std::endl;
  for(int i=-2; i< 3; i++ )
    {
     for(int j=-2; j< 3; j++ )
      {
      unsigned int index = outlier +
        (unsigned int) (i * (signed int)(dimLength*dimensions) + j);
      std::cout << params(index) << " ";
      if( params(index) == paramsFillValue )
        {
        std::cout << "Expected to read a smoothed value at this index."
                  << " Instead, read " << params(index) << std::endl;
        return EXIT_FAILURE;
        }
      }
    std::cout << std::endl;
    }


  /* Test UpdateTransformParameters */
  std::cout << "Testing UpdateTransformParameters..." << std::endl;
  /* fill with 0 */
  field->FillBuffer( zeroVector );
  DeformationTransformType::DerivativeType
                    derivative( deformationTransform->GetNumberOfParameters() );
  derivative.Fill(1.2);
  deformationTransform->UpdateTransformParameters( derivative );
  params = deformationTransform->GetParameters();
  //TODO create derivateTruth with zero's on boundaries and compare to
  // results from UpdateTransformParameters
  //std::cout  << "params: " << std::endl << params << std::endl;
             //<< "derivativeTruth: " << std::endl << derivative << std::endl

  /* Update with an uneven field to verify some smoothing is happening. */
  field->FillBuffer( zeroVector );
  derivative.Fill( 0 );
  derivative( dimLength*2 + dimLength / 2 ) = 123456789.0;
  deformationTransform->UpdateTransformParameters( derivative );
  params = deformationTransform->GetParameters();
  std::cout << "UpdateTransformParameters with uneven update: " << std::endl
            << "params: " << std::endl << params << std::endl;

  /* Test IsLinear()
   * Should always return false */
  if( deformationTransform->IsLinear() )
    {
    std::cout << "DeformationFieldTransform returned 'true' for IsLinear()."
      " Expected 'false'." << std::endl;
    return EXIT_FAILURE;
    }

    /* We haven't set an inverse deformation field for the inverse deformation
   * transform, so we should get a false return here */
  DeformationTransformType::Pointer inverseTransform
    = DeformationTransformType::New();
  if( deformationTransform->GetInverse( inverseTransform ) )
    {
    std::cout << "Expected GetInverse() to fail." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
