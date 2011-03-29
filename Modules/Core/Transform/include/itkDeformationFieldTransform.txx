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
#ifndef __itkDeformationFieldTransform_txx
#define __itkDeformationFieldTransform_txx

#include "itkDeformationFieldTransform.h"

#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_matrix_inverse.h"

namespace itk
{

/**
 * Constructor
 */
template<class TScalar, unsigned int NDimensions>
DeformationFieldTransform<TScalar, NDimensions>::
DeformationFieldTransform() : Superclass( NDimensions, 0 )
{
  this->m_DeformationField = NULL;
  this->m_InverseDeformationField = NULL;

  //Setup and assign default interpolator
  typedef VectorLinearInterpolateImageFunction<DeformationFieldType, ScalarType>
    DefaultInterpolatorType;
  typename DefaultInterpolatorType::Pointer interpolator
    = DefaultInterpolatorType::New();
  this->m_Interpolator = interpolator;

  //Setup and assign parameter helper. This will hold the deformation field
  // for access through the common TransformParameters interface.
  TransformParametersHelperType* helper = new TransformParametersHelperType;
  //After assigning this, m_Parametes will manage this,
  // deleting when appropriate.
  this->m_Parameters.SetHelper( helper );

  m_GaussianSmoothSigma = 3;
  m_DeformationFieldSetTime = 0;
  m_SmoothGaussTempFieldModifiedTime = 0;
  /** These are init'ed when SmoothDeformatFieldGauss is called, either for
   * the first time, or after a new deformation field has been assigned. */
  m_SmoothGaussTempField = NULL;
  m_SmoothGaussSmoother = NULL;

}

/**
 * Destructor
 */
template<class TScalar, unsigned int NDimensions>
DeformationFieldTransform<TScalar, NDimensions>::
~DeformationFieldTransform()
{
}

/**
 * Transform point
 */
template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::OutputPointType
DeformationFieldTransform<TScalar, NDimensions>
::TransformPoint( const InputPointType& inputPoint ) const
{
  if( !this->m_DeformationField )
    {
    itkExceptionMacro( "No deformation field is specified." );
    }
  if( !this->m_Interpolator )
    {
    itkExceptionMacro( "No interpolator is specified." );
    }

  typename InterpolatorType::ContinuousIndexType cidx;
  typename InterpolatorType::PointType point;
  point.CastFrom( inputPoint );

  OutputPointType outputPoint;
  outputPoint.CastFrom( inputPoint );

  if( this->m_Interpolator->IsInsideBuffer( point ) )
    {
    this->m_DeformationField->
      TransformPhysicalPointToContinuousIndex( point, cidx );
    typename InterpolatorType::OutputType displacement =
      this->m_Interpolator->EvaluateAtContinuousIndex( cidx );
    outputPoint += displacement;
    }
  else
    {
    ScalarType infinity;
    if( vcl_numeric_limits<ScalarType>::has_infinity )
      {
      infinity = vcl_numeric_limits<ScalarType>::infinity();
      }
    else
      {
      infinity = NumericTraits<ScalarType>::max();
      }
    outputPoint.Fill( infinity );
    }

  return outputPoint;
}

/**
 * Transform covariant vector
 */

template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::OutputCovariantVectorType
DeformationFieldTransform<TScalar, NDimensions>
::TransformCovariantVector( const InputCovariantVectorType& vector,
                            const InputPointType & point) const
{
  if( !this->m_DeformationField )
    {
    itkExceptionMacro( "No deformation field is specified." );
    }
  if( !this->m_Interpolator )
    {
    itkExceptionMacro( "No interpolator is specified." );
    }

  JacobianType jacobian;
  /* Get the inverse Jacobian directly for efficiency. It means we don't have
   * to compute an SVD inverse. */
  this->GetInverseJacobianOfForwardFieldWithRespectToPosition( point, jacobian );

  OutputCovariantVectorType result;

  for ( unsigned int i = 0; i < NDimensions; i++ )
    {
    result[i] = NumericTraits< ScalarType >::Zero;
    for ( unsigned int j = 0; j < NDimensions; j++ )
      {
      result[i] += jacobian[j][i] * vector[j]; // Inverse
                                               // transposed
      }
    }

  return result;
}

template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::OutputVectorPixelType
DeformationFieldTransform<TScalar, NDimensions>
::TransformCovariantVector( const InputVectorPixelType& vector,
                            const InputPointType & point) const
{
  if( !this->m_DeformationField )
    {
    itkExceptionMacro( "No deformation field is specified." );
    }
  if( !this->m_Interpolator )
    {
    itkExceptionMacro( "No interpolator is specified." );
    }

  JacobianType jacobian;
  /* Get the inverse Jacobian directly for efficiency. It means we don't have
   * to compute an SVD inverse. */
  this->GetInverseJacobianOfForwardFieldWithRespectToPosition( point, jacobian );

  const unsigned int numberOfComponents = NumericTraits< InputVectorPixelType >::GetLength( vector );

  OutputVectorPixelType result;     // Converted vector
  result.SetSize( numberOfComponents );

  JacobianType dataJacobian;
  dataJacobian.SetSize( numberOfComponents, numberOfComponents );

  for ( unsigned int i = 0; i < numberOfComponents; i++ )
    {
    if ( i < NDimensions )
      {
      result[i] = NumericTraits< ScalarType >::Zero;
      for ( unsigned int j = 0; j < NDimensions; j++ )
        {
        result[i] += jacobian[j][i] * vector[j];
        }
      }
    else
      {
      result[i] = vector[i];
      }
    }

  return result;
}


/**
 * Transform vector
 */
template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::OutputVectorType
DeformationFieldTransform<TScalar, NDimensions>
::TransformVector( const InputVectorType& vector, const InputPointType & point ) const
{
  if( !this->m_DeformationField )
    {
    itkExceptionMacro( "No deformation field is specified." );
    }
  if( !this->m_Interpolator )
    {
    itkExceptionMacro( "No interpolator is specified." );
    }

  JacobianType jacobian;
  this->GetJacobianWithRespectToPosition( point, jacobian );

  AffineTransformPointer localTransform = AffineTransformType::New();
  localTransform->SetIdentity();
  localTransform->SetMatrix( jacobian );
  return localTransform->TransformVector( vector );
}


/**
 * Transform vector
 */
template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::OutputVnlVectorType
DeformationFieldTransform<TScalar, NDimensions>
::TransformVector( const InputVnlVectorType& vector, const InputPointType & point ) const
{
  if( !this->m_DeformationField )
    {
    itkExceptionMacro( "No deformation field is specified." );
    }
  if( !this->m_Interpolator )
    {
    itkExceptionMacro( "No interpolator is specified." );
    }

  JacobianType jacobian;
  this->GetJacobianWithRespectToPosition( point, jacobian );

  AffineTransformPointer localTransform = AffineTransformType::New();
  localTransform->SetIdentity();
  localTransform->SetMatrix( jacobian );
  return localTransform->TransformVector( vector );
}

/**
 * Transform vector
 */
template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::OutputVectorPixelType
DeformationFieldTransform<TScalar, NDimensions>
::TransformVector( const InputVectorPixelType& vector, const InputPointType & point ) const
{
  if( !this->m_DeformationField )
    {
    itkExceptionMacro( "No deformation field is specified." );
    }
  if( !this->m_Interpolator )
    {
    itkExceptionMacro( "No interpolator is specified." );
    }

  JacobianType jacobian;
  this->GetJacobianWithRespectToPosition( point, jacobian );

  AffineTransformPointer localTransform = AffineTransformType::New();
  localTransform->SetIdentity();
  localTransform->SetMatrix( jacobian );
  return localTransform->TransformVector( vector );
}

/**
 * Transform tensor
 */
template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::OutputTensorType
DeformationFieldTransform<TScalar, NDimensions>
::TransformTensor( const InputTensorType& inputTensor, const InputPointType & point ) const
{
  if( !this->m_DeformationField )
    {
    itkExceptionMacro( "No deformation field is specified." );
    }
  if( !this->m_Interpolator )
    {
    itkExceptionMacro( "No interpolator is specified." );
    }

  //JacobianType jacobian;
  //this->GetJacobianWithRespectToPosition( point, jacobian );

  //Get Tensor-space version of local transform (i.e. always 3D)
  typedef MatrixOffsetTransformBase<ScalarType, InputTensorType::Dimension, InputTensorType::Dimension> EigenVectorTransformType;
  typename  EigenVectorTransformType::MatrixType matrix;
  typename  EigenVectorTransformType::MatrixType dMatrix;
  matrix.Fill(0.0);
  dMatrix.Fill(0.0);
  for (unsigned int i=0; i<InputTensorType::Dimension; i++)
    {
    matrix(i,i) = 1.0;
    dMatrix(i,i) = 1.0;
    }

  JacobianType invJacobian;
  this->GetInverseJacobianOfForwardFieldWithRespectToPosition( point, invJacobian );

  for (unsigned int i=0; i<NDimensions; i++)
    {
    for (unsigned int j=0; j<NDimensions; j++)
      {
      if ( (i < InputTensorType::Dimension) && (j < InputTensorType::Dimension))
        {
        matrix(i,j) = invJacobian(i,j);
        dMatrix(i,j) = this->GetDirectionChangeMatrix()(i,j);
        }
      }
    }

  typename InputTensorType::EigenValuesArrayType eigenValues;
  typename InputTensorType::EigenVectorsMatrixType eigenVectors;
  inputTensor.ComputeEigenAnalysis( eigenValues, eigenVectors );

  InputTensorEigenVectorType ev1;
  InputTensorEigenVectorType ev2;
  InputTensorEigenVectorType ev3;

  for (unsigned int i=0; i<InputTensorType::Dimension; i++)
    {
    ev1[i] = eigenVectors(2,i);
    ev2[i] = eigenVectors(1,i);
    }

  // Account for image direction changes between moving and fixed spaces
  ev1 = matrix * dMatrix * ev1;
  ev1.Normalize();

  // Get aspect of rotated e2 that is perpendicular to rotated e1
  ev2 = matrix * dMatrix * ev2;
  double dp = ev2 * ev1;
  if ( dp < 0 )
    {
    ev2 = ev2*(-1.0);
    dp = dp*(-1.0);
    }
  ev2 = ev2 - dp*ev1;
  ev2.Normalize();

  itk::CrossHelper<InputTensorEigenVectorType> vectorCross;
  ev3 = vectorCross( ev1, ev2 );

  // Outer product matrices
  typename EigenVectorTransformType::MatrixType e1;
  typename EigenVectorTransformType::MatrixType e2;
  typename EigenVectorTransformType::MatrixType e3;
  for (unsigned int i=0; i<InputTensorType::Dimension; i++)
    {
    for (unsigned int j=0; j<InputTensorType::Dimension; j++)
      {
      e1(i,j) = eigenValues[2] * ev1[i]*ev1[j];
      e2(i,j) = eigenValues[1] * ev2[i]*ev2[j];
      e3(i,j) = eigenValues[0] * ev3[i]*ev3[j];
      }
    }

  typename EigenVectorTransformType::MatrixType rotated = e1 + e2 + e3;

  OutputTensorType result;     // Converted vector
  result[0] = rotated(0,0);
  result[1] = rotated(0,1);
  result[2] = rotated(0,2);
  result[3] = rotated(1,1);
  result[4] = rotated(1,2);
  result[5] = rotated(2,2);

  return result;

  //AffineTransformPointer localTransform = AffineTransformType::New();
  //localTransform->SetIdentity();
  //localTransform->SetMatrix( jacobian );
  //localTransform->SetDirectionChangeMatrix( this->GetDirectionChangeMatrix() );
  //return localTransform->TransformTensor( inputTensor );

}

template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::OutputVectorPixelType
DeformationFieldTransform<TScalar, NDimensions>
::TransformTensor( const InputVectorPixelType& inputTensor, const InputPointType & point ) const
{
  if( !this->m_DeformationField )
    {
    itkExceptionMacro( "No deformation field is specified." );
    }
  if( !this->m_Interpolator )
    {
    itkExceptionMacro( "No interpolator is specified." );
    }

  OutputVectorPixelType result( InputTensorType::InternalDimension );     // Converted tensor
  result.Fill( 0.0 );

  InputTensorType dt(0.0);
  const unsigned int tDim = inputTensor.Size();
  for (unsigned int i=0; i<tDim; i++)
    {
    dt[i] = inputTensor[i];
    }

  OutputTensorType outDT = this->TransformTensor( dt, point );

  for (unsigned int i=0; i<InputTensorType::InternalDimension; i++)
    {
    result[i] = outDT[i];
    }

  return result;


  /*
  JacobianType jacobian;
  this->GetJacobianWithRespectToPosition( point, jacobian );

  AffineTransformPointer localTransform = AffineTransformType::New();
  localTransform->SetIdentity();
  localTransform->SetMatrix( jacobian );
  return localTransform->TransformTensor( inputTensor );
  */
}


/**
 * return an inverse transformation
 */
template<class TScalar, unsigned int NDimensions>
bool DeformationFieldTransform<TScalar, NDimensions>
::GetInverse( Self *inverse ) const
{
  if ( !inverse || !this->m_InverseDeformationField )
    {
    return false;
    }
  else
    {
    inverse->SetDeformationField( this->m_InverseDeformationField );
    inverse->SetInverseDeformationField( this->m_DeformationField );
    inverse->SetInterpolator( this->m_Interpolator );

    return true;
    }
}

// Return an inverse of this transform
template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::InverseTransformBasePointer
DeformationFieldTransform<TScalar, NDimensions>
::GetInverseTransform() const
{
  Pointer inverseTransform = New();
  if( this->GetInverse( inverseTransform ) )
    {
    return inverseTransform.GetPointer();
    }
  else
    {
    return NULL;
    }
}

/*
 * GetJacobian methods
 */

template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::JacobianType &
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobian( const InputPointType & point ) const
{
  itkExceptionMacro( "GetJacobian() not valid for DeformationFieldTransform. Use GetJacobianWithRespectToPosition()" );
}

template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobianWithRespectToPosition( const InputPointType & point,
                                      JacobianType & jacobian )
                                                                          const
{
  IndexType idx;
  this->m_DeformationField->TransformPhysicalPointToIndex( point, idx );
  this->GetJacobianWithRespectToPosition( idx, jacobian );
}

template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobianWithRespectToPosition( const IndexType & index,
                                      JacobianType & jacobian )
                                                                          const
{
  this->GetJacobianWithRespectToPositionInternal( index, jacobian, false );
}

template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetInverseJacobianOfForwardFieldWithRespectToPosition(
                                      const InputPointType & point,
                                      JacobianType & jacobian,
                                      bool useSVD )
                                                                          const
{
  IndexType idx;
  this->m_DeformationField->TransformPhysicalPointToIndex( point, idx );
  this->GetInverseJacobianOfForwardFieldWithRespectToPosition( idx, jacobian,
                                                               useSVD );
}

template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetInverseJacobianOfForwardFieldWithRespectToPosition(
                                      const IndexType & index,
                                      JacobianType & jacobian,
                                      bool useSVD )
                                                                          const
{
  if (useSVD)
    {
    this->GetJacobianWithRespectToPositionInternal( index, jacobian, false );
    AffineTransformPointer localTransform = AffineTransformType::New();
    localTransform->SetIdentity();
    localTransform->SetMatrix( jacobian );
    typename TransformType::Pointer invJacobian
      = localTransform->GetInverseTransform();

    for (unsigned int i=0; i<jacobian.rows(); i++)
      for (unsigned int j=0; j<jacobian.cols(); j++)
        jacobian(i,j) = invJacobian->GetParameters()[ j + i*jacobian.rows() ];

    }
  else
    {
    this->GetJacobianWithRespectToPositionInternal( index, jacobian, true );
    }
}

/*
 * GetJacobianWithRespectToPositionInternal. Worker method.
 */
template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobianWithRespectToPositionInternal( const IndexType & index,
                                      JacobianType & jacobian,
                                      bool doInverseJacobian )
                                                                          const
{
  jacobian.SetSize(NDimensions,NDimensions);
  //This may not be necessary. Double-check below.
  // jacobian.Fill(0.0);

  typename DeformationFieldType::SizeType size =
                this->m_DeformationField->GetLargestPossibleRegion().GetSize();
  typename DeformationFieldType::SpacingType spacing =
                                        this->m_DeformationField->GetSpacing();

  IndexType ddrindex;
  IndexType ddlindex;
  IndexType difIndex[NDimensions][2];

  unsigned int posoff=1;
  float difspace=1.0;
  float space=1.0;
  if (posoff == 0) difspace=1.0;
  float mindist=1.0;
  float dist=100.0;
  bool oktosample=true;
  float dPixSign = doInverseJacobian ? -1.0 : 1.0;

  for (unsigned int row=0; row<NDimensions; row++)
    {
    dist = fabs((float)index[row]);
    if (dist < mindist)
      {
      oktosample = false;
      }
    dist = fabs((float)size[row] - (float)index[row]);
    if (dist < mindist)
      {
      oktosample = false;
      }
    }

  if ( oktosample )
    {
    OutputVectorType cpix = this->m_DeformationField->GetPixel(index);
    m_DeformationField->TransformLocalVectorToPhysicalVector( cpix, cpix );
    //cpix = directionRaw->TransformVector( cpix );

    // itkCentralDifferenceImageFunction does not support vector images
    // so do this manually here
    for(unsigned int row=0; row< NDimensions;row++)
      {
      difIndex[row][0]=index;
      difIndex[row][1]=index;
      ddrindex=index;
      ddlindex=index;
      if ((int) index[row] < (int)(size[row]-2) )
        {
        difIndex[row][0][row] = index[row]+posoff;
        ddrindex[row] = index[row]+posoff*2;
        }
      if (index[row] > 1 )
        {
        difIndex[row][1][row] = index[row]-1;
        ddlindex[row] = index[row]-2;
        }

      float h=1;
      space=1.0; // should use image spacing here?

      OutputVectorType rpix = m_DeformationField->GetPixel( difIndex[row][1] );
      OutputVectorType lpix = m_DeformationField->GetPixel( difIndex[row][0] );
      OutputVectorType rrpix = m_DeformationField->GetPixel( ddrindex );
      OutputVectorType llpix = m_DeformationField->GetPixel( ddlindex );

      m_DeformationField->TransformLocalVectorToPhysicalVector( rpix, rpix );
      m_DeformationField->TransformLocalVectorToPhysicalVector( rrpix, rrpix );
      m_DeformationField->TransformLocalVectorToPhysicalVector( lpix, lpix );
      m_DeformationField->TransformLocalVectorToPhysicalVector( llpix, llpix );

      rpix =  rpix*h  + cpix*(1.-h);
      lpix =  lpix*h  + cpix*(1.-h);
      rrpix = rrpix*h + rpix*(1.-h);
      llpix = llpix*h + lpix*(1.-h);

      //4th order centered difference
      OutputVectorType dPix =
          ( lpix*8.0 + llpix - rrpix - rpix*8.0 ) * space / (12.0) * dPixSign;

      //typename DeformationFieldType::PixelType dPix=
      //      ( lpix - rpix )*space/(2.0*h); //2nd order centered difference

      for(unsigned int col=0; col< NDimensions; col++)
        {
        float val = dPix[col] / spacing[col];
        if (row == col)
          {
          val += 1.0;
          }
        jacobian(col,row) = val;
        }
      } // for row
    } //if oktosample

  for (unsigned int jx = 0; jx < NDimensions; jx++)
    {
    for (unsigned int jy = 0; jy < NDimensions; jy++)
      {
      if ( !vnl_math_isfinite(jacobian(jx,jy))  )
        {
        oktosample = false;
        }
      }
    }

  if ( !oktosample )
    {
    jacobian.Fill(0.0);
    for (unsigned int i=0; i<NDimensions; i++)
      {
      jacobian(i,i) = 1.0;
      }
    }
}

template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::UpdateTransformParameters( DerivativeType & update, ScalarType factor)
{
  //This simply adds the values.
  //TODO: This should be multi-threaded probably, via image add filter.
  Superclass::UpdateTransformParameters( update, factor );

  //Now we smooth the result. Not thread safe. Does it's own
  // threading.
  SmoothDeformationFieldGauss();
}

template<class TScalar, unsigned int NDimensions>
void DeformationFieldTransform<TScalar, NDimensions>
::SmoothDeformationFieldGauss()
{
  itkDebugMacro(" enter gauss smooth. m_GaussianSmoothSigma: "
                << m_GaussianSmoothSigma);
  if( this->m_GaussianSmoothSigma <= 0 )
    {
    return;
    }

  typename DeformationFieldType::Pointer field = this->GetDeformationField();

  /* Allocate temp field if new deformation field has been set.
   * We only want to allocate this field if this method is used */
  if( this->GetDeformationFieldSetTime() >
      this->m_SmoothGaussTempFieldModifiedTime )
    {
    this->m_SmoothGaussTempFieldModifiedTime = this->GetMTime();
    m_SmoothGaussTempField = DeformationFieldType::New();
    m_SmoothGaussTempField->SetSpacing( field->GetSpacing() );
    m_SmoothGaussTempField->SetOrigin( field->GetOrigin() );
    m_SmoothGaussTempField->SetDirection( field->GetDirection() );
    m_SmoothGaussTempField->SetLargestPossibleRegion(
                                          field->GetLargestPossibleRegion() );
    m_SmoothGaussTempField->SetRequestedRegion( field->GetRequestedRegion() );
    m_SmoothGaussTempField->SetBufferedRegion( field->GetBufferedRegion() );
    m_SmoothGaussTempField->Allocate();

    //This should only be allocated once as well, for efficiency.
    m_SmoothGaussSmoother = SmoothGaussSmootherType::New();
    }

  if( m_SmoothGaussTempField.IsNull() )
    {
    itkExceptionMacro("Expected m_SmoothGaussTempField to be allocated.");
    }

  typedef typename DeformationFieldType::PixelType    VectorType;
  typedef typename VectorType::ValueType              ScalarType;

  typedef typename DeformationFieldType::PixelContainerPointer
                                                        PixelContainerPointer;
  // I think we need to keep this as SmartPointer type, to preserve the
  // reference counting so we can assign the swapPtr to the main field and
  // not have to do a memory copy - this happens when image dimensions are odd.
  PixelContainerPointer swapPtr;

  // graft the output field onto the mini-pipeline
  m_SmoothGaussSmoother->GraftOutput( m_SmoothGaussTempField );

  for( unsigned int j = 0; j < Dimension; j++ )
    {
    // smooth along this dimension
    m_SmoothGaussOperator.SetDirection( j );
    m_SmoothGaussOperator.SetVariance( m_GaussianSmoothSigma );
    m_SmoothGaussOperator.SetMaximumError(0.001 );
    m_SmoothGaussOperator.SetMaximumKernelWidth( 256 );
    m_SmoothGaussOperator.CreateDirectional();

    // todo: make sure we only smooth within the buffered region
    m_SmoothGaussSmoother->SetOperator( m_SmoothGaussOperator );
    m_SmoothGaussSmoother->SetInput( field );
    try
      {
      m_SmoothGaussSmoother->Update();
      }
    catch( ExceptionObject & exc )
      {
      std::string msg("Caught exception: ");
      msg += exc.what();
      itkExceptionMacro( << msg );
      }

    if( j < Dimension - 1 )
      {
      // swap the containers
      swapPtr = m_SmoothGaussSmoother->GetOutput()->GetPixelContainer();
      m_SmoothGaussSmoother->GraftOutput( field );
      // SetPixelContainer does a smartpointer assignment, so the pixel
      // container won't be deleted if field  points to the
      // temporary field upon exiting this method.
      field->SetPixelContainer( swapPtr );
      m_SmoothGaussSmoother->Modified();
      }
    }

  if( Dimension % 2 == 0 )
    {
    // For even number of dimensions, the final pass writes the output
    // into field's original pixel container, so we just point back to that.
    // And point the temporary field back to its original container for next
    // time through.
    m_SmoothGaussTempField->SetPixelContainer( field->GetPixelContainer() );
    field->SetPixelContainer(
                    m_SmoothGaussSmoother->GetOutput()->GetPixelContainer() );
    }

  //make sure boundary does not move
  ScalarType weight = 1.0;
  if (m_GaussianSmoothSigma < 0.5)
    {
    weight=1.0 - 1.0 * ( this->m_GaussianSmoothSigma / 0.5);
    }
  ScalarType weight2 = 1.0 - weight;
  typedef ImageRegionIteratorWithIndex<DeformationFieldType> Iterator;
  typename DeformationFieldType::SizeType size =
                                field->GetLargestPossibleRegion().GetSize();
  Iterator outIter( field, field->GetLargestPossibleRegion() );
  for( outIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter )
  {
    bool onboundary=false;
    typename DeformationFieldType::IndexType index= outIter.GetIndex();
    for (int i=0; i < Dimension; i++)
      {
      if (index[i] < 1 || index[i] >= static_cast<int>( size[i] )-1 )
        {
        onboundary=true;
        }
      }
    if( onboundary )
      {
      VectorType vec;
      vec.Fill(0.0);
      outIter.Set(vec);
      }
    else
      {
      VectorType svec = m_SmoothGaussSmoother->GetOutput()->GetPixel( index );
      outIter.Set( svec * weight + outIter.Get() * weight2);
      }
  }

  itkDebugMacro("done gauss smooth ");
}

template<class TScalar, unsigned int NDimensions>
void DeformationFieldTransform<TScalar, NDimensions>
::SetDeformationField( DeformationFieldType* field )
{
  itkDebugMacro("setting DeformationField to " << field);
  if ( this->m_DeformationField != field )
    {
    this->m_DeformationField = field;
    this->Modified();
    /* Store this separately for use in smoothing because we only want
     * to know when the deformation field object has changed, not just
     * its contents. */
    this->m_DeformationFieldSetTime = this->GetMTime();
    if( ! this->m_Interpolator.IsNull() )
      {
      this->m_Interpolator->SetInputImage( this->m_DeformationField );
      }
    //Assign to parameters object
    this->m_Parameters.SetParametersObject( this->m_DeformationField );
    }
}

template<class TScalar, unsigned int NDimensions>
void DeformationFieldTransform<TScalar, NDimensions>
::SetInterpolator( InterpolatorType* interpolator )
{
  itkDebugMacro("setting Interpolator to " << interpolator);
  if ( this->m_Interpolator != interpolator )
    {
    this->m_Interpolator = interpolator;
    this->Modified();
    if( ! this->m_DeformationField.IsNull() )
      {
      this->m_Interpolator->SetInputImage( this->m_DeformationField );
      }
    }
}

template <class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os,indent );

  std::cout << indent << "Interpolator: " << std::endl;
  std::cout << indent << indent << this->m_Interpolator << std::endl;

  if( this->m_DeformationField )
    {
    std::cout << indent << "Deformation Field: " << std::endl;
    std::cout << indent << indent << this->m_DeformationField << std::endl;
    }

  if( this->m_InverseDeformationField )
    {
    std::cout << indent << "Inverse Deformation Field: " << std::endl;
    std::cout << indent << indent << this->m_InverseDeformationField << std::endl;
    }
}
} // namespace itk

#endif
