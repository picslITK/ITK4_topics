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
#include "vnl/algo/vnl_symmetric_eigensystem.h"

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
::TransformCovariantVector( const InputCovariantVectorType& vector, const InputPointType & point ) const
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
  OutputCovariantVectorType outputVector = localTransform->TransformCovariantVector( vector );

  return outputVector;

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
  OutputVectorType outputVector = localTransform->TransformVector( vector );

  return outputVector;

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

  JacobianType jacobian;
  this->GetJacobianWithRespectToPosition( point, jacobian );

  LocalTransformPointer localTransform = LocalTransformType::New();
  localTransform->SetIdentity();
  localTransform->SetMatrix( jacobian );

  OutputTensorType result = localTransform->TransformTensor( inputTensor );

  return result;
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


template<class TScalar, unsigned int NDimensions>
typename DeformationFieldTransform<TScalar, NDimensions>::JacobianType &
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobian( const InputPointType & point ) const
{
  itkExceptionMacro( "GetJacobian() not valid for DeformationFieldTransform. Use GetJacobianWithRespectToPosition()" );
}


/*
template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobianWithRespectToParameters( const InputPointType & point,
                                      JacobianType & jacobian ) const
{
  jacobian.SetSize(NDimensions,NDimensions);
  jacobian.Fill(0.0);
  for (unsigned int i=0; i<NDimensions; i++)
    {
    jacobian(i,i) = 1.0;
    }
}

template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobianWithRespectToParameters( const IndexType & index ,
                                      JacobianType & jacobian ) const
{
  jacobian.SetSize(NDimensions,NDimensions);
  jacobian.Fill(0.0);
  for (unsigned int i=0; i<NDimensions; i++)
    {
    jacobian(i,i) = 1.0;
    }
}
*/

template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobianWithRespectToPosition( const InputPointType & point,
                                      JacobianType & jacobian ) const
{
  IndexType idx;
  this->m_DeformationField->TransformPhysicalPointToIndex( point, idx );
  this->GetJacobianWithRespectToPosition( idx, jacobian );
}

template<class TScalar, unsigned int NDimensions>
void
DeformationFieldTransform<TScalar, NDimensions>
::GetJacobianWithRespectToPosition( const IndexType & index,
                                      JacobianType & jacobian ) const
{
  jacobian.SetSize(NDimensions,NDimensions);
  jacobian.Fill(0.0);

  AffineTransformPointer direction = AffineTransformType::New();
  direction->SetIdentity();
  direction->SetMatrix( this->m_DeformationField->GetDirection() );
  typename DeformationFieldType::SizeType size = this->m_DeformationField->GetLargestPossibleRegion().GetSize();
  typename DeformationFieldType::SpacingType spacing = this->m_DeformationField->GetSpacing();

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
    cpix = direction->TransformVector( cpix );

  // itkCentralDifferenceImageFunction does not support vector images so do this manually here
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

    //if (this->m_UseImageDirection)
    //{
    rpix = direction->TransformVector( rpix );
    lpix = direction->TransformVector( lpix );
    rrpix = direction->TransformVector( rrpix );
    llpix = direction->TransformVector( llpix );
    //}

    rpix = rpix*h+cpix*(1.-h);
    lpix = lpix*h+cpix*(1.-h);
    rrpix = rrpix*h+rpix*(1.-h);
    llpix = llpix*h+lpix*(1.-h);

    OutputVectorType dPix = ( lpix*8.0 + llpix - rrpix - rpix*8.0 )*space/(12.0); //4th order centered difference

    //typename DeformationFieldType::PixelType dPix=( lpix - rpix )*space/(2.0*h); //2nd order centered difference

    for(unsigned int col=0; col< NDimensions; col++)
      {
      float val = dPix[col] / spacing[col];
      if (row == col)
        {
        val += 1.0;
        }
      jacobian(col,row) = val;
      }
    }
  }

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
void DeformationFieldTransform<TScalar, NDimensions>
::SetDeformationField( DeformationFieldType* field )
{
  itkDebugMacro("setting DeformationField to " << field);
  if ( this->m_DeformationField != field )
    {
    this->m_DeformationField = field;
    this->Modified();
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
