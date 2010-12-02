/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDeformationFieldTransform.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDeformationFieldTransform_txx
#define __itkDeformationFieldTransform_txx

#include "itkDeformationFieldTransform.h"

#include "itkVectorLinearInterpolateImageFunction.h"

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
  this->m_PreviousDeformationFieldMTime = 0;
  this->m_PreviousInterpolatorMTime = 0;

  typedef VectorLinearInterpolateImageFunction<DeformationFieldType, ScalarType>
    DefaultInterpolatorType;
  typename DefaultInterpolatorType::Pointer interpolator
    = DefaultInterpolatorType::New();
  this->m_Interpolator = interpolator;
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
  /* Check if either the deformation field or iterpolatr have changed since
   * we were last in here. */
  if( this->m_DeformationField->GetMTime() >
        this->m_PreviousDeformationFieldMTime ||
      this->m_Interpolator->GetMTime() > this->m_PreviousInterpolatorMTime )
    {
    this->m_Interpolator->SetInputImage( this->m_DeformationField );
    }
  this->m_PreviousDeformationFieldMTime = this->GetMTime();
  this->m_PreviousInterpolatorMTime = this->GetMTime();

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
  /** We implement the Jacobian " by hand " to allow improved efficiency ( not well
    * achieved in this implementation ) and because we need the implementation to
    * use physical space finite differences.
    */
  /** FIXME --- here we compute the Jacobian only at a point because that is all it influences.
   *  However, we should store this Jacobian value at the right place in the full list of
   *  Jacobians or otherwise (somehow) avoid storing the full parameter Jacobian list at all
   *  and just return this local Jacobian. If we do this, then we are reducing the full, dense
   *  deformation field down to the local neighborhood and view it only from there.
   *  This strategy is good b/c it reduces storage but is bad b/c it requires computing the
   *  Jacobian each time it will be used.  Should we also implement the Hessian?
   */
  typedef double Scalar;
  enum { ImageDimension = DeformationFieldType::ImageDimension };
  typename DeformationFieldType::PointType difIndex[ImageDimension][2];

  /** FIXME explicitly resizing jacobian is necessary for now but should be done by base class */
  this->m_Jacobian.SetSize(ImageDimension,ImageDimension);
  this->m_Jacobian.Fill( 0.0 );
  for(unsigned int col=0; col< ImageDimension;col++) this->m_Jacobian[col][col]=1;

  /**FIXME should store these variables as member variables */
  typename DeformationFieldType::SizeType fieldSize=this->m_DeformationField->GetLargestPossibleRegion().GetSize();
  typename DeformationFieldType::SpacingType fieldSpacing=this->m_DeformationField->GetSpacing();

  /**FIXME below , need standard bounds checking , none performed here , should implement s.t. oktosample is replaced with bounds check at  point and point +/- 1 spacing value */
  bool oktosample=true;
  if (oktosample)
  {
    for(unsigned int row=0; row< ImageDimension;row++)
    {
      difIndex[row][0]=point;
      difIndex[row][1]=point;
      difIndex[row][0][row]=point[row]+fieldSpacing[row];
      difIndex[row][1][row]=point[row]-fieldSpacing[row];
      typename DeformationFieldType::PixelType rpix = const_cast< Self * >( this )->ReorientVectorByDirection(difIndex[row][1]);
      typename DeformationFieldType::PixelType lpix = const_cast< Self * >( this )->ReorientVectorByDirection(difIndex[row][0]);
      for(unsigned int col=0; col< ImageDimension;col++){
      /** here is the physical space central difference calculation */
      /** Note: this same Jacobian matrix is used to reorient vectors or tensors */
        Scalar val=(lpix[col]-rpix[col]);
        if (row == col) val=val/(2.0*fieldSpacing[col])+1.0;
        else val = val/(2.0*fieldSpacing[col]);
        this->m_Jacobian[col][row]=val;
      }
    }
    /** the determinant of the jacobian matrix, should probably calculate explicitly */
    Scalar det = vnl_determinant(this->m_Jacobian);
    if (det < 0.0)
    {
      /** FIXME throw warning, return identity matrix and/or remove det calculation
        * entirely b/c of computational cost
        */
      this->m_Jacobian.Fill( 0.0 );
      for(unsigned int col=0; col< ImageDimension;col++) this->m_Jacobian[col][col]=1;
    }
  }//oktosample if
  return this->m_Jacobian;
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
