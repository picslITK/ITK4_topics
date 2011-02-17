/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsImageToImageMetric.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDemonsImageToImageMetric_txx
#define __itkDemonsImageToImageMetric_txx

#include "itkDemonsImageToImageMetric.h"
#include "itkCovariantVector.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage >
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::DemonsImageToImageMetric()
{
//  this->SetComputeGradient(true);
  this->m_VirtualImage=NULL;
  this->m_ThreaderMSE = NULL;
  this->m_ThreaderDerivatives = NULL;
//  this->m_WithinThreadPreProcess = false;
//  this->m_WithinThreadPostProcess = false;

  m_Normalizer = 1.0;

  //  For backward compatibility, the default behavior is to use all the pixels
  //  in the fixed image.
  //  This should be fixed in ITKv4 so that this metric behaves as the others.
  // this->SetUseAllPixels( true );
}

template < class TFixedImage, class TMovingImage >
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::~DemonsImageToImageMetric()
{
  if(m_ThreaderMSE != NULL)
    {
    delete [] m_ThreaderMSE;
    }
  m_ThreaderMSE = NULL;

  if(m_ThreaderDerivatives != NULL)
    {
    delete [] m_ThreaderDerivatives;
    }
  m_ThreaderDerivatives = NULL;
}

/**
 * Print out internal information about this class
 */
template < class TFixedImage, class TMovingImage  >
void
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{

  Superclass::PrintSelf(os, indent);

}


/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage>
void
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(void) throw ( ExceptionObject )
{
//  this->Superclass::Initialize();
//  this->Superclass::MultiThreadingInitialize();
  // Compute the normalizer
  typename FixedImageType::SpacingType fixedImageSpacing =
    this->m_FixedImage->GetSpacing();

  m_Normalizer = 0.0;
  for ( unsigned int k = 0; k < MovingImageDimension; k++ )
    {
    m_Normalizer += fixedImageSpacing[k] * fixedImageSpacing[k];
    }
  m_Normalizer /= static_cast< double >( MovingImageDimension );

    if ( ! this->m_VirtualImage ) {
    std::cout <<" Allocate virtual image " << std::endl;
      RegionType region;
      region.SetSize(this->GetVirtualDomainSize() );
      region.SetIndex(this->GetVirtualDomainIndex() );
      this->m_VirtualImage = FixedImageType::New();
      this->m_VirtualImage->SetSpacing( this->GetVirtualDomainSpacing() );
      this->m_VirtualImage->SetOrigin( this->GetVirtualDomainOrigin() );
      this->m_VirtualImage->SetDirection( this->GetVirtualDomainDirection() );
      this->m_VirtualImage->SetRegions( region );
      this->m_VirtualImage->Allocate();
      this->m_VirtualImage->FillBuffer( 0 );
      this->m_FixedInterpolator=FixedInterpolatorType::New();
      this->m_MovingInterpolator=MovingInterpolatorType::New();
    }
    if ( this->m_FixedInterpolator ) 
      this->m_FixedInterpolator->SetInputImage(m_FixedImage);
    if ( this->m_MovingInterpolator ) 
      this->m_MovingInterpolator->SetInputImage(m_MovingImage);
    // std::cout <<" allocate-done " << std::endl;

}

/*
template < class TFixedImage, class TMovingImage  >
inline bool
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::GetValueThreadProcessSample( unsigned int threadID,
                               unsigned long fixedImageSample,
                               const MovingImagePointType & itkNotUsed(mappedPoint),
                               double movingImageValue) const
{
  double diff = movingImageValue - this->m_FixedImageSamples[fixedImageSample].value;

  m_ThreaderMSE[threadID] += diff*diff;

  return true;
}

template < class TFixedImage, class TMovingImage  >
typename DemonsImageToImageMetric<TFixedImage,TMovingImage>
::MeasureType
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::GetValue( const ParametersType & parameters ) const
{
  itkDebugMacro("GetValue( " << parameters << " ) ");

  if( !this->m_FixedImage )
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  memset( m_ThreaderMSE,
          0,
          this->m_NumberOfThreads * sizeof(MeasureType) );

  // Set up the parameters in the transform
  this->m_Transform->SetParameters( parameters );
  this->m_Parameters = parameters;

  // MUST BE CALLED TO INITIATE PROCESSING
  this->GetValueMultiThreadedInitiate();

  itkDebugMacro( "Ratio of voxels mapping into moving image buffer: "
                 << this->m_NumberOfPixelsCounted << " / "
                 << this->m_NumberOfFixedImageSamples
                 << std::endl );

  if( this->m_NumberOfPixelsCounted <
      this->m_NumberOfFixedImageSamples / 4 )
    {
    itkExceptionMacro( "Too many samples map outside moving image buffer: "
                       << this->m_NumberOfPixelsCounted << " / "
                       << this->m_NumberOfFixedImageSamples
                       << std::endl );
    }

  double mse = m_ThreaderMSE[0];
  for(unsigned int t=1; t<this->m_NumberOfThreads; t++)
    {
    mse += m_ThreaderMSE[t];
    }
  mse /= this->m_NumberOfPixelsCounted;

  return mse;
}


template < class TFixedImage, class TMovingImage  >
inline bool
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeThreadProcessSample( unsigned int threadID,
                                    unsigned long fixedImageSample,
                                    const MovingImagePointType & itkNotUsed(mappedPoint),
                                    double movingImageValue,
                                    const ImageDerivativesType &
                                    movingImageGradientValue ) const
{
  double diff = movingImageValue - this->m_FixedImageSamples[fixedImageSample].value;

  m_ThreaderMSE[threadID] += diff*diff;

  FixedImagePointType fixedImagePoint = this->m_FixedImageSamples[fixedImageSample].point;

  // Need to use one of the threader transforms if we're
  // not in thread 0.
  //
  // Use a raw pointer here to avoid the overhead of smart pointers.
  // For instance, Register and UnRegister have mutex locks around
  // the reference counts.
  TransformType* transform;

  if (threadID > 0)
    {
    transform = this->m_ThreaderTransform[threadID - 1];
    }
  else
    {
    transform = this->m_Transform;
    }

  // Jacobian should be evaluated at the unmapped (fixed image) point.
  const TransformJacobianType & jacobian = transform
                                               ->GetJacobian( fixedImagePoint );

  double gradientSquaredMagnitude = 0;
  for ( unsigned int j = 0; j < MovingImageDimension; j++ )
    {
    gradientSquaredMagnitude += vnl_math_sqr( movingImageGradientValue[j] );
    }

  const double denominator = vnl_math_sqr( diff ) / m_Normalizer +
    gradientSquaredMagnitude;

  for(unsigned int par=0; par<this->m_NumberOfParameters; par++)
    {
    double sum = 0.0;

    for(unsigned int dim=0; dim<MovingImageDimension; dim++)
      {
      sum += ( diff * movingImageGradientValue[dim] * jacobian( dim, par ) /
        denominator );
      }
    m_ThreaderDerivatives[threadID][par] += sum;
    }

  return true;
}
*/
/**
 * Get the both Value and Derivative Measure
template < class TFixedImage, class TMovingImage  >
void
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivative( const ParametersType & parameters,
                         MeasureType & value,
                         DerivativeType & derivative) const
{

  if( !this->m_FixedImage )
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  // Set up the parameters in the transform
  this->m_Transform->SetParameters( parameters );
  this->m_Parameters = parameters;

  // Reset the joint pdfs to zero
  memset( m_ThreaderMSE,
          0,
          this->m_NumberOfThreads * sizeof(MeasureType) );

  // Set output values to zero
  if(derivative.GetSize() != this->m_NumberOfParameters)
    {
    derivative = DerivativeType( this->m_NumberOfParameters );
    }
  memset( derivative.data_block(),
          0,
          this->m_NumberOfParameters * sizeof(double) );

  for( unsigned int threadID = 0; threadID<this->m_NumberOfThreads; threadID++ )
    {
    memset( m_ThreaderDerivatives[threadID].data_block(),
            0,
            this->m_NumberOfParameters * sizeof(double) );
    }

  // MUST BE CALLED TO INITIATE PROCESSING
  this->GetValueAndDerivativeMultiThreadedInitiate();

  itkDebugMacro( "Ratio of voxels mapping into moving image buffer: "
                 << this->m_NumberOfPixelsCounted << " / "
                 << this->m_NumberOfFixedImageSamples
                 << std::endl );

  if( this->m_NumberOfPixelsCounted <
      this->m_NumberOfFixedImageSamples / 4 )
    {
    itkExceptionMacro( "Too many samples map outside moving image buffer: "
                       << this->m_NumberOfPixelsCounted << " / "
                       << this->m_NumberOfFixedImageSamples
                       << std::endl );
    }

  value = 0;
  for(unsigned int t=0; t<this->m_NumberOfThreads; t++)
    {
    value += m_ThreaderMSE[t];
    for(unsigned int parameter = 0; parameter < this->m_NumberOfParameters;
        parameter++)
      {
      derivative[parameter] += m_ThreaderDerivatives[t][parameter];
      }
    }

  value /= this->m_NumberOfPixelsCounted;
  for(unsigned int parameter = 0; parameter < this->m_NumberOfParameters;
      parameter++)
    {
    derivative[parameter] /= this->m_NumberOfPixelsCounted;
    }
}

 */

/**
 * Get the match measure derivative
template < class TFixedImage, class TMovingImage  >
void
DemonsImageToImageMetric<TFixedImage,TMovingImage>
::GetDerivative( const ParametersType & parameters,
                 DerivativeType & derivative ) const
{
  if( !this->m_FixedImage )
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  MeasureType value;
  // call the combined version
  this->GetValueAndDerivative( parameters, value, derivative );
}
 */

} // end namespace itk


#endif
