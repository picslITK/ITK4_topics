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
  this->m_InputImageVectorLength=0;
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
//      this->m_VirtualImage->SetVectorLength(1);
      this->m_VirtualImage->Allocate();
      this->m_VirtualImage->FillBuffer( 0 );
      this->m_FixedInterpolator=FixedInterpolatorType::New();
      this->m_MovingInterpolator=MovingInterpolatorType::New();
    }
    if ( this->m_FixedInterpolator )
      this->m_FixedInterpolator->SetInputImage(this->m_FixedImage);
    if ( this->m_MovingInterpolator )
      this->m_MovingInterpolator->SetInputImage(this->m_MovingImage);
    // std::cout <<" allocate-done " << std::endl;
//  if ( this->m_FixedImage->GetVectorLength() != this->m_MovingImage->GetVectorLength() )
//    itkExceptionMacro( << "Fixed image vector length does not equal moving image vector length." );
//  this->m_InputImageVectorLength=this->m_FixedImage->GetVectorLength();
}




} // end namespace itk


#endif
