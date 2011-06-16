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

/** This function computes the local voxel-wise contribution of
 *  the metric to the global integral of the metric/derivative.
 */
/* NOTE - why aren't all params passed by reference ? */
double ComputeLocalContributionToMetricAndDerivative(PointType mappedFixedPoint, PointType mappedMovingPoint, ImageDerivativesType movingImageGradientValue, TransformJacobianType jacobian , DerivativeType& localDerivative )
{
  double metricval=0;
  /** Only the voxelwise contribution given the point pairs. */

  FixedImagePixelType fpix=this->m_FixedInterpolator->Evaluate(mappedFixedPoint);
  MovingImagePixelType mpix=this->m_MovingInterpolator->Evaluate(mappedMovingPoint);
  FixedImagePixelType diff = fpix - mpix ;
  metricval+=fabs(diff)/(double)FixedImageDimension;
  /* Stauff: this is ok. It will actually just return an identity matrix. */
  this->m_MovingImageTransform->GetJacobianWithRespectToParameters(
    mappedMovingPoint, jacobian);

  for ( unsigned int par = 0; par < this->m_MovingImageTransform->GetNumberOfLocalParameters(); par++ )
  {
    double sum = 0.0;
    for (unsigned int c=0; c < this->m_InputImageVectorLength; c++)
    {
      for ( unsigned int dim = 0; dim < MovingImageDimension; dim++ )
      {
  //    sum += 2.0 *diff[c]*jacobian(dim, par);// * movingImageGradientValue[dim];// VectorLength
       sum += 2.0 *diff*jacobian(dim, par)*movingImageGradientValue[dim];// VectorLength
      }
    }
    localDerivative[par]+=sum;
  }
  return metricval;
}

/** This function is calls the ComputeMetricAndDerivative() function
 *  over the domain of interest.
 */
double ComputeMetricAndDerivative(const ImageRegionType &thread_region, DerivativeType& derivative )
{
  /** This should be set if the image vector length is nonzero! */
  this->m_InputImageVectorLength=1;
  DerivativeType localDerivative(this->m_MovingImageTransform->GetNumberOfLocalParameters());
  std::cout <<"Demons: alloced local derivative of size " << localDerivative.Size()<< std::endl;
  typedef typename MovingImageType::OffsetValueType OffsetValueType;
  TransformJacobianType jacobian(FixedImageDimension,this->m_MovingImageTransform->GetNumberOfLocalParameters());
  jacobian.Fill(0);
  std::cout <<"Demons: alloced jacobian " << std::endl;
  /** TODO
   *  1. Define derivative type and how to access its entries
   *  2. How do we compute image gradients in both moving & fixed space for both
   *  random and dense derivative estimates?
   */

  /** For each location in the virtual domain, map to both the fixed and moving space
   *  and compute the values of the voxels in the corresponding locations.  There should
   *  be a transform between the virtual space and the fixed/moving space s.t. the images
   *  are interpolated in an unbiased manner.
   */
  double metric_sum=0;
  unsigned long ct=0;
  ImageRegionConstIteratorWithIndex<FixedImageType> ItV( this->m_VirtualImage,
      thread_region );

 /////Stauffer 6/6
 //Shouldn't need to warp to the virtual image if we're not going to calc image deriv from it.
//  std::cout << " warp the image " << std::endl;
  /* compute the image gradient */
//  ItV.GoToBegin();
//  while( !ItV.IsAtEnd() )
//  {
    /** use the fixed and moving transforms to compute the corresponding points.*/
/*    bool sampleOk = true;
    // convert the index to a point
    PointType mappedPoint;
    PointType mappedMovingPoint;
    this->m_VirtualImage->TransformIndexToPhysicalPoint(ItV.GetIndex(),mappedPoint);
    mappedMovingPoint = this->m_MovingImageTransform->TransformPoint(mappedPoint);
    if (  !this->m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint) )
         sampleOk=false;
    if ( sampleOk )
      {
        MovingImagePixelType mpix=this->m_MovingInterpolator->Evaluate(mappedMovingPoint);
        this->m_VirtualImage->SetPixel(ItV.GetIndex(),mpix);
      }
    ++ItV;
  }
  std::cout << " warp the image done " << std::endl;
 */
///////////

  typename DerivativeFunctionType::Pointer derivativeCalculator = DerivativeFunctionType::New();
  derivativeCalculator->UseImageDirectionOn();

  ////// Stauffer 6/6
  // Calc deriv from moving image (or the "in use" transform(s), I think, when we get that added)
  //  and transform result into virtual space to avoid threading issues that arise from boundary
  //  conditions when using a warp to virtual image
  // derivativeCalculator->SetInputImage(this->m_VirtualImage);
  derivativeCalculator->SetInputImage(this->m_MovingImageTransform);
  //////

  //    std::cout << " thread region  " << thread_region << std::endl;
  ItV.GoToBegin();
  while( !ItV.IsAtEnd() )
  {
    /** use the fixed and moving transforms to compute the corresponding points.*/
    bool sampleOk = true;
    // convert the index to a point
    PointType mappedPoint;
    PointType mappedFixedPoint;
    PointType mappedMovingPoint;
    this->m_VirtualImage->TransformIndexToPhysicalPoint(ItV.GetIndex(),mappedPoint);
    mappedFixedPoint = this->m_FixedImageTransform->TransformPoint(mappedPoint);
    mappedMovingPoint = this->m_MovingImageTransform->TransformPoint(mappedPoint);
    if ( !this->m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint) ||
         !this->m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint) )
         sampleOk=false;
    if ( sampleOk )
      {
        localDerivative.Fill(0);

        ////// Stauffer 6/6
        // ImageDerivativesType gradient = derivativeCalculator->Evaluate(mappedPoint);
        ImageDerivativesType gradient = derivativeCalculator->Evaluate( mappedMovingPoint );
        gradientWarped = this->m_MovingImageTransform->TransformCovariantVector( gradient );
        //////

        double metricval=this->ComputeLocalContributionToMetricAndDerivative(mappedFixedPoint,mappedMovingPoint,gradientWarped,jacobian,localDerivative);

        // std::cout << " locDer size " << localDerivative.Size() << " glodir size " << derivative.Size() <<  " nvox " << thread_region.GetNumberOfPixels() << std::endl;
        if ( ! this->m_MovingImageTransform->HasLocalSupport() ) {
          derivative+=localDerivative;
        }
        else {
          // update derivative at some index
          // this requires the moving image deformation field to be
          // same size as virtual image
          OffsetValueType offset=this->m_VirtualImage->ComputeOffset(ItV.GetIndex());
          offset*=this->m_MovingImageTransform->GetNumberOfLocalParameters();
          for (unsigned int i=0; i< this->m_MovingImageTransform->GetNumberOfLocalParameters(); i++) {
            //              std::cout <<" Put in " << offset+i << " ind " << ItV.GetIndex() << std::endl;
            derivative[offset+i]=localDerivative[i];
          }
        }
        metric_sum+=metricval;
        ct++;
      }
    ++ItV;
  }
  if ( ct > 0 ) {
    std::cout << " metric_sum " << metric_sum << " ct " << ct << " thread_region " << thread_region.GetIndex() << " sz " << thread_region.GetSize() << std::endl;
    return metric_sum;
  }
  else return 0;
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
