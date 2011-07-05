/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsImageToImageObjectMetric.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDemonsImageToImageObjectMetric_txx
#define __itkDemonsImageToImageObjectMetric_txx

#include "itkDemonsImageToImageObjectMetric.h"
#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TVirtualImage >
DemonsImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::DemonsImageToImageObjectMetric()
{
}

template < class TFixedImage, class TMovingImage, class TVirtualImage >
DemonsImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::~DemonsImageToImageObjectMetric()
{
}

/*
 * GetValueAndDerivative
 */
template < class TFixedImage, class TMovingImage, class TVirtualImage >
void
DemonsImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::GetValueAndDerivative( MeasureType & value, DerivativeType & derivative)
{
  // This starts threading, and will iterate over virtual image region and
  // call GetValueAndDerivativeProcessPoint.
  this->GetValueAndDerivativeMultiThreadedInitiate( derivative );

  // Sums up results from each thread, and optionally averages them.
  // Derivative results are written directly to \c derivative.
  this->GetValueAndDerivativeMultiThreadedPostProcess( true /*doAverage*/ );

  value = this->GetValueResult();
}

/** This function computes the local voxel-wise contribution of
 *  the metric to the global integral of the metric/derivative.
 */
template < class TFixedImage, class TMovingImage, class TVirtualImage >
bool
DemonsImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           virtualPoint,
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageDerivativesType &  fixedImageDerivatives,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageDerivativesType & movingImageDerivatives,
                    MeasureType &                      metricValueReturn,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType                       threadID)
{
  //localDerivativeReturn.Fill(0); //see direct assignment below of sum

  /** Only the voxelwise contribution given the point pairs. */
  FixedImagePixelType diff = fixedImageValue - movingImageValue;
  metricValueReturn =
    vcl_fabs( diff  ) / (double) this->FixedImageDimension;

  /** For dense transforms, this returns identity */
  this->m_MovingTransform->GetJacobianWithRespectToParameters(
                                              mappedMovingPoint, m_Jacobian);

  for ( unsigned int par = 0;
          par < this->GetNumberOfLocalParameters(); par++ )
  {
    double sum = 0.0;
    for ( unsigned int dim = 0; dim < this->MovingImageDimension; dim++ )
      {
        sum += 2.0 * diff * m_Jacobian(dim, par) * movingImageDerivatives[dim];
      }
    //localDerivativeReturn[par]+=sum;
    localDerivativeReturn[par] = sum;
  }
  //  std::cout << localDerivativeReturn << std::endl;
  // Return true if the point was used in evaluation
  return true;
}

/**
 * Print out internal information about this class
 */
template < class TFixedImage, class TMovingImage, class TVirtualImage  >
void
DemonsImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "m_Jacobian size: " << m_Jacobian.size();
}

/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
DemonsImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::Initialize(void) throw ( ExceptionObject )
{
  this->Superclass::Initialize();

  //Size the jacobian matrix here so we can do it
  // just once.
  m_Jacobian.SetSize( this->VirtualImageDimension,
                      this->GetNumberOfLocalParameters() );
}

} // end namespace itk


#endif
