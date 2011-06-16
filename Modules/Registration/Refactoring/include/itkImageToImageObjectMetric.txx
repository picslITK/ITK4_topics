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
#ifndef __itkImageToImageObjectMetric_txx
#define __itkImageToImageObjectMetric_txx

#include "itkImageToImageObjectMetric.h"
//#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkPixelTraits.h"

namespace itk
{

/*
 * Constructor
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ImageToImageObjectMetric
{
  this->m_ValueAndDerivativeThreader->SetThreadedGenerateData(
    Self::GetValueAndDerivativeMultiThreadedCallback );
  this->m_ValueAndDerivativeThreader->SetHolder( this );

  /* These will be instantiated if needed in Initialize */
  m_MovingGradientCalculator = NULL;
  m_FixedGradientCalculator = NULL;
}

/*
 * Initialize. One-time initializations, i.e. once before each
 * registration loop.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::Initialize()
{

  Verify that virtual domain and moving image are the same size, when
  moving image transform is dense (e.g. deformation field). Effects
  calc of offset in UpdateDerivativeResult.
  if( PixelTraits< VirtualImagePixelType >::Dimension != 1 )
    ...
    Or can I check this at compile time? ConceptChecking has a
    HasPixelTraits thingy, but looks like it just verifies that type T
    has PixelTraits associated with it.

  Verify virtual image pixel type is scalar. Effects calc of offset
  in UpdateDerivativeResult. Right??

  /* Initialize memory for per-thread results */

  /* Memory for per-thread measure results */
  this->m_MeasurePerThread.resize( this->m_NumberOfThreads );

  this->m_NumberOfValidPointsPerThread.resize( this->m_NumberOfThreads );

  /* Allocate and initialize memory for holding derivative results.*/
  this->m_DerivativesPerThread.resize( this->m_NumberOfThreads );
  /* This size always comes from the moving image */
  unsigned long globalDerivativeSize =
    this->m_MovingImageTransform->GetNumberOfParameters();
  itkDebugMacro("ImageToImageObjectMetric::Initialize: deriv size  "
                  << globalDerivativeSize << std::endl);
  /* NOTE: this does *not* get init'ed to 0. We want to be able to use
   * a single gradient array between multiple metrics. */
  this->m_GlobalDerivative.SetSize( globalDerivativeSize );
  /* For transforms with local support, e.g. deformation field,
   * use a single derivative container that's updated by region
   * in multiple threads. */
  if ( this->m_MovingImageTransform->HasLocalSupport() )
    {
    itkDebugMacro(
      "ImageToImageObjectMetric::Initialize: tx HAS local support\n");
    for (int i=0; i<this->m_NumberOfThreads; i++)
      {
      this->m_DerivativesPerThread[i].SetData(
                                    this->m_GlobalDerivative.data_block(),
                                    this->m_GlobalDerivative.Size(),
                                    false );
      }
    }
  else
    {
    itkDebugMacro(
      "ImageToImageObjectMetric::Initialize: tx does NOT have local support\n");
    /* Global transforms get a separate derivatives container for each thread
     * that holds the result for a particular region. */
    for (int i=0; i<this->m_NumberOfThreads; i++)
      {
      this->m_DerivativesPerThread[i].SetSize( globalDerivativeSize );
      }
    }
}

/*
 * Initiate threading.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetValueAndDerivativeMultiThreadedInitiate()
{
  maybe have flag in method args to zero the GlobalDerivatives/gradientAccumulator var, which defaults to true. For use from MultiVariate metric
  Actually we may also not want to perform averaging while collecting
  post-threading results, from the MultiVariate metric, but derived class
  may call the default post-processing with averaging enabled.
  Actually we want to set these options from GetValueAndDerivative() and
  GetValue() so MultiVariate metric can access. Maybe theres an
  m_accumulateMode flag thats false by default, and can be set to true
  on a metric before MultiVariate metric calls it, and it controls several
  options liking zeroing deriv and averaging results. Dont actually know
  if MultiVariate metric will care if per-metric results are averaged.
  Or pass accumulate flag in GetValue and GetValueAndDerivative

  // Initialize some threading values that require per-iteration
  // initialization.
  for (int i=0; i<this->m_NumberOfThreads; i++)
    {
    this->m_NumberOfValidPointsPerThread[i] = 0;
    this->m_MeasurePerThread[i] = 0;
    if ( ! this->m_MovingImageTransform->HasLocalSupport() )
      {
      /* Be sure to init to 0 here, because the threader may not use
       * all the threads if the region is better split into fewer
       * subregions. */
      this->m_DerivativesPerThread[i].Fill( 0 );
      }
    }

  // Do the threaded calculations
  this->m_GetValueAndDerivativeThreader->GenerateData();

  // Count number of valid points.
  // Other post-processing can be done by calling
  // GetValueAndDerivativeMultiThreadedPostProcess, or direclty
  // in the derived class.
  this->m_NumberOfValidPoints = 0;
  for (int i=0; i<this->m_NumberOfThreads; i++)
    {
    this->m_NumberOfValidPoints += this->m_NumberOfValidPointsPerThread[i];
    }
  itkDebugMacro( "ImageToImageObjectMetric: NumberOfValidPoints: "
                 << this->m_NumberOfValidPoints );
}

/*
 * Default post-processing for threaded GetValueAndDerivative calculation.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetValueAndDerivativeMultiThreadedPostProcess( bool doAverage )

  /* For global transforms, sum the derivatives from each region. */
  if ( ! this->m_MovingImageTransform->HasLocalSupport() )
    {
    for (int i=0; i<this->m_NumberOfThreads; i++)
      {
      this->m_GlobalDerivative += this->m_DerivativesPerThread[i];
      }
    }

  /* Accumulate the metric value from threads and store */
  this->m_Value =
    NumericTraits<InternalComputationValueType>::Zero;
  for(unsigned int i=0; i< this->m_MeasurePerThread.size(); i++)
    {
    this->m_Value += this->m_MeasurePerThread[i];
    }

  if( doAverage )
    {
    if ( ! this->m_MovingImageTransform->HasLocalSupport() )
      {
      this->m_GlobalDerivative /= this->m_NumberOfValidPoints;
      }
    this->m_Value /= this->m_NumberOfValidPoints;
    }
}

/*
 * Threader callback.
 * Iterates over image region and calls derived class for calculations.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetValueAndDerivativeMultiThreadedCallback(
                const ThreaderInputObjectType& virtualImageSubRegion,
                ThreadIdType threadId,
                Self * self)
{
  /* Create an iterator over the virtual sub region */
  ImageRegionConstIteratorWithIndex<VirtualImageType>
      ItV( this->m_VirtualDomainImage, virtualImageSubRegion );

  VirtualPointType            virtualPoint;
  FixedOutputPointType        mappedFixedPoint;
  FixedImagePixelType         fixedImageValue;
  FixedImageDerivativesType   fixedImageDerivatives;
  MovingOutputPointType       mappedMovingPoint;
  MovingImagePixelType        movingImageValue;
  MovingImageDerivativesType  movingImageDerivatives;
  bool                        pointIsValid;
  MeasureType                 metricValue;
  MeasureType                 metricValueSum = 0;

  DerivativeType            localDerivative( self->GetLocalDerivativeSize() );

  /* Iterate over the sub region */
  ItV.GoToBegin();
  while( !ItV.IsAtEnd() )
  {

    self->m_VirtualImage->TransformIndexToPhysicalPoint(
                                              ItV.GetIndex(), virtualPoint);

    /* Transform the point into fixed and moving spaces, and evaluate.
     * These methods will check that the point lies within the mask if
     * one has been set, and then verify they lie in the fixed or moving
     * space as appropriate.
     * If both tests pass, the point is evaluated and 'true' is returned. */
    self->TransformAndEvaluateFixedPoint( virtualPoint,
                                          mappedFixedPoint,
                                          pointIsValid,
                                          fixedImageValue,
                                          true /*compute gradient*/,
                                          fixedImageDerivatives,
                                          threadID );
    if( pointIsValid )
      {
      self->TransformAndEvaluateMovingPoint( virtualPoint,
                                            mappedMovingPoint,
                                            pointIsValid,
                                            movingImageValue,
                                            true /*compute gradient*/,
                                            movingImageDerivatives,
                                            threadID );
      }

    /* Call the user method in derived classes to do the specific
     * calculations for value and derivative. */
    if( pointIsValid )
      {
      pointIsValid = self->GetValueAndDerivativeProcessPoint(
                   virtualPoint,
                   mappedFixedPoint, fixedImageValue, fixedImageDerivatives,
                   mappedMovingPoint, movingImageValue, movingImageDerivatives,
                   metricValue, localDerivative, threadID );
      }

    if( pointIsValid )
      {
      this->m_NumberOfValidPointsPerThread[i]++;
      metricValueSum += metricValue;
      /* Store the localDerivative. This depends on what type of
       * transform is being used. */
      UpdateDerivativeResult( localDerivative,
                              ItV.GetIndex(), threadID );
      }

    //next index
    ++ItV;
  }

  /* Store metric value result for this thread. */
  this->m_MeasurePerThread[threadID] = metricValueSum;
}

/*
 * Store the derivative results from a single point.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::UpdateDerivativeResult(  DerivativeType & localDerivative,
                           VirtualIndexType & virtualIndex,
                           ThreadIdType threadID )
{
  check this...
  DerivativeType & derivativeResult = this->m_DerivativesPerThread[threadID];
  if ( ! this->m_MovingImageTransform->HasLocalSupport() )
    {
    derivativeResult += localDerivative;
    }
  else
    {
    // update derivative at some index
    // this requires the moving image deformation field to be
    // same size as virtual image, and that VirtualImage PixelType
    // is scalar.
    OffsetValueType offset=this->m_VirtualImage->ComputeOffset(ItV.GetIndex());
    offset*=this->m_MovingImageTransform->GetNumberOfLocalParameters();
    for (unsigned int i=0;
          i < this->m_MovingImageTransform->GetNumberOfLocalParameters(); i++)
      {
      derivativeResult[offset+i]=localDerivative[i];
      }
    }
}

/*
 * Transform a point from VirtualImage domain to FixedImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateFixedPoint( VirtualPointType & point,
                      FixedImagePointType & mappedFixedPoint,
                      bool & pointIsValid,
                      FixedImagePixelType & fixedImageValue,
                      bool computeImageGradient,
                      FixedImageDerivativesType & fixedGradient,
                      unsigned int threadID) const
{
  pointIsValid = true;
  fixedImageValue = 0;

  // map the point into fixed space
  mappedFixedPoint = m_FixedImageTransform->TransformPoint( point );

  // If user provided a mask over the fixed image
  if ( m_FixedImageMask )
    {
    // Check if mapped point is within the support region of the moving image
    // mask
    pointIsValid = m_FixedImageMask->IsInside( mappedFixedPoint );
    }

  if( pointIsValid )
    {
    if ( m_FixedInterpolatorIsBSpline )
      {
      // Check if mapped point inside image buffer
      pointIsValid =
        m_FixedBSplineInterpolator->IsInsideBuffer( mappedFixedPoint );
      if ( pointIsValid )
        {
        if( computeImageGradient )
          {
          /* The assumption here is that calling EvaluateValueAndDerivative
           * is faster than separately calling Evaluate and EvaluateDerivative,
           * because of setup for the bspline evaluation. */
          m_FixedBSplineInterpolator->EvaluateValueAndDerivative(
                                                        mappedFixedPoint,
                                                        fixedImageValue,
                                                        fixedGradient,
                                                        threadID);
          }
        else
          {
          fixedImageValue =
            m_FixedBSplineInterpolator->Evaluate(mappedFixedPoint, threadID);
          }
        }
      }
    else
      {
      // Check if mapped point is inside image buffer
      pointIsValid = m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint);
      if ( pointIsValid )
        {
        fixedImageValue = m_FixedInterpolator->Evaluate(mappedFixedPoint);
        if( computeImageGradient )
          {
          this->ComputeFixedImageDerivatives(mappedFixedPoint,
                                             fixedGradient,
                                             threadID);
          }
        }
      }

    if( pointIsValid && computeImageGradient )
      {
      //Transform the gradient into the virtual domain. We compute gradient
      // in the fixed and moving domains and then transform to virtual, because
      // first warping the fixed and moving images to virtual domain and then
      // calculating gradient would create threading issues with the creation
      // of the warped image, unless they were created during initilization.
      // But, pre-warping the images would be inefficient when a mask or
      // sampling is used to compute only a subset of points.
      fixedGradient =
        m_FixedImageTransform->TransformCovariantVector( fixedGradient );
      }
    }
}

/*
 * Transform a point from VirtualImage domain to MovingImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateMovingPoint( VirtualPointType & point,
                      MovingImagePointType & mappedMovingPoint,
                      bool & pointIsValid,
                      MovingImagePixelType & movingImageValue,
                      bool computeImageGradient,
                      MovingImageDerivativesType & movingGradient,
                      unsigned int threadID) const
{
  pointIsValid = true;
  movingImageValue = 0;
  mappedMovingPoint = m_MovingImageTransform->TransformPoint( point );

  // If user provided a mask over the Moving image
  if ( m_MovingImageMask )
    {
    // Check if mapped point is within the support region of the moving image
    // mask
    pointIsValid = m_MovingImageMask->IsInside( mappedMovingPoint );
    }

  if( pointIsValid )
    {
    if ( m_MovingInterpolatorIsBSpline )
      {
      // Check if mapped point inside image buffer
      pointIsValid =
        m_MovingBSplineInterpolator->IsInsideBuffer( mappedMovingPoint );
      if ( pointIsValid )
        {
        if( computeImageGradient )
          {
          /* The assumption here is that calling EvaluateValueAndDerivative
           * is faster than separately calling Evaluate and EvaluateDerivative,
           * because of setup for the bspline evaluation. */
          m_MovingBSplineInterpolator->EvaluateValueAndDerivative(
                                                        mappedMovingPoint,
                                                        movingImageValue,
                                                        movingGradient,
                                                        threadID);
          // Transform into the virtual space
          gradient =
            m_MovingImageTransform->TransformCovariantVector( gradient );
          }
        else
          {
          movingImageValue =
            m_MovingBSplineInterpolator->Evaluate(mappedMovingPoint, threadID);
          }
        }
      }
    else
      {
      // Check if mapped point inside image buffer
      pointIsValid = m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint);
      if ( pointIsValid )
        {
        movingImageValue = m_MovingInterpolator->Evaluate(mappedMovingPoint);
        if( computeImageGradient )
          {
          this->ComputeMovingImageDerivatives(mappedMovingPoint,
                                              movingGradient,
                                              threadID);
          }
        }
      }
    if( pointIsValid && computeImageGradient )
      {
      // Transform into the virtual space
      movingGradient =
        m_MovingImageTransform->TransformCovariantVector( movingGradient );
      }
    }
}

/*
 * Compute image derivatives for a Fixed point.
 * NOTE: This doesn't transform result into virtual space. For that,
 * see TransformAndEvaluateFixedPoint
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeFixedImageDerivatives(
                              const FixedImagePointType & mappedPoint,
                              FixedImageDerivativesType & gradient,
                              unsigned int threadID) const
{
  if ( m_FixedInterpolatorIsBSpline )
    {
    // Computed Fixed image gradient using derivative BSpline kernel.
    gradient = m_FixedBSplineInterpolator->EvaluateDerivative(mappedPoint,
                                                              threadID);
    }
  else
    {
    if ( m_UseGaussianGradient )
      {
      ContinuousIndex< double, FixedImageDimension > tempIndex;
      m_FixedImage->TransformPhysicalPointToContinuousIndex(mappedPoint,
                                                             tempIndex);
      FixedImageIndexType mappedIndex;
      mappedIndex.CopyWithRound(tempIndex);
      gradient = m_FixedGaussianGradientImage->GetPixel(mappedIndex);
      }
    else
      {
      // if not using the gradient image
      gradient = m_FixedGradientCalculator->Evaluate(mappedPoint);
      }
    }
}

/*
 * Compute image derivatives for a moving point.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeMovingImageDerivatives(
                              const MovingImagePointType & mappedPoint,
                              MovingImageDerivativesType & gradient,
                              unsigned int threadID) const
{
  if ( m_MovingInterpolatorIsBSpline )
    {
    // Computed moving image gradient using derivative BSpline kernel.
    gradient = m_MovingBSplineInterpolator->EvaluateDerivative(mappedPoint,
                                                               threadID);
    }
  else
    {
    if ( m_UseGaussianGradient )
      {
      ContinuousIndex< double, MovingImageDimension > tempIndex;
      m_MovingImage->TransformPhysicalPointToContinuousIndex(mappedPoint,
                                                             tempIndex);
      MovingImageIndexType mappedIndex;
      mappedIndex.CopyWithRound(tempIndex);
      gradient = m_MovingGaussianGradientImage->GetPixel(mappedIndex);
      }
    else
      {
      // if not using the gradient image
      gradient = m_MovingGradientCalculator->Evaluate(mappedPoint);
      }
    }

  //Transform the gradient into the virtual domain.
  gradient = m_MovingImageTransform->TransformCovariantVector( gradient );
}

/*
 * Get the size of localDerivative
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
SizeValueType
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetLocalDerivativeSize()
{
  switch( this->m_DerivativeSource )
    {
    case Moving:
      return this->m_MovingImageTransform->GetNumberOfLocalParameters() );
      break;
    case Fixed:
      return this->m_FixedImageTransform->GetNumberOfLocalParameters() );
      break;
    case Both:
      /* TODO: this isn't completely settled, what to do here. Derived
       * classes may need to override to provide different behavior. */
      return VirtualImageDimension;
      break;
    default:
      itkExceptionMacro("Invalid DerivativeSource enum: "
                         << this->m_DerivativeSource );
    }
  return 0;
}

}//namespace itk
