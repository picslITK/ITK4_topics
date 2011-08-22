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
#ifndef __itkImageToImageObjectMetric_hxx
#define __itkImageToImageObjectMetric_hxx

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
::ImageToImageObjectMetric()
{
  /* Instantiate the default threader, and set the callback and user-data
   * holder for threading */
  this->m_ValueAndDerivativeThreader = ValueAndDerivativeThreaderType::New();
  this->m_ValueAndDerivativeThreader->SetThreadedGenerateData(
    Self::GetValueAndDerivativeMultiThreadedCallback );
  this->m_ValueAndDerivativeThreader->SetHolder( this );
  this->m_NumberOfThreads =
      this->m_ValueAndDerivativeThreader->GetNumberOfThreads();
  this->m_ThreadingMemoryHasBeenInitialized = false;

  m_FixedImage = NULL;
  m_MovingImage = NULL;
  m_VirtualDomainImage = NULL;
  m_VirtualDomainRegionHasBeenSet = false;

  /* Both transforms default to an identity transform */
  m_FixedTransform = FixedIdentityTransformType::New();
  m_MovingTransform = MovingIdentityTransformType::New();

  /* Interpolators. Default to linear. */
  m_FixedInterpolator = FixedLinearInterpolatorType::New();
  m_MovingInterpolator = MovingLinearInterpolatorType::New();

  m_PrecomputeImageGradient = true;
  /* These will be instantiated if needed in Initialize */
  m_FixedGaussianGradientImage = NULL;
  m_MovingGaussianGradientImage = NULL;
  m_MovingGradientCalculator = NULL;
  m_FixedGradientCalculator = NULL;

  m_PreWarpImages = false;
  m_FixedWarpedImage = NULL;
  m_MovingWarpedImage = NULL;

  m_FixedImageMask = NULL;
  m_MovingImageMask = NULL;

  m_DerivativeResult = NULL;

  m_FixedTransformCanUseTransformIndex = false;
  m_MovingTransformCanUseTransformIndex = false;
}

/*
 * Initialize. One-time initializations, i.e. once before each
 * registration loop.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::Initialize() throw ( itk::ExceptionObject )
{

  /* Verify things are connected */
  if ( !m_FixedImage )
    {
    itkExceptionMacro(<< "FixedImage is not present");
    }
  if ( !m_MovingImage )
    {
    itkExceptionMacro(<< "MovingImage is not present");
    }
  if ( !m_FixedTransform )
    {
    itkExceptionMacro(<< "FixedTransform is not present");
    }
  if ( !m_MovingTransform )
    {
    itkExceptionMacro(<< "MovingTransform is not present");
    }

  // If the image is provided by a source, update the source.
  if ( m_MovingImage->GetSource() )
    {
    m_MovingImage->GetSource()->Update();
    }

  // If the image is provided by a source, update the source.
  if ( m_FixedImage->GetSource() )
    {
    m_FixedImage->GetSource()->Update();
    }

  /* If a virtual image has not been set, create one from fixed image */
  if( ! m_VirtualDomainImage )
    {
    itkDebugMacro("Creating VirtualDomainImage from FixedImage");
    /* This instantiation will fail at compilation if user has provided
     * a different type for VirtualImage in the template parameters. */
    m_VirtualDomainImage = FixedImageType::New();
    /* Graft the virtual image onto the fixed, to conserve memory. */
    m_VirtualDomainImage->Graft( m_FixedImage );
    /* If user hasn't already provided a region, get the buffered region
     * from fixed image. */
    if( ! m_VirtualDomainRegionHasBeenSet )
      {
      /* Make sure we set this before assigning it to threader below */
      this->SetVirtualDomainRegion( m_VirtualDomainImage->GetBufferedRegion());
      }
    }

  if( m_MovingTransform->HasLocalSupport() )
    {
    /* Verify that virtual domain and displacement field are the same size
    * and in the same physical space. Handles CompositeTransform by checking
    * if first applied transform is DisplacementFieldTransform */
    this->VerifyDisplacementFieldSizeAndPhysicalSpace();

    /* Verify virtual image pixel type is scalar. Effects calc of offset
    in StoreDerivativeResult.
    NOTE:  Can this be checked at compile time? ConceptChecking has a
    HasPixelTraits class, but looks like it just verifies that type T
    has PixelTraits associated with it, and not a particular value. */
    if( PixelTraits< VirtualImagePixelType >::Dimension != 1 )
      {
      itkExceptionMacro("VirtualImagePixelType must be scalar. "
                        "Dimensionality is " <<
                        PixelTraits< VirtualImagePixelType >::Dimension );
      }
    }//Moving transform is dense

  /*
   * Determine number of threads that will be used
   */

  /* Assign the virtual image region to the threader. Do this before
   * calling DetermineNumberOfThreasToUse. */
  this->m_ValueAndDerivativeThreader->SetOverallRegion(
                                            this->m_VirtualDomainRegion );

  /* Determine how many threads will be used, given the set OveralRegion.
   * The threader uses
   * its SplitRequestedObject method to split the image region over threads,
   * and it may decide to that using fewer than the number of available
   * threads is more effective. */
  this->m_NumberOfThreads = this->m_ValueAndDerivativeThreader->
                                      DetermineNumberOfThreadsToUse();

  /* Clear this flag to force initialization of threading memory
   * in GetValueAndDerivativeMultiThreadedInitiate. */
  this->m_ThreadingMemoryHasBeenInitialized = false;

  /*
   * Interpolators and image gradients.
   * Do this AFTER we have the
   * proper value in m_NumberOfThreads.
   */

  /* Inititialize interpolators.
   * Note that if image pre-warping is enabled, interpolators will be
   * updated as needed after pre-warped images are created. */
  m_FixedInterpolator->SetInputImage( m_FixedImage );
  m_MovingInterpolator->SetInputImage( m_MovingImage );

  /* Setup for image gradient calculations.
   * Instantiate a central difference derivative calculator
   * if appropriate. */

  /* Fixed */
    if( !m_PrecomputeImageGradient )
      {
      m_FixedGradientCalculator = FixedGradientCalculatorType::New();
      m_FixedGradientCalculator->UseImageDirectionOn();
      m_FixedGradientCalculator->SetInputImage(this->m_FixedImage);
      }
  /* Moving */
    if( ! m_PrecomputeImageGradient )
      {
      m_MovingGradientCalculator = MovingGradientCalculatorType::New();
      m_MovingGradientCalculator->UseImageDirectionOn();
      m_MovingGradientCalculator->SetInputImage(this->m_MovingImage);
      }

  /* If user set to use pre-calculated gradient image,
   * then we need to calculate the gradient images. */
  if ( m_PrecomputeImageGradient )
    {
    ComputeGaussianGradient();
    }

  /* Initialize resample image filters for pre-warping images if
   * option is set. */
  if( this->m_PreWarpImages )
    {
    if( this->m_FixedImageMask || this->m_MovingImageMask )
      {
      itkExceptionMacro("Use of m_PreWarpImages with image masks is not "
                        "yet supported." );
      }
    if( this->m_PrecomputeImageGradient )
      {
      itkExceptionMacro("Use of m_PreWarpImages with m_PrecomputeImageGradient "
                        "is not supported (currently, at least). ");
      }
    m_MovingWarpResampleImageFilter = MovingWarpResampleImageFilterType::New();
    m_MovingWarpResampleImageFilter->SetUseReferenceImage( true );
    m_MovingWarpResampleImageFilter->SetReferenceImage(
                                               this->GetVirtualDomainImage() );
    m_MovingWarpResampleImageFilter->SetNumberOfThreads( this->m_NumberOfThreads );
    m_MovingWarpResampleImageFilter->SetTransform( this->GetMovingTransform() );
    m_MovingWarpResampleImageFilter->SetInput( this->GetMovingImage() );

    m_FixedWarpResampleImageFilter = FixedWarpResampleImageFilterType::New();
    m_FixedWarpResampleImageFilter->SetUseReferenceImage( true );
    m_FixedWarpResampleImageFilter->SetReferenceImage(
                                               this->GetVirtualDomainImage() );
    m_FixedWarpResampleImageFilter->SetNumberOfThreads( this->m_NumberOfThreads );
    m_FixedWarpResampleImageFilter->SetTransform( this->GetFixedTransform() );
    m_FixedWarpResampleImageFilter->SetInput( this->GetFixedImage() );
    }
  else
    {
    /* Free memory if allocated from a previous run */
    m_MovingWarpedImage = NULL;
    m_FixedWarpedImage = NULL;
    }

  /* Check if the transforms can use TransformIndex.
   * Determine this during initialization and store
   * result because the check can involve several comparisons. */
  /*
   TODO sort out. Interface will likely be simplified...
  m_FixedTransformCanUseTransformIndex =
    m_FixedTransform->CanUseTransformIndex(
      this->GetVirtualDomainImage()->GetRequestedRegion(),
      this->GetVirtualDomainImage()->GetOrigin(),
      this->GetVirtualDomainImage()->GetSpacing(),
      this->GetVirtualDomainImage()->GetDirection() );

  m_MovingTransformCanUseTransformIndex =
    m_MovingTransform->CanUseTransformIndex(
      this->GetVirtualDomainImage()->GetRequestedRegion(),
      this->GetVirtualDomainImage()->GetOrigin(),
      this->GetVirtualDomainImage()->GetSpacing(),
      this->GetVirtualDomainImage()->GetDirection() );
  */
}

/*
 * Initiate threaded evaluation.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetValueAndDerivativeMultiThreadedInitiate( DerivativeType & derivativeReturn )
{
  //Initialize threading memory if this is the first time
  // in here since a call to Initialize, or if user has passed
  // in a different object for the results (why might they do that?).
  if( ! this->m_ThreadingMemoryHasBeenInitialized ||
      &derivativeReturn != this->m_DerivativeResult )
    {
    this->InitializeThreadingMemory( derivativeReturn );
    }

  //Initialization required for each iteration.
  InitializeForIteration();

  // Do the threaded evaluation. This will
  // call GetValueAndDerivativeMultiThreadedCallback, which
  // iterates over virtual domain region and calls derived class'
  // GetValueAndDerivativeProcessPoint.
  this->m_ValueAndDerivativeThreader->GenerateData();

  // Determine the total number of points used during calculations.
  CollectNumberOfValidPoints();

  // To collect the results from each thread into final values
  // the derived class can call GetValueAndDerivativeMultiThreadedPostProcess,
  // or do their own processing.
}

/*
 * InitializeThreadingMemory
 * We have a separate method that's called during the first call to evaluate
 * the metric after Initialize has been called. This is so we can handle
 * m_DerivativeResult as a raw pointer, obtained from user input.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::InitializeThreadingMemory( DerivativeType & derivativeReturn )
{
  /* Point our results object to the object provided by user. */
  this->m_DerivativeResult = &derivativeReturn;

  /* Per-thread results */
  this->m_MeasurePerThread.resize( this->m_NumberOfThreads );
  this->m_NumberOfValidPointsPerThread.resize( this->m_NumberOfThreads );
  this->m_DerivativesPerThread.resize( this->m_NumberOfThreads );
  /* This one is intermediary, for getting per-point results. */
  this->m_LocalDerivativesPerThread.resize( this->m_NumberOfThreads );
  /* Per-thread pre-allocated Jacobian objects for efficiency */
  this->m_MovingTransformJacobianPerThread.resize( this->m_NumberOfThreads );

  /* This size always comes from the moving image */
  unsigned long globalDerivativeSize =
    this->m_MovingTransform->GetNumberOfParameters();
  itkDebugMacro("ImageToImageObjectMetric::Initialize: deriv size  "
                  << globalDerivativeSize << std::endl);
  /* NOTE: this does *not* get init'ed to 0 here. */
  if( this->m_DerivativeResult->GetSize() != globalDerivativeSize )
    {
    this->m_DerivativeResult->SetSize( globalDerivativeSize );
    }
  for (ThreadIdType i=0; i<this->m_NumberOfThreads; i++)
    {
    /* Allocate intermediary per-thread storage used to get results from
     * derived classes */
    this->m_LocalDerivativesPerThread[i].SetSize(
                                          this->GetNumberOfLocalParameters() );
    this->m_MovingTransformJacobianPerThread[i].SetSize(
                                          this->VirtualImageDimension,
                                          this->GetNumberOfLocalParameters() );
    /* For transforms with local support, e.g. displacement field,
     * use a single derivative container that's updated by region
     * in multiple threads. */
    if ( this->m_MovingTransform->HasLocalSupport() )
      {
      itkDebugMacro(
        "ImageToImageObjectMetric::Initialize: tx HAS local support\n");
        /* Set each per-thread object to point to m_DerivativeResult */
        this->m_DerivativesPerThread[i].SetData(
                                      this->m_DerivativeResult->data_block(),
                                      this->m_DerivativeResult->Size(),
                                      false );
      }
    else
      {
      itkDebugMacro(
      "ImageToImageObjectMetric::Initialize: tx does NOT have local support\n");
      /* Global transforms get a separate derivatives container for each thread
       * that holds the result over a particular image region. */
        this->m_DerivativesPerThread[i].SetSize( globalDerivativeSize );
      }
    }
  /* This will be true until next call to Initialize */
  this->m_ThreadingMemoryHasBeenInitialized = true;
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::InitializeForIteration()
{
  /* Initialize some threading values that require per-iteration
   * initialization. */
  for (ThreadIdType i=0; i<this->m_NumberOfThreads; i++)
    {
    this->m_NumberOfValidPointsPerThread[i] = 0;
    this->m_MeasurePerThread[i] = 0;
    if ( ! this->m_MovingTransform->HasLocalSupport() )
      {
      /* Be sure to init to 0 here, because the threader may not use
       * all the threads if the region is better split into fewer
       * subregions. */
      this->m_DerivativesPerThread[i].Fill( 0 );
      }
    }

  /* Clear derivative final result. This will
   * require an option to skip for use with multivariate metric. */
  this->m_DerivativeResult->Fill( 0 );

  /* Pre-warp the images if set to do so. */
  if( this->m_PreWarpImages )
    {
    this->PreWarpImages();
    }
}

/*
 * Collect number of valid points, after threading is completed.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::CollectNumberOfValidPoints()
{
  /* Count number of valid points.
   * Other post-processing can be done by calling
   * GetValueAndDerivativeMultiThreadedPostProcess, or direclty
   * in the derived class. */
  this->m_NumberOfValidPoints = 0;
  for (ThreadIdType i=0; i<this->m_NumberOfThreads; i++)
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
{
  /* For global transforms, sum the derivatives from each region. */
  if ( ! this->m_MovingTransform->HasLocalSupport() )
    {
    for (ThreadIdType i=0; i<this->m_NumberOfThreads; i++)
      {
      *(this->m_DerivativeResult) += this->m_DerivativesPerThread[i];
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
    if ( ! this->m_MovingTransform->HasLocalSupport() )
      {
      *(this->m_DerivativeResult) /= this->m_NumberOfValidPoints;
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
                ThreadIdType threadID,
                Self * self)
{
  /* Create an iterator over the virtual sub region */
  ImageRegionConstIteratorWithIndex<VirtualImageType>
      ItV( self->m_VirtualDomainImage, virtualImageSubRegion );

  VirtualPointType            virtualPoint;
  FixedOutputPointType        mappedFixedPoint;
  FixedImagePixelType         mappedFixedPixelValue;
  FixedImageGradientType      mappedFixedImageGradient;
  MovingOutputPointType       mappedMovingPoint;
  MovingImagePixelType        mappedMovingPixelValue;
  MovingImageGradientType     mappedMovingImageGradient;
  bool                        pointIsValid = false;
  MeasureType                 metricValueResult;
  MeasureType                 metricValueSum = 0;

  /* Get pre-allocated local results object. This actually provides very
   * little benefit, since this only gets called once for each thread. However
   * if we get up to the hundres of threads, it might have an impact */
  DerivativeType & localDerivativeResult =
                                   self->m_LocalDerivativesPerThread[threadID];

  /* Iterate over the sub region */
  //ItV.GoToBegin();
  //while( !ItV.IsAtEnd() )
  for( ItV.GoToBegin(); !ItV.IsAtEnd(); ++ItV )
  {
    /* Get the virtual point */
    self->m_VirtualDomainImage->TransformIndexToPhysicalPoint(
                                              ItV.GetIndex(), virtualPoint);

    /* Transform the point into fixed and moving spaces, and evaluate.
     * Different behavior with pre-warping enabled is handled transparently.
     * Do this in a try block to catch exceptions and print more useful info
     * then we otherwise get when exceptions are caught in MultiThreader. */
    try
      {
      self->TransformAndEvaluateFixedPoint( ItV.GetIndex(),
                                        virtualPoint,
                                        self->GetGradientSourceIncludesFixed(),
                                        mappedFixedPoint,
                                        mappedFixedPixelValue,
                                        mappedFixedImageGradient,
                                        pointIsValid );
      }
    catch( ExceptionObject & exc )
      {
      //NOTE: there must be a cleaner way to do this:
      std::string msg("Caught exception: \n");
      msg += exc.what();
      ExceptionObject err(__FILE__, __LINE__, msg);
      throw err;
      }

    if( !pointIsValid )
      {
      continue;
      }

    try
      {
      self->TransformAndEvaluateMovingPoint( ItV.GetIndex(),
                                      virtualPoint,
                                      self->GetGradientSourceIncludesMoving(),
                                      mappedMovingPoint,
                                      mappedMovingPixelValue,
                                      mappedMovingImageGradient,
                                      pointIsValid );
      }
    catch( ExceptionObject & exc )
      {
      std::string msg("Caught exception: \n");
      msg += exc.what();
      ExceptionObject err(__FILE__, __LINE__, msg);
      throw err;
      }

    if( !pointIsValid )
      {
      continue;
      }

    /* Call the user method in derived classes to do the specific
     * calculations for value and derivative. */
    try
      {
      pointIsValid = self->GetValueAndDerivativeProcessPoint(
                                     virtualPoint,
                                     mappedFixedPoint, mappedFixedPixelValue,
                                     mappedFixedImageGradient,
                                     mappedMovingPoint, mappedMovingPixelValue,
                                     mappedMovingImageGradient,
                                     metricValueResult, localDerivativeResult,
                                     threadID );
      }
    catch( ExceptionObject & exc )
      {
      //NOTE: there must be a cleaner way to do this:
      std::string msg("Exception in GetValueAndDerivativeProcessPoint:\n");
      msg += exc.what();
      ExceptionObject err(__FILE__, __LINE__, msg);
      throw err;
      }

    if( !pointIsValid )
      {
      continue;
      }

    /* Assign the results */
    self->m_NumberOfValidPointsPerThread[ threadID ]++;
    metricValueSum += metricValueResult;
    /* Store the result. The behavior depends on what type of
     * transform is being used. */
    self->StoreDerivativeResult( localDerivativeResult,
                            ItV.GetIndex(), threadID );
  } //loop over region

  /* Store metric value result for this thread. */
  self->m_MeasurePerThread[threadID] = metricValueSum;
}

/*
 * Store the derivative results from a single point.
 * \warning See warning in definition, regarding overriding
 * or not using this class.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::StoreDerivativeResult(  DerivativeType & derivative,
                           const VirtualIndexType & virtualIndex,
                           ThreadIdType threadID )
{
  DerivativeType & derivativeResult = this->m_DerivativesPerThread[threadID];
  if ( ! this->m_MovingTransform->HasLocalSupport() )
    {
    derivativeResult += derivative;
    }
  else
    {
    // Dense transform, e.g. displacement field.
    // update derivative at some index
    // this requires the moving image displacement field to be
    // same size as virtual image, and that VirtualImage PixelType
    // is scalar.
    try
      {
      OffsetValueType offset =
        this->m_VirtualDomainImage->ComputeOffset(virtualIndex);
      offset *= this->m_MovingTransform->GetNumberOfLocalParameters();
      for (unsigned int i=0;
            i < this->m_MovingTransform->GetNumberOfLocalParameters(); i++)
        {
        /* Be sure to *add* here and not assign. Required for proper behavior
         * with multi-variate metric. */
        derivativeResult[offset+i] += derivative[i];
        }
      }
    catch( ExceptionObject & exc )
      {
      std::string msg("Caught exception: \n");
      msg += exc.what();
      ExceptionObject err(__FILE__, __LINE__, msg);
      throw err;
      }
    }
}

#if 0 //(***************************************************************
/*
 * Evaluate at an index within pre-warped images.
 */
virtual void TransformAndEvaluateWarpedFixedImageAtIndex(
                         const VirtualIndexType & index,
                         const VirtualPointType & point,
                         const bool computeImageGradient,
                         FixedImagePointType & mappedFixedPoint,
                         FixedImagePixelType & mappedFixedPixelValue,
                         FixedImageGradientType & mappedFixedImageGradient,
                         bool & pointIsValid ) const
{
  pointIsValid = true;

  /* Get the point in moving and fixed space for use below */
  if( m_FixedTransformCanUseTransformIndex )
    {
    mappedFixedPoint =
              self->m_FixedTransform->TransformIndex( index );
    }
  else
    {
    mappedFixedPoint =
                  self->m_FixedTransform->TransformPoint( virtualPoint );
    }

  // If user provided a mask over the fixed image
  if ( m_FixedImageMask )
    {
    // Check if mapped point is within the support region of the fixed image
    // mask
    pointIsValid = m_FixedImageMask->IsInside( mappedFixedPoint );
    }

  if( ! pointIsValid )
    {
    return;
    }

  /* Get the pixel values at this index */
  mappedFixedPixelValue = m_FixedWarpedImage->GetPixel( index );

  if( computeImageGradient )
    {
    ComputeFixedImageGradientAtIndex( index, mappedFixedImageGradient );
    }
}

/*
 * Evaluate at an index within pre-warped images.
 */
virtual void TransformAndEvaluateWarpedMovingImageAtIndex(
                         const VirtualIndexType & index,
                         const VirtualPointType & virtualPoint,
                         const bool computeImageGradient,
                         MovingImagePointType & mappedMovingPoint,
                         MovingImagePixelType & mappedMovingPixelValue,
                         MovingImageGradientType & mappedMovingImageGradient,
                         bool & pointIsValid ) const;
{
  pointIsValid = true;

  if( m_MovingTransformCanUseTransformIndex )
    {
    mappedMovingPoint =
             self->m_MovingTransform->TransformIndex( index );
    }
  else
    {
    mappedMovingPoint =
                  self->m_MovingTransform->TransformPoint( virtualPoint );
    }

  /* Check mask */
  if ( m_MovingImageMask )
    {
    // Check if mapped point is within the support region of the moving image
    // mask
    pointIsValid = m_MovingImageMask->IsInside( mappedMovingPoint );
    }

  if( ! pointIsValid )
    {
    return;
    }

  /* Get the pixel values at this index */
  mappedMovingPixelValue = m_MovingWarpedImage->GetPixel( index );

  if( computeImageGradient )
    {
    ComputeMovingImageGradientAtIndex( index, mappedMovingImageGradient );
    }
}

#endif //***************************************************************

/*
 * Transform a point from VirtualImage domain to FixedImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateFixedPoint(
                         const VirtualIndexType & index,
                         const VirtualPointType & point,
                         const bool computeImageGradient,
                         FixedImagePointType & mappedFixedPoint,
                         FixedImagePixelType & mappedFixedPixelValue,
                         FixedImageGradientType & mappedFixedImageGradient,
                         bool & pointIsValid ) const
{
  pointIsValid = true;
  mappedFixedPixelValue = NumericTraits<FixedImagePixelType>::Zero;

  // map the point into fixed space
  if( m_FixedTransformCanUseTransformIndex )
    {
    mappedFixedPoint = m_FixedTransform->TransformIndex( index );
    }
  else
    {
    mappedFixedPoint = m_FixedTransform->TransformPoint( point );
    }

  // check against the mask if one is assigned
  if ( m_FixedImageMask )
    {
    // Check if mapped point is within the support region of the fixed image
    // mask
    pointIsValid = m_FixedImageMask->IsInside( mappedFixedPoint );
    if( ! pointIsValid )
      {
      return;
      }
    }

  // Check if mapped point is inside image buffer
  pointIsValid = m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint);
  if( ! pointIsValid )
    {
    return;
    }

  if( m_PreWarpImages )
    {
    /* Get the pixel values at this index */
    mappedFixedPixelValue = m_FixedWarpedImage->GetPixel( index );

    if( computeImageGradient )
      {
      ComputeFixedImageGradientAtIndex( index, mappedFixedImageGradient );
      }
    }
  else
    {
    mappedFixedPixelValue = m_FixedInterpolator->Evaluate(mappedFixedPoint);
    if( computeImageGradient )
      {
      this->ComputeFixedImageGradient( mappedFixedPoint,
                                       mappedFixedImageGradient );
      //Transform the gradient into the virtual domain. We compute gradient
      // in the fixed and moving domains and then transform to virtual, because
      // first warping the fixed and moving images to virtual domain and then
      // calculating gradient would create threading issues with the creation
      // of the warped image, unless they were created during initilization.
      // But, pre-warping the images would be inefficient when a mask or
      // sampling is used to compute only a subset of points.
      mappedFixedImageGradient =
        m_FixedTransform->TransformCovariantVector( mappedFixedImageGradient,
                                                         mappedFixedPoint );
      }
    }
}

/*
 * Transform a point from VirtualImage domain to MovingImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateMovingPoint(
                         const VirtualIndexType & index,
                         const VirtualPointType & point,
                         const bool computeImageGradient,
                         MovingImagePointType & mappedMovingPoint,
                         MovingImagePixelType & mappedMovingPixelValue,
                         MovingImageGradientType & mappedMovingImageGradient,
                         bool & pointIsValid ) const
{
  pointIsValid = true;
  mappedMovingPixelValue = NumericTraits<MovingImagePixelType>::Zero;

  // map the point into fixed space
  if( m_MovingTransformCanUseTransformIndex )
    {
    mappedMovingPoint = m_MovingTransform->TransformIndex( index );
    }
  else
    {
    mappedMovingPoint = m_MovingTransform->TransformPoint( point );
    }

  // check against the mask if one is assigned
  if ( m_MovingImageMask )
    {
    // Check if mapped point is within the support region of the fixed image
    // mask
    pointIsValid = m_MovingImageMask->IsInside( mappedMovingPoint );
    if( ! pointIsValid )
      {
      return;
      }
    }

  // Check if mapped point is inside image buffer
  pointIsValid = m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint);
  if( ! pointIsValid )
    {
    return;
    }

  if( m_PreWarpImages )
    {
    /* Get the pixel values at this index */
    mappedMovingPixelValue = m_MovingWarpedImage->GetPixel( index );

    if( computeImageGradient )
      {
      ComputeMovingImageGradientAtIndex( index, mappedMovingImageGradient );
      }
    }
  else
    {
    mappedMovingPixelValue = m_MovingInterpolator->Evaluate(mappedMovingPoint);
    if( computeImageGradient )
      {
      this->ComputeMovingImageGradient( mappedMovingPoint,
                                       mappedMovingImageGradient );
      mappedMovingImageGradient =
        m_MovingTransform->TransformCovariantVector( mappedMovingImageGradient,
                                                         mappedMovingPoint );
      }
    }
}

#if 0 //*************************************************************
/*
 * Transform a point from VirtualImage domain to MovingImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateMovingPoint(
                      const VirtualIndexType & index,
                      const VirtualPointType & point,
                      MovingImagePointType & mappedMovingPoint,
                      bool & pointIsValid,
                      MovingImagePixelType & movingImageValue,
                      bool computeImageGradient,
                      MovingImageGradientType & movingImageGradient ) const
{
  pointIsValid = true;
  movingImageValue = 0;

  if( this->m_MovingTransform->HasLocalSupport() )
    {
    mappedMovingPoint = m_MovingTransform->TransformIndex( index );
    }
  else
    {
    mappedMovingPoint = m_MovingTransform->TransformPoint( point );
    }

  // If user provided a mask over the Moving image
  if ( m_MovingImageMask )
    {
    // Check if mapped point is within the support region of the moving image
    // mask
    pointIsValid = m_MovingImageMask->IsInside( mappedMovingPoint );
    }

  if( ! pointIsValid )
    {
    return;
    }

  // Check if mapped point inside image buffer
  pointIsValid = m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint);
  if ( pointIsValid )
    {
    movingImageValue = m_MovingInterpolator->Evaluate(mappedMovingPoint);
    if( computeImageGradient )
      {
      this->ComputeMovingImageGradient( mappedMovingPoint,
                                        movingImageGradient );
      }
    }

  if( pointIsValid && computeImageGradient )
    {
    // Transform into the virtual space. See TransformAndEvaluateFixedPoint.
    movingImageGradient =
      m_MovingTransform->TransformCovariantVector( movingImageGradient,
                                                        mappedMovingPoint );
    }
}


#endif //******************************************************************/

/*
 * Compute image derivatives for a Fixed point.
 * NOTE: This doesn't transform result into virtual space. For that,
 * see TransformAndEvaluateFixedPoint
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeFixedImageGradient(
                              const FixedImagePointType & mappedPoint,
                              FixedImageGradientType & gradient ) const
{
  if ( m_PrecomputeImageGradient )
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

/*
 * Compute image derivatives for a moving point.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeMovingImageGradient(
                              const MovingImagePointType & mappedPoint,
                              MovingImageGradientType & gradient ) const
{
  if ( m_PrecomputeImageGradient )
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

/*
 * Compute fixed warped image derivatives for an index at virtual domain.
 * NOTE: This doesn't transform result into virtual space. For that,
 * see TransformAndEvaluateFixedPoint
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeFixedImageGradientAtIndex(
                              const VirtualIndexType & index,
                              FixedImageGradientType & gradient ) const
{
  if ( m_PrecomputeImageGradient )
    {
      gradient = m_FixedGaussianGradientImage->GetPixel(index);
    }
  else
    {
    // if not using the gradient image
    gradient = m_FixedGradientCalculator->EvaluateAtIndex(index);
    }
}

/*
 * Compute image derivatives for a moving point.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeMovingImageGradientAtIndex(
                              const VirtualIndexType & index,
                              MovingImageGradientType & gradient ) const
{
  if ( m_PrecomputeImageGradient )
    {
    gradient = m_MovingGaussianGradientImage->GetPixel(index);
    }
  else
    {
    // if not using the gradient image
    gradient = m_MovingGradientCalculator->EvaluateAtIndex(index);
    }
}

/*
 * Pre-warp images.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::PreWarpImages()
{
  /* Call Modified to make sure the filter recalculates the output. We haven't
   * changed any settings, but we assume the transform parameters have changed,
   * e.g. while used during registration. */
  m_MovingWarpResampleImageFilter->Modified();
  m_MovingWarpResampleImageFilter->Update();
  m_MovingWarpedImage = m_MovingWarpResampleImageFilter->GetOutput();

  m_FixedWarpResampleImageFilter->Modified();
  m_FixedWarpResampleImageFilter->Update();
  m_FixedWarpedImage = m_FixedWarpResampleImageFilter->GetOutput();

  /* Point the interpolators to the warped images.
   * We should try to skip this for efficiency. Will be possible if
   * ResampleImageFilter always returns the same image pointer after
   * its first update, or if it can be set to allocate output during init. */
  /* No need to call Modified here on the calculators */
  m_FixedGradientCalculator->SetInputImage( m_FixedWarpedImage );
  m_MovingGradientCalculator->SetInputImage( m_MovingWarpedImage );
}

/*
 * ComputeGaussianGradient
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeGaussianGradient()
{
  /* Fixed */
  {
  typename FixedGradientImageFilterType::Pointer
    gradientFilter = FixedGradientImageFilterType::New();
  gradientFilter->SetInput(m_FixedImage);
  const typename FixedImageType::SpacingType & spacing = m_FixedImage
                                                          ->GetSpacing();
  double maximumSpacing = 0.0;
  for ( unsigned int i = 0; i < FixedImageDimension; i++ )
    {
    if ( spacing[i] > maximumSpacing )
      {
      maximumSpacing = spacing[i];
      }
    }
  gradientFilter->SetSigma(maximumSpacing);
  gradientFilter->SetNormalizeAcrossScale(true);
  gradientFilter->SetNumberOfThreads(m_NumberOfThreads);
  gradientFilter->SetUseImageDirection(true);
  gradientFilter->Update();

  m_FixedGaussianGradientImage = gradientFilter->GetOutput();
  }
  /* Moving */
  {
  typename MovingGradientImageFilterType::Pointer
    gradientFilter = MovingGradientImageFilterType::New();

  gradientFilter->SetInput(m_MovingImage);

  const typename MovingImageType::SpacingType & spacing = m_MovingImage
                                                          ->GetSpacing();
  double maximumSpacing = 0.0;
  for ( unsigned int i = 0; i < MovingImageDimension; i++ )
    {
    if ( spacing[i] > maximumSpacing )
      {
      maximumSpacing = spacing[i];
      }
    }
  gradientFilter->SetSigma(maximumSpacing);
  gradientFilter->SetNormalizeAcrossScale(true);
  gradientFilter->SetNumberOfThreads(m_NumberOfThreads);
  gradientFilter->SetUseImageDirection(true);
  gradientFilter->Update();

  m_MovingGaussianGradientImage = gradientFilter->GetOutput();
  }
}

/*
 * UpdateParameters
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::UpdateTransformParameters( DerivativeType & derivative,
                             ParametersValueType factor )
{
  if( derivative.GetSize() != this->GetNumberOfParameters() )
    {
    itkExceptionMacro("derivative is not the proper size");
    }
  this->m_MovingTransform->UpdateTransformParameters( derivative, factor );
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::SetNumberOfThreads( ThreadIdType number )
{
  if( number != this->m_NumberOfThreads )
    {
    this->m_NumberOfThreads = number;
    this->m_ValueAndDerivativeThreader->SetNumberOfThreads( number );
    this->Modified();
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::CreateVirtualDomainImage( VirtualSpacingType & spacing,
                            VirtualOriginType & origin,
                            VirtualDirectionType & direction,
                            VirtualSizeType & size,
                            VirtualIndexType & index )
{
  VirtualRegionType region;
  region.SetSize( size );
  region.SetIndex( index );
  CreateVirtualDomainImage( spacing, origin, direction, region );
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::CreateVirtualDomainImage( VirtualSpacingType & spacing,
                            VirtualOriginType & origin,
                            VirtualDirectionType & direction,
                            VirtualRegionType & region )
{
  this->m_VirtualDomainImage = VirtualImageType::New();
  this->m_VirtualDomainImage->SetSpacing( spacing );
  this->m_VirtualDomainImage->SetOrigin( origin );
  this->m_VirtualDomainImage->SetDirection( direction );
  this->m_VirtualDomainImage->SetRegions( region );
  this->m_VirtualDomainImage->Allocate();
  this->m_VirtualDomainImage->FillBuffer( 0 );
  this->SetVirtualDomainRegion(region);
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
const typename ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >::VirtualSpacingType
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetVirtualDomainSpacing( void )
{
  if( this->m_VirtualDomainImage )
    {
    return this->m_VirtualDomainImage->GetSpacing();
    }
  else
    {
    itkExceptionMacro("m_VirtualDomainImage is undefined. Cannot "
                      " return spacing. ");
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
const typename ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >::VirtualDirectionType
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetVirtualDomainDirection( void )
{
  if( this->m_VirtualDomainImage )
    {
    return this->m_VirtualDomainImage->GetDirection();
    }
  else
    {
    itkExceptionMacro("m_VirtualDomainImage is undefined. Cannot "
                      " return direction. ");
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
const typename ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >::VirtualOriginType
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetVirtualDomainOrigin( void )
{
  if( this->m_VirtualDomainImage )
    {
    return this->m_VirtualDomainImage->GetOrigin();
    }
  else
    {
    itkExceptionMacro("m_VirtualDomainImage is undefined. Cannot "
                      " return origin. ");
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::SetVirtualDomainImage(VirtualImageType* image)
{
  this->m_VirtualDomainImage = image;
  this->SetVirtualDomainRegion( image->GetBufferedRegion() );
}

/*
 * Verify a displacement field and virtual image are in the same space.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::VerifyDisplacementFieldSizeAndPhysicalSpace()
{

//TODO replace with common external method. Possibly use Transform::CanUseTransformIndex as stop-gap.

#if 0
  /* Verify that virtual domain and displacement field are the same size
   * and in the same physical space.
   * Effects transformation, and calc of offset in StoreDerivativeResult.
   * If it's a composite transform and the displacement field is the first
   * to be applied (i.e. the most recently added), then it has to be
   * of the same size, otherwise not. But actually at this point, if
   * a CompositeTransform has local support, it means all its sub-transforms
   * have local support. So they should all be displacement fields, so just
   * verify that the first one is at least.
   * Eventually we'll want a method in Transform something like a
   * GetInputDomainSize to check this cleanly. */
  MovingTransformType* transform;
  transform = this->m_MovingTransform.GetPointer();
  /* If it's a CompositeTransform, get the last transform (1st applied). */
  MovingCompositeTransformType* comptx =
               dynamic_cast< MovingCompositeTransformType * > ( transform );
  if( comptx != NULL )
    {
    transform = comptx->GetBackTransform().GetPointer();
    }
  /* Check that it's a DisplacementField type, the only type we expect
   * at this point */
  MovingDisplacementFieldTransformType* deftx =
          dynamic_cast< MovingDisplacementFieldTransformType * >( transform );
  if( deftx == NULL )
    {
    itkExceptionMacro("Expected m_MovingTransform to be of type "
                      "DisplacementFieldTransform" );
    }
  typedef typename MovingDisplacementFieldTransformType::DisplacementFieldType
                                                                    FieldType;
  typename FieldType::Pointer field = deftx->GetDisplacementField();
  typename FieldType::RegionType
    fieldRegion = field->GetBufferedRegion();
  VirtualRegionType virtualRegion =
                            m_VirtualDomainImage->GetBufferedRegion();
  if( virtualRegion.GetSize() != fieldRegion.GetSize() ||
      virtualRegion.GetIndex() != fieldRegion.GetIndex() )
    {
    itkExceptionMacro("Virtual domain and moving transform displacement field"
                      " must have the same size and index for "
                      " LargestPossibleRegion."
                      << std::endl << "Virtual size/index: "
                      << virtualRegion.GetSize() << " / "
                      << virtualRegion.GetIndex() << std::endl
                      << "Displacement field size/index: "
                      << fieldRegion.GetSize() << " / "
                      << fieldRegion.GetIndex() << std::endl );
    }

    /* check that the image occupy the same physical space, and that
     * each index is at the same physical location.
     * this code is from ImageToImageFilter */

    /* tolerance for origin and spacing depends on the size of pixel
     * tolerance for directions a fraction of the unit cube. */
    const double coordinateTol
      = 1.0e-6 * m_VirtualDomainImage->GetSpacing()[0]; // use first dimension spacing
    const double directionTol = 1.0e-6;

    if ( !m_VirtualDomainImage->GetOrigin().GetVnlVector().
               is_equal( field->GetOrigin().GetVnlVector(), coordinateTol ) ||
         !m_VirtualDomainImage->GetSpacing().GetVnlVector().
               is_equal( field->GetSpacing().GetVnlVector(), coordinateTol ) ||
         !m_VirtualDomainImage->GetDirection().GetVnlMatrix().as_ref().
               is_equal( field->GetDirection().GetVnlMatrix(), directionTol ) )
      {
      std::ostringstream originString, spacingString, directionString;
      originString << "m_VirtualDomainImage Origin: "
                   << m_VirtualDomainImage->GetOrigin()
                   << ", DisplacementField Origin: " << field->GetOrigin()
                   << std::endl;
      spacingString << "m_VirtualDomainImage Spacing: "
                    << m_VirtualDomainImage->GetSpacing()
                    << ", DisplacementField Spacing: "
                    << field->GetSpacing() << std::endl;
      directionString << "m_VirtualDomainImage Direction: "
                      << m_VirtualDomainImage->GetDirection()
                      << ", DisplacementField Direction: "
                      << field->GetDirection() << std::endl;
      itkExceptionMacro(<< "m_VirtualDomainImage and DisplacementField do not "
                        << "occupy the same physical space! "
                        << std::endl
                        << originString.str() << spacingString.str()
                        << directionString.str() );
      }
#endif
}

}//namespace itk

#endif
