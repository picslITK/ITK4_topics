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
#ifndef __itkVectorImageToVectorImageObjectMetric_hxx
#define __itkVectorImageToVectorImageObjectMetric_hxx

#include "itkVectorImageToVectorImageObjectMetric.h"
//#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkPixelTraits.h"

namespace itk
{

/*
 * Constructor
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::VectorImageToVectorImageObjectMetric()
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
  /* m_[Fixed|Moving]Interpolator get assigned to this during Initialize
   * if they're found to be bspline interpolators. */
  //m_FixedBSplineInterpolator = NULL;
  //m_MovingBSplineInterpolator = NULL;

  m_PrecomputeImageGradient = false;
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
}

/*
 * Initialize. One-time initializations, i.e. once before each
 * registration loop.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
    /* Verify that virtual domain and displacement field are the same size.
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
    MovingCompositeTransformType* comptx =
                 dynamic_cast< MovingCompositeTransformType * > ( transform );
    if( comptx != NULL )
      {
      transform = comptx->GetBackTransform().GetPointer();
      }
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

  /* Setup for image gradient calculations and possible BSpline interpolation.
    Check if the interpolator is of type BSplineInterpolateImageFunction.
    If so, we can make use of its EvaluateDerivatives method.
    Otherwise, we instantiate an external central difference
    derivative calculator.
  */
  /* Fixed */
  //FixedBSplineInterpolatorType *fixedTestPtr =
  //      dynamic_cast< FixedBSplineInterpolatorType * >
  //        ( this->m_FixedInterpolator.GetPointer() );
  //if ( !fixedTestPtr )
    //{
    /* It's a BSpline interpolator */
    m_FixedInterpolatorIsBSpline = false;
    //m_FixedBSplineInterpolator = NULL;
    if( !m_PrecomputeImageGradient )
      {
      m_FixedGradientCalculator = FixedGradientCalculatorType::New();
      m_FixedGradientCalculator->UseImageDirectionOn();
      m_FixedGradientCalculator->SetInputImage(this->m_FixedImage);
      }
    itkDebugMacro("Fixed Interpolator is not BSpline");
    //}
  //else
  //  {
    /* It's not a BSpline interpolator */
    //m_FixedInterpolatorIsBSpline = true;
    //m_FixedBSplineInterpolator = fixedTestPtr;
    //m_FixedBSplineInterpolator->SetNumberOfThreads(m_NumberOfThreads);
    //m_FixedBSplineInterpolator->UseImageDirectionOn();
    //m_FixedGradientCalculator = NULL;
    //itkDebugMacro("Fixed Interpolator is BSpline");
    //}
  /* Moving */
  //MovingBSplineInterpolatorType *movingTestPtr =
  //      dynamic_cast< MovingBSplineInterpolatorType * >
  //        ( this->m_MovingInterpolator.GetPointer() );
  //if ( !movingTestPtr )
    //{
    m_MovingInterpolatorIsBSpline = false;
    //m_MovingBSplineInterpolator = NULL;
    if( ! m_PrecomputeImageGradient )
      {
      m_MovingGradientCalculator = MovingGradientCalculatorType::New();
      m_MovingGradientCalculator->UseImageDirectionOn();
      m_MovingGradientCalculator->SetInputImage(this->m_MovingImage);
      }
    itkDebugMacro("Moving Interpolator is not BSpline");
    //}
  //else
  //  {
  //  m_MovingInterpolatorIsBSpline = true;
  //  m_MovingBSplineInterpolator = movingTestPtr;
  //  m_MovingBSplineInterpolator->SetNumberOfThreads(m_NumberOfThreads);
  //  m_MovingBSplineInterpolator->UseImageDirectionOn();
  //  m_MovingGradientCalculator = NULL;
  //  itkDebugMacro("Moving Interpolator is BSpline");
  //  }
  /* If user set to use pre-calculated gradient image, and either interpolator
   * isn't bspline, then we need to calculate the gradient images. */
  if ( m_PrecomputeImageGradient &&
       !( m_FixedInterpolatorIsBSpline && m_MovingInterpolatorIsBSpline ) )
    {
    //NOTE: This could be broken into separate fixed- and moving-gaussian
    // calculations if users will be having one interpolator as bspline
    // and not the other.
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
    if( this->m_MovingInterpolatorIsBSpline || m_FixedInterpolatorIsBSpline )
      {
      /* Using BSpline with pre-warp option requires handling the two different
       * ways bspline interpolators are used: 1) image gradient calcs; 2)
       * interpolation. */
      itkExceptionMacro("Use of BSpline interpolators is not currently "
                        " with m_PreWarpImages option." );
      }
    m_MovingWarpResampleImageFilter = MovingWarpResampleImageFilterType::New();
    //m_MovingWarpResampleImageFilter->SetUseReferenceImage( true );
    //m_MovingWarpResampleImageFilter->SetReferenceImage(
    //                                           this->GetVirtualDomainImage() );
    m_MovingWarpResampleImageFilter->SetSize(
      this->GetVirtualDomainImage()->GetLargestPossibleRegion().GetSize() );
    m_MovingWarpResampleImageFilter->SetOutputOrigin(
      this->GetVirtualDomainImage()->GetOrigin() );
    m_MovingWarpResampleImageFilter->SetOutputSpacing(
      this->GetVirtualDomainImage()->GetSpacing() );
    m_MovingWarpResampleImageFilter->SetOutputDirection(
      this->GetVirtualDomainImage()->GetDirection() );
    m_MovingWarpResampleImageFilter->SetNumberOfThreads( this->m_NumberOfThreads );
    m_MovingWarpResampleImageFilter->SetTransform( this->GetMovingTransform() );
    m_MovingWarpResampleImageFilter->SetInput( this->GetMovingImage() );

    m_FixedWarpResampleImageFilter = FixedWarpResampleImageFilterType::New();
    //m_FixedWarpResampleImageFilter->SetUseReferenceImage( true );
    //m_FixedWarpResampleImageFilter->SetReferenceImage(
    //                                           this->GetVirtualDomainImage() );
    m_FixedWarpResampleImageFilter->SetSize(
      this->GetVirtualDomainImage()->GetLargestPossibleRegion().GetSize() );
    m_FixedWarpResampleImageFilter->SetOutputOrigin(
      this->GetVirtualDomainImage()->GetOrigin() );
    m_FixedWarpResampleImageFilter->SetOutputSpacing(
      this->GetVirtualDomainImage()->GetSpacing() );
    m_FixedWarpResampleImageFilter->SetOutputDirection(
      this->GetVirtualDomainImage()->GetDirection() );
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
}

/*
 * Initiate threaded evaluation.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
  this->m_ValueAndDerivativeThreader->StartThreadedExecution();
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
  itkDebugMacro("VectorImageToVectorImageObjectMetric::Initialize: deriv size  "
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
        "VectorImageToVectorImageObjectMetric::Initialize: tx HAS local support\n");
        /* Set each per-thread object to point to m_DerivativeResult */
        this->m_DerivativesPerThread[i].SetData(
                                      this->m_DerivativeResult->data_block(),
                                      this->m_DerivativeResult->Size(),
                                      false );
      }
    else
      {
      itkDebugMacro(
      "VectorImageToVectorImageObjectMetric::Initialize: tx does NOT have local support\n");
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
  itkDebugMacro( "VectorImageToVectorImageObjectMetric: NumberOfValidPoints: "
                 << this->m_NumberOfValidPoints );
}

/*
 * Default post-processing for threaded GetValueAndDerivative calculation.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
  FixedImagePixelType         fixedImageValue;
  FixedImageGradientType   fixedImageGradient;
  MovingOutputPointType       mappedMovingPoint;
  MovingImagePixelType        movingImageValue;
  MovingImageGradientType  movingImageGradient;
  bool                        pointIsValid = false;
  MeasureType                 metricValueResult;
  MeasureType                 metricValueSum = 0;

  /* Get pre-allocated local results object. This actually provides very
   * little benefit, since this only gets called once for each thread. However
   * if we get up to the hundres of threads, it might have an impact */
  DerivativeType & localDerivativeResult =
                                   self->m_LocalDerivativesPerThread[threadID];

  /* Iterate over the sub region */
  ItV.GoToBegin();
  while( !ItV.IsAtEnd() )
  {
    /* Get the virtual point */
    self->m_VirtualDomainImage->TransformIndexToPhysicalPoint(
                                              ItV.GetIndex(), virtualPoint);

    if( self->m_PreWarpImages )
      {
      try
        {
        self->EvaluatePreWarpedImagesAtIndex( ItV.GetIndex(),
                                        virtualPoint,
                                        true, /* compute gradient */
                                        fixedImageValue,
                                        movingImageValue,
                                        fixedImageGradient,
                                        movingImageGradient,
                                        pointIsValid,
                                        threadID );
        /* Get the point in moving and fixed space for use below */
        mappedFixedPoint =
                      self->m_FixedTransform->TransformPoint( virtualPoint );
        mappedMovingPoint =
                      self->m_MovingTransform->TransformPoint( virtualPoint );
        }
      catch( ExceptionObject & exc )
        {
        //NOTE: there must be a cleaner way to do this. We want to add
        // the this filename and line number to give user more useful
        // information about where the exception was generated.
        std::string msg("Caught exception: \n");
        msg += exc.what();
        ExceptionObject err(__FILE__, __LINE__, msg);
        throw err;
        }
      }
    else
      {
      /* Transform the point into fixed and moving spaces, and evaluate.
       * These methods will check that the point lies within the mask if
       * one has been set, and then verify they lie in the fixed or moving
       * space as appropriate.
       * If both tests pass, the point is evaluated and pointIsValid is
       * returned as \c true.
       * Do this in a try block to catch exceptions and print more useful info
       * then we otherwise get when exceptions are caught in MultiThreader. */
      try
        {
        self->TransformAndEvaluateFixedPoint( virtualPoint,
                                              mappedFixedPoint,
                                              pointIsValid,
                                              fixedImageValue,
                                              true /*compute gradient*/,
                                              fixedImageGradient,
                                              threadID );
        if( pointIsValid )
          {
          self->TransformAndEvaluateMovingPoint( virtualPoint,
                                                mappedMovingPoint,
                                                pointIsValid,
                                                movingImageValue,
                                                true /*compute gradient*/,
                                                movingImageGradient,
                                                threadID );
          }
        }
      catch( ExceptionObject & exc )
        {
        //NOTE: there must be a cleaner way to do this:
        std::string msg("Caught exception: \n");
        msg += exc.what();
        ExceptionObject err(__FILE__, __LINE__, msg);
        throw err;
        }
      }


    /* Call the user method in derived classes to do the specific
     * calculations for value and derivative. */
    try
      {
      if( pointIsValid )
        {
        pointIsValid = self->GetValueAndDerivativeProcessPoint(
                 virtualPoint,
                 mappedFixedPoint, fixedImageValue, fixedImageGradient,
                 mappedMovingPoint, movingImageValue, movingImageGradient,
                 metricValueResult, localDerivativeResult, threadID );
        }
      }
    catch( ExceptionObject & exc )
      {
      //NOTE: there must be a cleaner way to do this:
      std::string msg("Exception in GetValueAndDerivativeProcessPoint:\n");
      msg += exc.what();
      ExceptionObject err(__FILE__, __LINE__, msg);
      throw err;
      }

    /* Assign the results */
    if( pointIsValid )
      {
      self->m_NumberOfValidPointsPerThread[ threadID ]++;
      metricValueSum += metricValueResult;
      /* Store the result. This depends on what type of
       * transform is being used. */
      self->StoreDerivativeResult( localDerivativeResult,
                              ItV.GetIndex(), threadID );
      }

    //next index
    ++ItV;
  }

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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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

/*
 * Evaluate at an index within pre-warped images.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::EvaluatePreWarpedImagesAtIndex( const VirtualIndexType & index,
                                 const VirtualPointType & virtualPoint,
                                 const bool computeImageGradient,
                                 FixedImagePixelType & fixedImageValue,
                                 MovingImagePixelType & movingImageValue,
                                 FixedImageGradientType & fixedImageGradient,
                                 MovingImageGradientType & movingImageGradient,
                                 bool & pointIsValid,
                                 const ThreadIdType threadID ) const
{
  /* For now, this is always true. When we enable mask usage with pre-warping,
   * we'll check the mask here. */
  pointIsValid = true;

  /* TODO check mask */

  /* Get the pixel values at this index */
  fixedImageValue = m_FixedWarpedImage->GetPixel( index );
  movingImageValue = m_MovingWarpedImage->GetPixel( index );

  if( computeImageGradient )
    {
    ComputeFixedImageGradient( virtualPoint, fixedImageGradient, threadID );
    ComputeMovingImageGradient( virtualPoint, movingImageGradient, threadID );
    }
}


/*
 * Transform a point from VirtualImage domain to FixedImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateFixedPoint(const VirtualPointType & point,
                      FixedImagePointType & mappedFixedPoint,
                      bool & pointIsValid,
                      FixedImagePixelType & fixedImageValue,
                      bool computeImageGradient,
                      FixedImageGradientType & fixedImageGradient,
                      ThreadIdType threadID) const
{
  pointIsValid = true;
  fixedImageValue = 0;

  // map the point into fixed space
  mappedFixedPoint = m_FixedTransform->TransformPoint( point );

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

    // Check if mapped point is inside image buffer
    pointIsValid = m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint);
    if ( pointIsValid )
      {
      fixedImageValue = m_FixedInterpolator->Evaluate(mappedFixedPoint);
      if( computeImageGradient )
        {
        this->ComputeFixedImageGradient(mappedFixedPoint,
                                           fixedImageGradient,
                                           threadID);
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

    for (unsigned int com=0;
      com < m_FixedImage->GetNumberOfComponentsPerPixel(); com++)
      {
      FixedGradientVectorPixelType cGrad;
      for (unsigned int dim=0; dim < FixedImageDimension; dim++)
        {
        cGrad[dim] = fixedImageGradient(com,dim);
        }
      FixedGradientVectorPixelType tGrad
        = m_FixedTransform->TransformCovariantVector( cGrad, mappedFixedPoint );

      for (unsigned int dim=0; dim < FixedImageDimension; dim++)
        {
        fixedImageGradient(com,dim) = tGrad[dim];
        }
      }
    }
}

/*
 * Transform a point from VirtualImage domain to MovingImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateMovingPoint(const VirtualPointType & point,
                      MovingImagePointType & mappedMovingPoint,
                      bool & pointIsValid,
                      MovingImagePixelType & movingImageValue,
                      bool computeImageGradient,
                      MovingImageGradientType & movingImageGradient,
                      ThreadIdType threadID) const
{
  pointIsValid = true;
  movingImageValue = 0;

  mappedMovingPoint = m_MovingTransform->TransformPoint( point );

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
        this->ComputeMovingImageGradient(mappedMovingPoint,
                                            movingImageGradient,
                                            threadID);
        }
      }

  if( pointIsValid && computeImageGradient )
    {
    // Transform into the virtual space. See TransformAndEvaluateFixedPoint.
    for (unsigned int com=0;
      com < m_MovingImage->GetNumberOfComponentsPerPixel(); com++)
      {
      MovingGradientVectorPixelType cGrad;
      for (unsigned int dim=0; dim < MovingImageDimension; dim++)
        {
        cGrad[dim] = movingImageGradient(com,dim);
        }
      MovingGradientVectorPixelType tGrad
        = m_MovingTransform->TransformCovariantVector( cGrad, mappedMovingPoint );

      for (unsigned int dim=0; dim < MovingImageDimension; dim++)
        {
        movingImageGradient(com,dim) = tGrad[dim];
        }
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeFixedImageGradient(
                              const FixedImagePointType & mappedPoint,
                              FixedImageGradientType & gradient,
                              ThreadIdType ) const
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeMovingImageGradient(
                              const MovingImagePointType & mappedPoint,
                              MovingImageGradientType & gradient,
                              ThreadIdType ) const
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
 * Pre-warp images.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
  if( m_FixedInterpolatorIsBSpline )
    {
    itkExceptionMacro(
                "BSpline interpolators with pre-warping is not yet supported.");
    }
  else
    {
    /* No need to call Modified here on the calculators */
    m_FixedGradientCalculator->SetInputImage( m_FixedWarpedImage );
    m_MovingGradientCalculator->SetInputImage( m_MovingWarpedImage );
    }
}

/*
 * ComputeGaussianGradient
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::ComputeGaussianGradient()
{

  /* Fixed */
  /*
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
  */

  /* Moving */
  /*
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
  */
}

/*
 * UpdateParameters
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
const typename VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >::VirtualSpacingType
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
const typename VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >::VirtualDirectionType
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
const typename VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >::VirtualOriginType
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
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
VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::SetVirtualDomainImage(VirtualImageType* image)
{
  this->m_VirtualDomainImage = image;
  this->SetVirtualDomainRegion( image->GetBufferedRegion() );
}

/**
 * Compute the maximum voxel shift when a change of deltaParameters is applied to
 * the current transform.
 *
 * isMovingTransform specifies whether it is about the transform from VirtualImage
 * to MovingImage or to FixedImage.
 *
 */
template< class TFixedImage, class TMovingImage, class TVirtualImage >
double
VectorImageToVectorImageObjectMetric< TFixedImage, TMovingImage, TVirtualImage >
::ComputeMaximumVoxelShift(bool isMovingTransform, ParametersType deltaParameters)
{
  double voxelShift;
  double distance;
  unsigned int numSamples = 0;

  this->ComputeVirtualImageCornerPoints();
  numSamples = m_VirtualImageCornerPoints.size();

  ParametersType oldParameters(deltaParameters.size());
  ParametersType newParameters(deltaParameters.size());

  //itkExceptionMacro("Verify that MovingTransformType and FixedTransformType are correct below...");

  if (isMovingTransform)
    {
    oldParameters = m_MovingTransform->GetParameters();
    for (unsigned int p=0; p<oldParameters.size(); p++)
      {
      newParameters[p] = oldParameters[p] + deltaParameters[p];
      }

    MovingImagePointType point, oldMappedPoint, newMappedPoint;
    ContinuousIndex<double, MovingImageDimension> oldMappedIndex, newMappedIndex, diffIndex;

    // find max shift by checking each sample point
    voxelShift = 0.0;
    for (unsigned int c=0; c<numSamples; c++)
      {
      point = m_VirtualImageCornerPoints[c];

      oldMappedPoint = m_MovingTransform->TransformPoint(point);

      m_MovingTransform->SetParameters(newParameters);
      newMappedPoint = m_MovingTransform->TransformPoint(point);
      m_MovingTransform->SetParameters(oldParameters); //restore parameters

      this->m_MovingImage->TransformPhysicalPointToContinuousIndex(oldMappedPoint, oldMappedIndex);
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex(newMappedPoint, newMappedIndex);

      distance = 0.0;
      //itkExceptionMacro("Verify MovingImageDimension below...");
      for (unsigned int d=0; d<MovingImageDimension; d++)
        {
        diffIndex[d] = oldMappedIndex[d] - newMappedIndex[d];
        distance += diffIndex[d] * diffIndex[d];
        }
      distance = vcl_sqrt(distance);

      if ( voxelShift < distance )
        {
        voxelShift = distance;
        }
      } // end for numSamples
    } // end of if isMovingTransform
  else
    {
    oldParameters = m_FixedTransform->GetParameters();
    for (unsigned int p=0; p<oldParameters.size(); p++)
      {
      newParameters[p] = oldParameters[p] + deltaParameters[p];
      }
    MovingTransformPointer newFixedTransform = NULL; //FixedTransformType::New();
    newFixedTransform->SetParameters(newParameters);

    FixedImagePointType point, oldMappedPoint, newMappedPoint;
    ContinuousIndex<double, FixedImageDimension> oldMappedIndex, newMappedIndex, diffIndex;

    // find max shift by checking each sample point
    voxelShift = 0.0;
    for (unsigned int c=0; c<numSamples; c++)
      {
      point = m_VirtualImageCornerPoints[c];

      oldMappedPoint = m_FixedTransform->TransformPoint(point);

      m_FixedTransform->SetParameters(newParameters);
      newMappedPoint = m_FixedTransform->TransformPoint(point);
      m_FixedTransform->SetParameters(oldParameters); //restore parameters

      this->m_FixedImage->TransformPhysicalPointToContinuousIndex(oldMappedPoint, oldMappedIndex);
      this->m_FixedImage->TransformPhysicalPointToContinuousIndex(newMappedPoint, newMappedIndex);

      distance = 0.0;
      //itkExceptionMacro("Verify FixedImageDimension below...");
      for (unsigned int d=0; d<FixedImageDimension; d++)
        {
        diffIndex[d] = oldMappedIndex[d] - newMappedIndex[d];
        distance += diffIndex[d] * diffIndex[d];
        }
      distance = vcl_sqrt(distance);

      if ( voxelShift < distance )
        {
        voxelShift = distance;
        }
      } // end of for numSamples
    } // end of else isMovingTransform

  return voxelShift;
}

/**
 * Get the physical coordinates of image corners
 */
template< class TFixedImage, class TMovingImage, class TVirtualImage >
void
VectorImageToVectorImageObjectMetric< TFixedImage, TMovingImage, TVirtualImage >
::ComputeVirtualImageCornerPoints( )
{
  if (m_VirtualImageCornerPoints.size() > 0) return;

  m_VirtualImageCornerPoints.clear();

  //itkExceptionMacro("Verify that VirtualImageDimension is correct below...");

  VirtualRegionType region = m_VirtualDomainImage->GetLargestPossibleRegion();
  VirtualIndexType firstCorner = region.GetIndex();
  VirtualIndexType corner;
  VirtualPointType point;

  VirtualSizeType size = region.GetSize();
  int cornerNumber = 1 << VirtualImageDimension; // 2^ImageDimension

  for( int i=0; i<cornerNumber; i++ )
    {
    int bit;
    for (int d=0; d<VirtualImageDimension; d++)
      {
      bit = (int) (( i & (1 << d) ) != 0); // 0 or 1
      corner[d] = firstCorner[d] + bit * (size[d] - 1);
      }

    m_VirtualDomainImage->TransformIndexToPhysicalPoint(corner, point);
    m_VirtualImageCornerPoints.push_back(point);
    }
}

/** Compute parameter scales
 *
 * Define a transform as Y = T(P,X) where X, Y are coordinates, and P are
 * parameters. Let p be one of the parameters. Introduce an intermediate
 * variable q = s * p where s is the scaling factor. The p column of the
 * transform Jacobian is written as d_Y/d_p.
 *
 * The scales are estimated in a way such that a unit change of each
 * parameter would yield roughly similar voxel shifts. That is, d_Y/d_q
 * is roughly same for all scaled parameters.
 *
 * After applying scale s, d_Y/d_q = d_Y/d_p * (1/s). Therefore we may
 * choose s = d_Y/d_p.

 * And we have the metric derivative d_metric/d_q = d_metric/d_p * (1/s).
 *
 * Let us consider the impact of s to the scales in ITK optimizers. In the
 * ITK class GradientDescentOptimizer, the step of parameter p is
 *   step_p = (d_metric/d_p) / old_opt_scale * learningRate .
 *
 * With q = s * p, the step for parameter q will be
 *   step_q = (d_metric/d_q) / old_opt_scale * learningRate
 *          = (d_metric/d_p) / s / old_opt_scale * learningRate .
 *
 * Since step_p = step_q / s, we have
 *   step_p = (d_metric/d_p) / s / s / old_opt_scale * learningRate .
 *
 * With parameter scale s, the optimizer scale for p is changed to
 *   new_opt_scale = s * s * old_opt_scale .
 * where s = d_Y/d_p.
 *
 * Here we use s = 1/MaximumVoxelShift for each parameter p.
 */
template< class TFixedImage, class TMovingImage, class TVirtualImage >
void
VectorImageToVectorImageObjectMetric< TFixedImage, TMovingImage, TVirtualImage >
::EstimateScales(bool isMovingTransform, ParametersType &parameterScales)
{
  double maxShift;
  unsigned int numPara = parameterScales.size();

  ParametersType deltaParameters(numPara);
  for (unsigned int j=0; j<numPara; j++)
    {
    deltaParameters[j] = 0;
    }

  // To avoid division-by-zero, assign a small value for parameters
  // whose impact on voxel shift is zero. This small value may be
  // the minimum non-zero shift.
  double minNonZeroShift = NumericTraits<double>::max();

  // compute variation for each transform parameter
  for (unsigned int i=0; i<numPara; i++)
    {
    deltaParameters[i] = 1;
    maxShift = this->ComputeMaximumVoxelShift(isMovingTransform, deltaParameters);
    deltaParameters[i] = 0;

    parameterScales[i] = maxShift;
    if ( maxShift != 0 && maxShift < minNonZeroShift )
      {
      minNonZeroShift = maxShift;
      }
    }

  for (unsigned int i=0; i<numPara; i++)
    {
    if (parameterScales[i] == 0)
      {
      parameterScales[i] = minNonZeroShift * minNonZeroShift;
      //use NumericTraits<double>::max() if we don't want to
      //optimize this parameter at all?
      }
    else
      {
      parameterScales[i] *= parameterScales[i];
      }
    }

}

}//namespace itk

#endif
