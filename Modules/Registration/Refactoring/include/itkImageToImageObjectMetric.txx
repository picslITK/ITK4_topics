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
  m_FixedImageTransform = NULL;
  m_MovingImage = NULL;
  m_MovingImageTransform = NULL;
  m_VirtualDomainImage = NULL;
  m_VirtualDomainRegionHasBeenSet = false;

  /* Interpolators. Default to linear. */
  m_FixedInterpolator = FixedLinearInterpolatorType::New();
  m_MovingInterpolator = MovingLinearInterpolatorType::New();
  /* m_[Fixed|Moving]Interpolator get assigned to this during Initialize
   * if they're found to be bspline interpolators. */
  m_FixedBSplineInterpolator = NULL;
  m_MovingBSplineInterpolator = NULL;

  m_PrecomputeImageGradient = true;
  /* These will be instantiated if needed in Initialize */
  m_FixedGaussianGradientImage = NULL;
  m_MovingGaussianGradientImage = NULL;
  m_MovingGradientCalculator = NULL;
  m_FixedGradientCalculator = NULL;

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
  if ( !m_FixedImageTransform )
    {
    itkExceptionMacro(<< "FixedImageTransform is not present");
    }
  if ( !m_MovingImageTransform )
    {
    itkExceptionMacro(<< "MovingImageTransform is not present");
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
    /* If user hasn't already provided a region, get the requested region
     * from fixed image. */
    if( ! m_VirtualDomainRegionHasBeenSet )
      {
      /* Make sure we set this before assigning it to threader below */
      this->SetVirtualDomainRegion( m_VirtualDomainImage->GetBufferedRegion());
      }
    }

  if( m_MovingImageTransform->HasLocalSupport() )
    {
    /* Verify that virtual domain and deformation field are the same size.
     * Effects transformation, and calc of offset in StoreDerivativeResult.
     * If it's a composite transform and the deformation field is the first
     * to be applied (i.e. the most recently added), then it has to be
     * of the same size, otherwise not. But actually at this point, if
     * a CompositeTransform has local support, it means all its sub-transforms
     * have local support. So they should all be deformation fields, so just
     * verify that the first one is at least.
     * Eventually we'll want a method in Transform something like a
     * GetInputDomainSize to check this cleanly. */
    MovingTransformType* transform;
    transform = this->m_MovingImageTransform.GetPointer();
    MovingCompositeTransformType* comptx =
                 dynamic_cast< MovingCompositeTransformType * > ( transform );
    if( comptx != NULL )
      {
      transform = comptx->GetBackTransform().GetPointer();
      }
    MovingDeformationFieldTransformType* deftx =
            dynamic_cast< MovingDeformationFieldTransformType * >( transform );
    if( deftx == NULL )
      {
      itkExceptionMacro("Expected m_MovingImageTransform to be of type "
                        "DeformationFieldTransform" );
      }
    typedef typename MovingDeformationFieldTransformType::DeformationFieldType
                                                                      FieldType;
    typename FieldType::Pointer field = deftx->GetDeformationField();
    typename FieldType::RegionType
      fieldRegion = field->GetBufferedRegion();
    VirtualRegionType virtualRegion =
                              m_VirtualDomainImage->GetBufferedRegion();
    if( virtualRegion.GetSize()  != fieldRegion.GetSize() ||
        virtualRegion.GetIndex() != fieldRegion.GetIndex() )
      {
      itkExceptionMacro("Virtual domain and moving transform deformation field"
                        " must have the same size and index for "
                        " LargestPossibleRegion."
                        << std::endl << "Virtual size/index: "
                        << virtualRegion.GetSize() << " / "
                        << virtualRegion.GetIndex() << std::endl
                        << "Deformation field size/index: "
                        << fieldRegion.GetSize() << " / "
                        << fieldRegion.GetIndex() << std::endl );
      }
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

  /* Inititialize interpolators */
  m_FixedInterpolator->SetInputImage( m_FixedImage );
  m_MovingInterpolator->SetInputImage( m_MovingImage );

  /* Setup for image gradient calculations and possible BSpline interpolation.
    Check if the interpolator is of type BSplineInterpolateImageFunction.
    If so, we can make use of its EvaluateDerivatives method.
    Otherwise, we instantiate an external central difference
    derivative calculator.
  */
  /* Fixed */
  FixedBSplineInterpolatorType *fixedTestPtr =
        dynamic_cast< FixedBSplineInterpolatorType * >
          ( this->m_FixedInterpolator.GetPointer() );
  if ( !fixedTestPtr )
    {
    m_FixedInterpolatorIsBSpline = false;
    m_FixedBSplineInterpolator = NULL;
    if( !m_PrecomputeImageGradient )
      {
      m_FixedGradientCalculator = FixedGradientCalculatorType::New();
      m_FixedGradientCalculator->UseImageDirectionOn();
      m_FixedGradientCalculator->SetInputImage(this->m_FixedImage);
      }
    itkDebugMacro("Fixed Interpolator is not BSpline");
    }
  else
    {
    m_FixedInterpolatorIsBSpline = true;
    m_FixedBSplineInterpolator = fixedTestPtr;
    m_FixedBSplineInterpolator->SetNumberOfThreads(m_NumberOfThreads);
    m_FixedBSplineInterpolator->UseImageDirectionOn();
    m_FixedGradientCalculator = NULL;
    itkDebugMacro("Fixed Interpolator is BSpline");
    }
  /* Moving */
  MovingBSplineInterpolatorType *movingTestPtr =
        dynamic_cast< MovingBSplineInterpolatorType * >
          ( this->m_MovingInterpolator.GetPointer() );
  if ( !movingTestPtr )
    {
    m_MovingInterpolatorIsBSpline = false;
    m_MovingBSplineInterpolator = NULL;
    if( ! m_PrecomputeImageGradient )
      {
      m_MovingGradientCalculator = MovingGradientCalculatorType::New();
      m_MovingGradientCalculator->UseImageDirectionOn();
      m_MovingGradientCalculator->SetInputImage(this->m_MovingImage);
      }
    itkDebugMacro("Moving Interpolator is not BSpline");
    }
  else
    {
    m_MovingInterpolatorIsBSpline = true;
    m_MovingBSplineInterpolator = movingTestPtr;
    m_MovingBSplineInterpolator->SetNumberOfThreads(m_NumberOfThreads);
    m_MovingBSplineInterpolator->UseImageDirectionOn();
    m_MovingGradientCalculator = NULL;
    itkDebugMacro("Moving Interpolator is BSpline");
    }
  /* If user set to use pre-calculated gradient image, and both interpolators
   * aren't bspline, then we need to calculate the gradient images. */
  if ( m_PrecomputeImageGradient &&
       !( m_FixedInterpolatorIsBSpline && m_MovingInterpolatorIsBSpline ) )
    {
    //NOTE: This could be broken into separate fixed- and moving-gaussian
    // calculations if users will be having one interpolator as bspline
    // and not the other.
    ComputeGaussianGradient();
    }
}

/*
 * Initiate thread evaluation.
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

  /* Allocate and initialize memory for holding derivative results.*/
  this->m_DerivativesPerThread.resize( this->m_NumberOfThreads );
  /* This size always comes from the moving image */
  unsigned long globalDerivativeSize =
    this->m_MovingImageTransform->GetNumberOfParameters();
  itkDebugMacro("ImageToImageObjectMetric::Initialize: deriv size  "
                  << globalDerivativeSize << std::endl);
  /* NOTE: this does *not* get init'ed to 0 here. */
  if( this->m_DerivativeResult->GetSize() != globalDerivativeSize )
    {
    this->m_DerivativeResult->SetSize( globalDerivativeSize );
    }
  /* For transforms with local support, e.g. deformation field,
   * use a single derivative container that's updated by region
   * in multiple threads. */
  if ( this->m_MovingImageTransform->HasLocalSupport() )
    {
    itkDebugMacro(
      "ImageToImageObjectMetric::Initialize: tx HAS local support\n");
    for (ThreadIdType i=0; i<this->m_NumberOfThreads; i++)
      {
      this->m_DerivativesPerThread[i].SetData(
                                    this->m_DerivativeResult->data_block(),
                                    this->m_DerivativeResult->Size(),
                                    false );
      }
    }
  else
    {
    itkDebugMacro(
      "ImageToImageObjectMetric::Initialize: tx does NOT have local support\n");
    /* Global transforms get a separate derivatives container for each thread
     * that holds the result for a particular region. */
    for (ThreadIdType i=0; i<this->m_NumberOfThreads; i++)
      {
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
  // Initialize some threading values that require per-iteration
  // initialization.
  for (ThreadIdType i=0; i<this->m_NumberOfThreads; i++)
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
  /* Clear derivative final result. This will
   * require an option to skip for use with multivariate metric. */
  this->m_DerivativeResult->Fill( 0 );
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
  if ( ! this->m_MovingImageTransform->HasLocalSupport() )
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
    if ( ! this->m_MovingImageTransform->HasLocalSupport() )
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
  FixedImagePixelType         fixedImageValue;
  FixedImageDerivativesType   fixedImageDerivatives;
  MovingOutputPointType       mappedMovingPoint;
  MovingImagePixelType        movingImageValue;
  MovingImageDerivativesType  movingImageDerivatives;
  bool                        pointIsValid;
  MeasureType                 metricValueResult;
  MeasureType                 metricValueSum = 0;

  DerivativeType   localDerivativeResult( self->GetNumberOfLocalParameters() );

  /* Iterate over the sub region */
  ItV.GoToBegin();
  while( !ItV.IsAtEnd() )
  {
    /* Get the virtual point */
    self->m_VirtualDomainImage->TransformIndexToPhysicalPoint(
                                              ItV.GetIndex(), virtualPoint);

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
      }
    catch( ExceptionObject & exc )
      {
      //NOTE: there must be a cleaner way to do this:
      std::string msg("Caught exception: \n");
      msg += exc.what();
      ExceptionObject err(__FILE__, __LINE__, msg);
      throw err;
      }

    /* Call the user method in derived classes to do the specific
     * calculations for value and derivative. */
    try
      {
      if( pointIsValid )
        {
        pointIsValid = self->GetValueAndDerivativeProcessPoint(
                 virtualPoint,
                 mappedFixedPoint, fixedImageValue, fixedImageDerivatives,
                 mappedMovingPoint, movingImageValue, movingImageDerivatives,
                 metricValueResult, localDerivativeResult, threadID );
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
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::StoreDerivativeResult(  DerivativeType & derivative,
                           const VirtualIndexType & virtualIndex,
                           ThreadIdType threadID )
{
  DerivativeType & derivativeResult = this->m_DerivativesPerThread[threadID];
  if ( ! this->m_MovingImageTransform->HasLocalSupport() )
    {
    derivativeResult += derivative;
    }
  else
    {
    // Dense transform, e.g. deformation field.
    // update derivative at some index
    // this requires the moving image deformation field to be
    // same size as virtual image, and that VirtualImage PixelType
    // is scalar.
    try
      {
      OffsetValueType offset =
        this->m_VirtualDomainImage->ComputeOffset(virtualIndex);
      offset*=this->m_MovingImageTransform->GetNumberOfLocalParameters();
      for (unsigned int i=0;
            i < this->m_MovingImageTransform->GetNumberOfLocalParameters(); i++)
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
 * Transform a point from VirtualImage domain to FixedImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateFixedPoint(const VirtualPointType & point,
                      FixedImagePointType & mappedFixedPoint,
                      bool & pointIsValid,
                      FixedImagePixelType & fixedImageValue,
                      bool computeImageGradient,
                      FixedImageDerivativesType & fixedGradient,
                      ThreadIdType threadID) const
{
  pointIsValid = true;
  fixedImageValue = 0;

  // map the point into fixed space
  mappedFixedPoint = m_FixedImageTransform->TransformPoint( point );

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
      m_FixedImageTransform->TransformCovariantVector( fixedGradient,
                                                        mappedFixedPoint );
    }
}

/*
 * Transform a point from VirtualImage domain to MovingImage domain.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateMovingPoint(const VirtualPointType & point,
                      MovingImagePointType & mappedMovingPoint,
                      bool & pointIsValid,
                      MovingImagePixelType & movingImageValue,
                      bool computeImageGradient,
                      MovingImageDerivativesType & movingGradient,
                      ThreadIdType threadID) const
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

  if( ! pointIsValid )
    {
    return;
    }

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
    // Transform into the virtual space. See TransformAndEvaluateFixedPoint.
    movingGradient =
      m_MovingImageTransform->TransformCovariantVector( movingGradient,
                                                        mappedMovingPoint );
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
                              ThreadIdType threadID) const
{
  if ( m_FixedInterpolatorIsBSpline )
    {
    // Computed Fixed image gradient using derivative BSpline kernel.
    gradient = m_FixedBSplineInterpolator->EvaluateDerivative(mappedPoint,
                                                              threadID);
    }
  else
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
                              ThreadIdType threadID) const
{
  if ( m_MovingInterpolatorIsBSpline )
    {
    // Computed moving image gradient using derivative BSpline kernel.
    gradient = m_MovingBSplineInterpolator->EvaluateDerivative(mappedPoint,
                                                               threadID);
    }
  else
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
}

/*
 * Get the size of derivative result
 */

/* See .h comments...

template<class TFixedImage,class TMovingImage,class TVirtualImage>
SizeValueType
ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage >
::GetLocalDerivativeSize() const
{
  switch( this->m_DerivativeSource )
    {
    case Self::Moving:
      // For global transforms, this just returns GetNumberOfParameters.
      // For local/dense transforms, this returns number of local parameters.
      return this->m_MovingImageTransform->GetNumberOfLocalParameters();
      break;
    case Self::Fixed:
      return this->m_FixedImageTransform->GetNumberOfLocalParameters();
      break;
    case Self::Both:
      // TODO: this isn't completely settled, what to do here. Derived
      // classes may need to override to provide different behavior.
      return VirtualImageDimension;
      break;
    default:
      itkExceptionMacro("Invalid DerivativeSource enum: "
                         << this->m_DerivativeSource );
    }
  return 0;
}

*/
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
::UpdateTransformParameters( DerivativeType & derivative )
{
  if( derivative.GetSize() != this->GetNumberOfParameters() )
    {
    itkExceptionMacro("derivative is not the proper size");
    }
  this->m_MovingImageTransform->UpdateTransformParameters( derivative );
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

}//namespace itk

#endif
