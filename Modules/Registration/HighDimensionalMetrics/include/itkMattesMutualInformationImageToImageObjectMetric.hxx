/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMattesMutualInformationImageToImageObjectMetric.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkMattesMutualInformationImageToImageObjectMetric_hxx
#define __itkMattesMutualInformationImageToImageObjectMetric_hxx

#include "itkMattesMutualInformationImageToImageObjectMetric.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageIterator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "vnl/vnl_math.h"
#include "itkImageFileWriter.h"
namespace itk
{
/**
  Constructor
 */
template <class TFixedImage,class TMovingImage,class TVirtualImage>
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::MattesMutualInformationImageToImageObjectMetric()
{
  // Initialize the Marginal PDFs
  this->m_FixedImageMarginalPDF=NULL;
  this->m_MovingImageMarginalPDF=NULL;

  // Initialize the joint PDF to NULL
  this->m_JointPDF=NULL;

  // Initialize histogram properties
  this->m_NumberOfHistogramBins=32;
  this->m_FixedImageTrueMin=(0.0);
  this->m_FixedImageTrueMax=(0.0);
  this->m_MovingImageTrueMin=(0.0);
  this->m_MovingImageTrueMax=(0.0);
  this->m_FixedImageBinSize=(0.0);
  this->m_MovingImageBinSize=(0.0);
  this->m_Padding=2;
  this->m_JointPDFSum=(0.0);
  this->m_ThreaderJointPDFInterpolator=NULL;
  this->m_ThreaderMovingImageMarginalPDFInterpolator=NULL;
  this->m_ThreaderFixedImageMarginalPDFInterpolator=NULL;
  this->m_Log2=vcl_log(2.0);
  // Threading variables
  //  this->m_ThreaderJointPDF=(NULL);
  //  this->m_ThreaderJointPDFStartBin=(NULL);
  //  this->m_ThreaderJointPDFEndBin=(NULL);
  //  this->m_ThreaderJointPDFSum=(NULL);
}

/**
 * Destructor
 */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::~MattesMutualInformationImageToImageObjectMetric()
{
  if ( this->m_ThreaderFixedImageMarginalPDFInterpolator != NULL )
    {
    delete[] this->m_ThreaderFixedImageMarginalPDFInterpolator;
    }
  if ( this->m_ThreaderMovingImageMarginalPDFInterpolator != NULL )
    {
    delete[] this->m_ThreaderMovingImageMarginalPDFInterpolator;
    }
  if ( this->m_ThreaderJointPDFInterpolator != NULL )
    {
    delete[] this->m_ThreaderJointPDFInterpolator;
    }

}

/** Initialize the metric */
template <class TFixedImage,class TMovingImage,class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::Initialize() throw (itk::ExceptionObject)
{
  Superclass::Initialize();
  /** Get the fixed and moving image true max's and mins.
   *  Initialize them to the PixelType min and max. */
  this->m_FixedImageTrueMin =
                    vcl_numeric_limits<typename TFixedImage::PixelType>::max();
  this->m_FixedImageTrueMax =
                    vcl_numeric_limits<typename TFixedImage::PixelType>::min();
  this->m_MovingImageTrueMin =
                    vcl_numeric_limits<typename TMovingImage::PixelType>::max();
  this->m_MovingImageTrueMax =
                    vcl_numeric_limits<typename TMovingImage::PixelType>::min();

  /** Iterate through the fixed image and set the true
   *  max and min for the fixed image. */
  itk::ImageRegionConstIteratorWithIndex<TFixedImage>
                fi(this->m_FixedImage,this->m_FixedImage->GetRequestedRegion());

  while( !fi.IsAtEnd() )
    {
    typename TFixedImage::PointType fixedSpacePhysicalPoint;
    this->m_FixedImage->TransformIndexToPhysicalPoint(fi.GetIndex(), fixedSpacePhysicalPoint);
    if ( this->m_FixedImageMask.IsNull()  /* A null mask implies entire space is to be used.*/
         || this->m_FixedImageMask->IsInside(fixedSpacePhysicalPoint) )
       {
       const typename TFixedImage::PixelType currentValue = fi.Get();
       // update the Fixed Image true min accordingly
       if ( currentValue < this->m_FixedImageTrueMin )
         {
         this->m_FixedImageTrueMin = currentValue;
         }
       // update the Fixed Image true max accordingly
       if ( currentValue > this->m_FixedImageTrueMax )
         {
         this->m_FixedImageTrueMax = currentValue;
         }
       }
      ++fi;
    }
  /** Iterate through the moving image and set the true
   * max and min for the moving image. */
  itk::ImageRegionConstIteratorWithIndex<TMovingImage>
              mi(this->m_MovingImage,this->m_MovingImage->GetBufferedRegion());

  while( !mi.IsAtEnd() )
    {
    typename TMovingImage::PointType movingSpacePhysicalPoint;
    this->m_MovingImage->TransformIndexToPhysicalPoint
                                      (mi.GetIndex(), movingSpacePhysicalPoint);

    if ( this->m_MovingImageMask.IsNull() /* A null mask implies entire space is to be used.*/
         || this->m_MovingImageMask->IsInside(movingSpacePhysicalPoint) )
       {
       const typename TMovingImage::PixelType currentValue=mi.Get();
       // update the Moving Image true min accordingly
       if ( currentValue < this->m_MovingImageTrueMin )
         {
         this->m_MovingImageTrueMin = currentValue;
         }
       // update the Moving Image true max accordingly
       if ( currentValue > this->m_MovingImageTrueMax )
         {
         this->m_MovingImageTrueMax = currentValue;
         }
       }
      ++mi;
    }
  itkDebugMacro(" FixedImageMin: " << this->m_FixedImageTrueMin
                                   << " FixedImageMax: "
                                   << this->m_FixedImageTrueMax << std::endl);
  itkDebugMacro(" MovingImageMin: " << this->m_MovingImageTrueMin
                                    << " MovingImageMax: "
                                    << this->m_MovingImageTrueMax << std::endl);


  // Allocate memory for the joint PDF.
  this->m_JointPDF = JointPDFType::New();

  // Instantiate a region, index, size
  JointPDFRegionType jointPDFRegion;
  JointPDFIndexType  jointPDFIndex;
  JointPDFSizeType   jointPDFSize;

  // the jointPDF is of size NumberOfBins x NumberOfBins
  jointPDFSize.Fill(m_NumberOfHistogramBins);
  jointPDFIndex.Fill(0);
  jointPDFRegion.SetIndex(jointPDFIndex);
  jointPDFRegion.SetSize(jointPDFSize);

  // Set the regions and allocate
  this->m_JointPDF->SetRegions(jointPDFRegion);

  //By setting these values, the joint histogram physical locations will correspond to intensity values.
  JointPDFSpacingType spacing;
  spacing[0]=1/(InternalComputationValueType)(this->m_NumberOfHistogramBins-(InternalComputationValueType)this->m_Padding*2-1);
  spacing[1]=spacing[0];
  this->m_JointPDF->SetSpacing(spacing);
  this->m_JointPDFSpacing=this->m_JointPDF->GetSpacing();
  JointPDFPointType origin;
  origin[0]=this->m_JointPDFSpacing[0]*(InternalComputationValueType)this->m_Padding*(-1.0);
  origin[1]=origin[0];
  this->m_JointPDF->SetOrigin(origin);
  this->m_JointPDF->Allocate();

  // do the same thing for the marginal pdfs
  this->m_FixedImageMarginalPDF = MarginalPDFType::New();
  this->m_MovingImageMarginalPDF = MarginalPDFType::New();

  // Instantiate a region, index, size
  typedef typename MarginalPDFType::RegionType  MarginalPDFRegionType;
  typedef typename MarginalPDFType::SizeType    MarginalPDFSizeType;
  MarginalPDFRegionType marginalPDFRegion;
  MarginalPDFIndexType  marginalPDFIndex;
  MarginalPDFSizeType   marginalPDFSize;

  // the marginalPDF is of size NumberOfBins x NumberOfBins
  marginalPDFSize.Fill(m_NumberOfHistogramBins);
  marginalPDFIndex.Fill(0);
  marginalPDFRegion.SetIndex(marginalPDFIndex);
  marginalPDFRegion.SetSize(marginalPDFSize);

  // Set the regions and allocate
  this->m_FixedImageMarginalPDF->SetRegions(marginalPDFRegion);
  this->m_MovingImageMarginalPDF->SetRegions(marginalPDFRegion);

  //By setting these values, the marginal histogram physical locations will correspond to intensity values.
  typename MarginalPDFType::PointType fixedorigin;
  typename MarginalPDFType::PointType movingorigin;
  fixedorigin[0]=origin[0];
  movingorigin[0]=origin[1];
  this->m_FixedImageMarginalPDF->SetOrigin(fixedorigin);
  this->m_MovingImageMarginalPDF->SetOrigin(movingorigin);
  typename MarginalPDFType::SpacingType mspacing;
  mspacing[0]=spacing[0];
  this->m_FixedImageMarginalPDF->SetSpacing(mspacing);
  mspacing[0]=spacing[1];
  this->m_MovingImageMarginalPDF->SetSpacing(mspacing);
  this->m_FixedImageMarginalPDF->Allocate();
  this->m_MovingImageMarginalPDF->Allocate();

   // Threaded data
   if (this->m_ThreaderJointPDFInterpolator != NULL)
    {
    delete[] this->m_ThreaderJointPDFInterpolator;
    }
   this->m_ThreaderJointPDFInterpolator = new JointPDFInterpolatorPointer[this->GetNumberOfThreads() ];
   if (this->m_ThreaderFixedImageMarginalPDFInterpolator != NULL)
    {
    delete[] this->m_ThreaderFixedImageMarginalPDFInterpolator;
    }
   this->m_ThreaderFixedImageMarginalPDFInterpolator = new MarginalPDFInterpolatorPointer[this->GetNumberOfThreads() ];

   if (this->m_ThreaderMovingImageMarginalPDFInterpolator != NULL)
    {
    delete[] this->m_ThreaderMovingImageMarginalPDFInterpolator;
    }
   this->m_ThreaderMovingImageMarginalPDFInterpolator = new MarginalPDFInterpolatorPointer[this->GetNumberOfThreads() ];

   for (ThreadIdType threadID = 0; threadID < this->GetNumberOfThreads(); threadID++)
     {
     this->m_ThreaderJointPDFInterpolator[threadID] = JointPDFInterpolatorType::New();
     this->m_ThreaderJointPDFInterpolator[threadID]->SetInputImage(this->m_JointPDF);
     //     this->m_ThreaderJointPDFInterpolator[threadID]->SetSplineOrder(3);
     this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID] = MarginalPDFInterpolatorType::New();
     this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->SetInputImage(this->m_FixedImageMarginalPDF);
     //     this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->SetSplineOrder(3);
     this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID] = MarginalPDFInterpolatorType::New();
     this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->SetInputImage(this->m_MovingImageMarginalPDF);
     //     this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->SetSplineOrder(3);
     }

}


/**  */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::EnforceJointHistogramBoundaryConditions()
{
  typename JointPDFType::IndexType jind;
  jind.Fill(0);
  for( SizeValueType i = 0; i < this->m_NumberOfHistogramBins; i++)
  {
    // first column
    jind[0] = i;
    jind[1] = this->m_Padding;
    PDFValueType val = this->m_JointPDF->GetPixel(jind);
    for (  int j = (int)this->m_Padding-1; j >= 0; j-- ) {
      jind[1] = j;
      this->m_JointPDF->SetPixel(jind,val);
    }
    // first row
    jind[1] = i;
    jind[0] = this->m_Padding;
    val = this->m_JointPDF->GetPixel(jind);
    for (  int j = (int)this->m_Padding-1; j >= 0; j-- ) {
      jind[0] = j;
      this->m_JointPDF->SetPixel(jind,val);
    }
    // last column
    jind[0] = i;
    jind[1] = this->m_NumberOfHistogramBins-this->m_Padding-1;
    val = this->m_JointPDF->GetPixel(jind);
    for ( unsigned int j = 1; j <= this->m_Padding; j++  ) {
      jind[1] = this->m_NumberOfHistogramBins-this->m_Padding-1+j;
      this->m_JointPDF->SetPixel(jind,val);
    }
    // last row
    jind[1] = i;
    jind[0] = this->m_NumberOfHistogramBins-this->m_Padding-1;
    val = this->m_JointPDF->GetPixel(jind);
    for ( unsigned int j = 1; j <= this->m_Padding; j++  ) {
      jind[0] = this->m_NumberOfHistogramBins-this->m_Padding-1+j;
      this->m_JointPDF->SetPixel(jind,val);
    }
  }

}

/** Prepare histograms for use in GetValueAndDerivative */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::UpdateHistograms()
{
  // Initialize the joint pdf and the fixed and moving image marginal pdfs
  PDFValueType pdfzero = NumericTraits< PDFValueType >::Zero;
  this->m_JointPDF->FillBuffer(pdfzero);
  this->m_FixedImageMarginalPDF->FillBuffer(pdfzero);
  this->m_MovingImageMarginalPDF->FillBuffer(pdfzero);

  /**
   * First, we compute the joint histogram
   */

  /* Create an iterator over the virtual sub region */
  ImageRegionConstIteratorWithIndex<typename Superclass::VirtualImageType>
    ItV( this->GetVirtualDomainImage(), this->GetVirtualDomainRegion() );

  typename Superclass::VirtualPointType            virtualPoint;
  typename Superclass::FixedOutputPointType        mappedFixedPoint;
  typename Superclass::FixedImagePixelType         fixedImageValue;
  FixedImageGradientsType   fixedImageGradients;
  typename Superclass::MovingOutputPointType       mappedMovingPoint;
  typename Superclass::MovingImagePixelType        movingImageValue;
  MovingImageGradientsType  movingImageGradients;
  bool                        pointIsValid = false;

  /* Iterate over the sub region */
  ItV.GoToBegin();
  while( !ItV.IsAtEnd() )
  {
    /* Get the virtual point */
    this->GetVirtualDomainImage()->TransformIndexToPhysicalPoint(
                                              ItV.GetIndex(), virtualPoint);
      try
        {
        this->TransformAndEvaluateFixedPoint( ItV.GetIndex(),
                                              virtualPoint,
                                              false /*compute gradient*/,
                                              mappedFixedPoint,
                                              fixedImageValue,
                                              fixedImageGradients,
                                              pointIsValid );
        if( pointIsValid )
          {
          this->TransformAndEvaluateMovingPoint( ItV.GetIndex(),
                                                virtualPoint,
                                                false /*compute gradient*/,
                                                mappedMovingPoint,
                                                movingImageValue,
                                                movingImageGradients,
                                                pointIsValid );
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
      /** add the paired intensity points to the joint histogram */
      JointPDFPointType jointPDFpoint;
      this->ComputeJointPDFPoint(fixedImageValue,movingImageValue, jointPDFpoint,0);
      JointPDFIndexType  jointPDFIndex;
      jointPDFIndex.Fill(0);
      this->m_JointPDF->TransformPhysicalPointToIndex(jointPDFpoint,jointPDFIndex);
      this->m_JointPDF->SetPixel(jointPDFIndex,this->m_JointPDF->GetPixel(jointPDFIndex)+1);
    //next index
    ++ItV;
  }
  /** here we copy the edges of the joint histogram to the padded region */
  //  this->EnforceJointHistogramBoundaryConditions();
  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator ( m_JointPDF, m_JointPDF->GetBufferedRegion() );

  // Compute joint PDF normalization factor (to ensure joint PDF sum adds to 1.0)
  InternalComputationValueType jointPDFSum = 0.0;
  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
      float temp = jointPDFIterator.Get();
      jointPDFSum += temp;
      ++jointPDFIterator;
    }

// of derivatives
  if ( jointPDFSum == 0.0 )
    {
    itkExceptionMacro( "Joint PDF summed to zero" );
    }


  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }

  bool smoothjh = true;
  if (smoothjh)
    {
      typedef DiscreteGaussianImageFilter<JointPDFType,JointPDFType> DgType;
      typename DgType::Pointer dg = DgType::New();
      dg->SetInput(this->m_JointPDF);
      dg->SetVariance(1.5);
      dg->SetUseImageSpacingOff();
      dg->SetMaximumError(.01f);
      dg->Update();
      this->m_JointPDF = dg->GetOutput();
    }


  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter(m_JointPDF, m_JointPDF->GetBufferedRegion() );
  linearIter.SetDirection( 0 );
  linearIter.GoToBegin();
  unsigned int fixedIndex = 0;
  while( !linearIter.IsAtEnd() )
    {
      InternalComputationValueType sum = 0.0;
      while( !linearIter.IsAtEndOfLine() )
      {
        sum += linearIter.Get();
        ++linearIter;
      }
      MarginalPDFIndexType mind;
      mind[0] = fixedIndex;
      m_FixedImageMarginalPDF->SetPixel(mind,static_cast<PDFValueType>(sum));
      linearIter.NextLine();
      ++fixedIndex;
    }

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;
  while( !linearIter.IsAtEnd() )
    {
      InternalComputationValueType sum = 0.0;
      while( !linearIter.IsAtEndOfLine() )
      {
        sum += linearIter.Get();
        ++linearIter;
      }
      MarginalPDFIndexType mind;
      mind[0] = movingIndex;
      m_MovingImageMarginalPDF->SetPixel(mind,static_cast<PDFValueType>(sum));
      linearIter.NextLine();
      ++movingIndex;
    }


  /** optionally renormalize the marginal
  for (unsigned int i=0; i<this->m_NumberOfHistogramBins; i++) {
    mind[0]=i;
    this->m_MovingImageMarginalPDF->SetPixel(  mind,this->m_MovingImageMarginalPDF->GetPixel(mind)/margsum );
  }
  */

  /** get the gradient images */
  /*  InternalComputationValueType sig=3*this->m_JointPDFSpacing[0];
  typename JPDFGradientFilterType::Pointer filter = JPDFGradientFilterType::New();
  filter->SetInput( this->m_JointPDF );
  filter->SetSigma( sig );
  filter->Update();
  this->m_JPDFGradientImage = filter->GetOutput();
  typename MarginalGradientFilterType::Pointer mfilter = MarginalGradientFilterType::New();
  mfilter->SetInput( this->m_MovingImageMarginalPDF );
  mfilter->SetSigma( sig );
  mfilter->Update();
  this->m_MarginalGradientImage = mfilter->GetOutput();
  */

  /* debug code for writing out jh over iterations
  static int iter=0;
  typedef itk::ImageFileWriter<JointPDFType> writertype;
  std::string s;
  std::stringstream out;
  out << iter;
  s = out.str();
  std::string outname=std::string("jh"+s+".nii.gz");
  typename writertype::Pointer writer = writertype::New();
  writer->SetFileName(outname.c_str());
  writer->SetInput( this->m_JointPDF );
  writer->Write();
  iter++;
  */
}

/** Get the value and derivative */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::GetValueAndDerivative(MeasureType & value, DerivativeType & derivative)
{
  // Prepare the histograms
  this->UpdateHistograms();

  // Calculate value
  this->m_Value = this->GetValue();
  std::cout <<" Mutual information value " << this->m_Value << std::endl;

  // Multithreaded initiate and process sample.
  // This will put results in 'derivative'.
  this->GetValueAndDerivativeMultiThreadedInitiate( derivative );

  // Post processing
  this->GetValueAndDerivativeMultiThreadedPostProcess( true /*doAverage*/ );

  // Return value.
  value = this->m_Value;
}


/** Process the sample point*/
template <class TFixedImage,class TMovingImage,class TVirtualImage>
bool
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &,
                    const FixedImagePointType &,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageGradientsType &,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageGradientsType &   movingImageGradients,
                    MeasureType &,
                    DerivativeType &                   localDerivativeReturn,
                    const ThreadIdType                 threadID)
{
  // check that the moving image sample is within the range of the true min
  // and max, hence being within the moving image mask
  if ( movingImageValue < this->m_MovingImageTrueMin )
    {
    return false;
    }
  else if ( movingImageValue > this->m_MovingImageTrueMax )
    {
    return false;
    }
  /** the scalingfactor is the MI specific scaling of the image gradient and jacobian terms */
  InternalComputationValueType scalingfactor = 0; // for scaling the jacobian terms

  JointPDFPointType jointPDFpoint;
  bool pointok = this->ComputeJointPDFPoint(
                      fixedImageValue,movingImageValue, jointPDFpoint,threadID);
  if ( !pointok )
    {
    return false;
    }
  InternalComputationValueType jointPDFValue =
    this->m_ThreaderJointPDFInterpolator[threadID]->Evaluate(jointPDFpoint);
  SizeValueType ind = 1;
  InternalComputationValueType dJPDF =
    this->ComputeJointPDFDerivative( jointPDFpoint, threadID , ind );
  typename MarginalPDFType::PointType mind;
  mind[0] = jointPDFpoint[ind];
  InternalComputationValueType movingImagePDFValue =
    this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->Evaluate(mind);
  InternalComputationValueType dMmPDF =
    this->ComputeMovingImageMarginalPDFDerivative( mind , threadID );

  InternalComputationValueType term1 = 0;
  InternalComputationValueType term2 = 0;
  InternalComputationValueType eps = 1.e-16;
  if( jointPDFValue > eps &&  (movingImagePDFValue) > eps )
    {
    const InternalComputationValueType pRatio =
                            vcl_log(jointPDFValue)-vcl_log(movingImagePDFValue);
    term1 = dJPDF*pRatio;
    term2 = vcl_log(2.0) * dMmPDF * jointPDFValue / movingImagePDFValue;
    scalingfactor =  ( term2 - term1 );
    }  // end if-block to check non-zero bin contribution
  else
    {
    scalingfactor = 0;
    }

  /* Use a pre-allocated jacobian object for efficiency */
  typename Superclass::JacobianType & jacobian =
                            this->m_MovingTransformJacobianPerThread[threadID];

  /** For dense transforms, this returns identity */
  this->m_MovingTransform->ComputeJacobianWithRespectToParameters(
                                                            mappedMovingPoint,
                                                            jacobian);

  // this correction is necessary for consistent derivatives across N threads
  typedef typename DerivativeType::ValueType    DerivativeValueType;
  DerivativeValueType floatingpointcorrectionresolution = 10000.0;
  // NOTE: change 'unsigned int' here when we have NumberOfParametersType
  // defined in metric base.
  for ( unsigned int par = 0; par < this->GetNumberOfLocalParameters(); par++ )
    {
    InternalComputationValueType sum = 0.0;
    for ( SizeValueType dim = 0; dim < this->MovingImageDimension; dim++ )
      {
      sum += scalingfactor * jacobian(dim, par) * movingImageGradients[dim];
      }
    localDerivativeReturn[par] = sum;
    intmax_t test = static_cast<intmax_t>
             ( localDerivativeReturn[par] * floatingpointcorrectionresolution );
    localDerivativeReturn[par] = static_cast<DerivativeValueType>
                                   ( test / floatingpointcorrectionresolution );
    }
  return true;
}


// Get the value
template <class TFixedImage, class TMovingImage, class TVirtualImage>
typename MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>::MeasureType
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::GetValue()
{
  /**
  1- The padding is 2 in this implementation.
  2- The MI energy is bounded in the range of [0  min(H(x),H(y))].
  3- The ComputeMutualInformation() iterator range should cover the entire PDF.
  4- The normalization is done based on NumberOfHistogramBins-1 instead of NumberOfHistogramBins. */
  InternalComputationValueType px,py,pxy;
  InternalComputationValueType total_mi = 0;
  InternalComputationValueType local_mi;
  InternalComputationValueType eps =
                        NumericTraits<InternalComputationValueType>::epsilon();
  typename JointPDFType::IndexType index;
  for (SizeValueType ii = 0; ii<m_NumberOfHistogramBins; ii++)
  {
    MarginalPDFIndexType mind;
    mind[0] = ii;
    px = this->m_FixedImageMarginalPDF->GetPixel(mind);
    for (SizeValueType jj = 0; jj<m_NumberOfHistogramBins; jj++)
    {
      mind[0] = jj;
      py = this->m_MovingImageMarginalPDF->GetPixel(mind);
      InternalComputationValueType denom = px * py;
      index[0] = ii;
      index[1] = jj;
      pxy = m_JointPDF->GetPixel(index);
      local_mi = 0;
      if ( fabs(denom) > eps )
      {
        if (pxy / denom > eps )
        {
          //the classic mi calculation
          local_mi = pxy * vcl_log( pxy / denom );
        }
      }
      total_mi += local_mi;
    } // over jh bins 2
  } // over jh bins 1
  return ( -1.0 * total_mi / this->m_Log2  );

}

template <class TFixedImage, class TMovingImage, class TVirtualImage>
bool
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::ComputeJointPDFPoint( const FixedImagePixelType fixedImageValue,
                        const MovingImagePixelType movingImageValue,
                        JointPDFPointType& jointPDFpoint,
                        const ThreadIdType threadID )
{
    InternalComputationValueType a =
        ( fixedImageValue - this->m_FixedImageTrueMin ) /
          ( this->m_FixedImageTrueMax - this->m_FixedImageTrueMin );
    InternalComputationValueType b =
        ( movingImageValue - this->m_MovingImageTrueMin ) /
           ( this->m_MovingImageTrueMax - this->m_MovingImageTrueMin );
    jointPDFpoint[0] = a;
    jointPDFpoint[1] = b;
    bool isInsideBuffer = this->m_ThreaderJointPDFInterpolator[threadID]->
                                                IsInsideBuffer(jointPDFpoint );
    return isInsideBuffer;
}

template <class TFixedImage, class TMovingImage, class TVirtualImage>
typename MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>::InternalComputationValueType
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::ComputeFixedImageMarginalPDFDerivative(
                                        const MarginalPDFPointType margPDFpoint,
                                        const ThreadIdType threadID )
{
  InternalComputationValueType offset = 0.5*this->m_JointPDFSpacing[0];
  InternalComputationValueType eps = this->m_JointPDFSpacing[0];
  MarginalPDFPointType         leftpoint = margPDFpoint;
  leftpoint[0] -= offset;
  MarginalPDFPointType  rightpoint = margPDFpoint;
  rightpoint[0] += offset;
  if (leftpoint[0] < eps ) leftpoint[0] = eps;
  if (rightpoint[0] < eps ) rightpoint[0] = eps;
  if (leftpoint[0] > 1 ) leftpoint[0] = 1;
  if (rightpoint[0] > 1  ) rightpoint[0] = 1;
  InternalComputationValueType delta = rightpoint[0]-leftpoint[0];
  if ( delta > 0 )
    {
    InternalComputationValueType deriv = this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->Evaluate(rightpoint) -
      this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->Evaluate(leftpoint);
    return deriv/delta;
    }
  else
    {
    return 0;
    }
}

template <class TFixedImage, class TMovingImage, class TVirtualImage>
typename MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>::InternalComputationValueType
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::ComputeMovingImageMarginalPDFDerivative(
                                        const MarginalPDFPointType margPDFpoint,
                                        const ThreadIdType threadID )
{
  InternalComputationValueType offset = 0.5*this->m_JointPDFSpacing[0];
  InternalComputationValueType eps = this->m_JointPDFSpacing[0];
  MarginalPDFPointType  leftpoint = margPDFpoint;
  leftpoint[0] -= offset;
  MarginalPDFPointType  rightpoint = margPDFpoint;
  rightpoint[0] += offset;
  if( leftpoint[0] < eps )
    {
    leftpoint[0] = eps;
    }
  if( rightpoint[0] < eps )
    {
    rightpoint[0] = eps;
    }
  if( leftpoint[0] > 1 )
    {
    leftpoint[0] = 1;
    }
  if( rightpoint[0] > 1  )
    {
    rightpoint[0] = 1;
    }
  InternalComputationValueType delta = rightpoint[0] - leftpoint[0];
  if ( delta > 0 )
    {
    InternalComputationValueType deriv =
      this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->Evaluate(rightpoint) -
      this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->Evaluate(leftpoint);
    return deriv/delta;
    }
  else
    {
    return 0;
    }
}

template <class TFixedImage, class TMovingImage, class TVirtualImage>
typename MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>::InternalComputationValueType
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::ComputeJointPDFDerivative( const JointPDFPointType jointPDFpoint,
                             const ThreadIdType threadID,
                             const SizeValueType ind )
{
  InternalComputationValueType offset = 0.5*this->m_JointPDFSpacing[ind];
  InternalComputationValueType eps = this->m_JointPDFSpacing[ind];
  JointPDFPointType  leftpoint = jointPDFpoint;
  leftpoint[ind] -= offset;
  JointPDFPointType  rightpoint = jointPDFpoint;
  rightpoint[ind] += offset;
  if (leftpoint[ind] < eps ) leftpoint[ind] = eps;
  if (rightpoint[ind] < eps ) rightpoint[ind] = eps;
  if (leftpoint[ind] > 1 ) leftpoint[ind] = 1;
  if (rightpoint[ind] > 1 ) rightpoint[ind] = 1;
  InternalComputationValueType delta = rightpoint[ind] - leftpoint[ind];
  InternalComputationValueType deriv = 0;
  if ( delta > 0 )
    {
    deriv = this->m_ThreaderJointPDFInterpolator[threadID]->Evaluate(rightpoint)-
          this->m_ThreaderJointPDFInterpolator[threadID]->Evaluate(leftpoint);
    return deriv/delta;
    }
  else
    {
    return deriv;
    }
}

/** Print function */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::PrintSelf (std::ostream & os, Indent indent) const
{
  // Print the superclass
  Superclass::Print(os);

  os << indent << "NumberOfHistogramBins: ";
  os << this->m_NumberOfHistogramBins << std::endl;
  os << indent << "MovingImageTrueMin: ";
  os << this->m_MovingImageTrueMin << std::endl;
  os << indent << "MovingImageTrueMax: ";
  os << this->m_MovingImageTrueMax << std::endl;
  os << indent << "FixedImageBinSize: ";
  os << this->m_FixedImageBinSize << std::endl;
  os << indent << "MovingImageBinSize: ";
  os << this->m_MovingImageBinSize << std::endl;

  if( this->m_JointPDF.IsNotNull() )
    {
    os << indent << "JointPDF: ";
    os << this->m_JointPDF << std::endl;
    }
}

} // end namespace itk

#endif
