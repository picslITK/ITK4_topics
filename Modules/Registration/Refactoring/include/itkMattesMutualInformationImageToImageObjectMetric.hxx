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

#ifndef __itkMattesMutualInformationImageToImageObjectMetric_txx
#define __itkMattesMutualInformationImageToImageObjectMetric__txx

#include "itkMattesMutualInformationImageToImageObjectMetric.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageIterator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "vnl/vnl_math.h"

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
  this->m_NumberOfHistogramBins=50;
  this->m_MovingImageNormalizedMin=0.0;
  this->m_FixedImageNormalizedMin=0.0;
  this->m_FixedImageTrueMin=(0.0);
  this->m_FixedImageTrueMax=(0.0);
  this->m_MovingImageTrueMin=(0.0);
  this->m_MovingImageTrueMax=(0.0);
  this->m_FixedImageBinSize=(0.0);
  this->m_MovingImageBinSize=(0.0);

  this->m_JointPDFSum=(0.0);
  this->m_ThreaderJointPDFInterpolator=NULL;
  this->m_ThreaderMovingImageMarginalPDFInterpolator=NULL;
  this->m_ThreaderFixedImageMarginalPDFInterpolator=NULL;
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
    delete[] m_ThreaderFixedImageMarginalPDFInterpolator;
    }
  if ( this->m_ThreaderMovingImageMarginalPDFInterpolator != NULL )
    {
    delete[] m_ThreaderMovingImageMarginalPDFInterpolator;
    }
  if ( this->m_ThreaderJointPDFInterpolator != NULL )
    {
    delete[] m_ThreaderJointPDFInterpolator;
    }

  /*
  if ( this->m_FixedImageMarginalPDF != NULL )
    {
    delete[] this->m_FixedImageMarginalPDF;
    }
  this->m_FixedImageMarginalPDF = NULL;

  if ( this->m_MovingImageMarginalPDF != NULL )
    {
    delete[] this->m_MovingImageMarginalPDF;
    }
  this->m_MovingImageMarginalPDF = NULL;

  if ( this->m_ThreaderJointPDF != NULL )
    {
    delete[] this->m_ThreaderJointPDF;
    }
  this->m_ThreaderJointPDF = NULL;

  if ( this->m_ThreaderJointPDFStartBin != NULL )
    {
    delete[] this->m_ThreaderJointPDFStartBin;
    }
  this->m_ThreaderJointPDFStartBin = NULL;

  if ( this->m_ThreaderJointPDFEndBin != NULL )
    {
    delete[] this->m_ThreaderJointPDFEndBin;
    }
  this->m_ThreaderJointPDFEndBin = NULL;

  if ( this->m_ThreaderJointPDFSum != NULL )
    {
    delete[] this->m_ThreaderJointPDFSum;
    }
  this->m_ThreaderJointPDFSum = NULL;
  */
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
  os << indent << "FixedImageNormalizedMin: ";
  os << this->m_FixedImageNormalizedMin << std::endl;
  os << indent << "MovingImageNormalizedMin: ";
  os << this->m_MovingImageNormalizedMin << std::endl;
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

  /**
   * Compute binsize for the histograms.
   *
   * The binsize for the image intensities needs to be adjusted so that
   * we can avoid dealing with boundary conditions using the cubic
   * spline as the Parzen window.  We do this by increasing the size
   * of the bins so that the joint histogram becomes "padded" at the
   * borders. Because we are changing the binsize,
   * we also need to shift the minimum by the padded amount in order to
   * avoid minimum values filling in our padded region.
   *
   * Note that there can still be non-zero bin values in the padded region,
   * it's just that these bins will never be a central bin for the Parzen
   * window.
   *
   */
  const int padding = 0;  // this will pad by 0 bins

  this->m_FixedImageBinSize = ( this->m_FixedImageTrueMax - this->m_FixedImageTrueMin )
                        / static_cast< double >( this->m_NumberOfHistogramBins
                                                 - 2 * padding );
  this->m_FixedImageNormalizedMin = this->m_FixedImageTrueMin / this->m_FixedImageBinSize
                              - static_cast< double >( padding );

  this->m_MovingImageBinSize = ( this->m_MovingImageTrueMax - this->m_MovingImageTrueMin )
                         / static_cast< double >( this->m_NumberOfHistogramBins
                                                  - 2 * padding );
  this->m_MovingImageNormalizedMin = this->m_MovingImageTrueMin / this->m_MovingImageBinSize
                               - static_cast< double >( padding );

  itkDebugMacro("FixedImageNormalizedMin: "  << this->m_FixedImageNormalizedMin);
  itkDebugMacro("MovingImageNormalizedMin: " << this->m_MovingImageNormalizedMin);
  itkDebugMacro("FixedImageBinSize: "        << this->m_FixedImageBinSize);
  itkDebugMacro("MovingImageBinSize; "       << this->m_MovingImageBinSize);

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
  typename JointPDFType::PointType origin;
  origin[0]=0;
  origin[1]=0;
  this->m_JointPDF->SetOrigin(origin);
  typename JointPDFType::SpacingType spacing;
  spacing[0]=1/(double)this->m_NumberOfHistogramBins;
  spacing[1]=1/(double)this->m_NumberOfHistogramBins;
  this->m_JointPDF->SetSpacing(spacing);
  this->m_JointPDF->Allocate();
  this->m_JointPDFBufferSize = jointPDFSize[0] * jointPDFSize[1] * sizeof( PDFValueType );

  // do the same thing for the marginal pdfs
  this->m_FixedImageMarginalPDF = MarginalPDFType::New();
  this->m_MovingImageMarginalPDF = MarginalPDFType::New();

  // Instantiate a region, index, size
  typedef typename MarginalPDFType::RegionType MarginalPDFRegionType;
  typedef typename MarginalPDFType::SizeType MarginalPDFSizeType;
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
  //  this->m_MarginalPDFBufferSize = marginalPDFSize[0] * sizeof( PDFValueType );

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
     this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID] = MarginalPDFInterpolatorType::New();
     this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->SetInputImage(this->m_FixedImageMarginalPDF);
     this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID] = MarginalPDFInterpolatorType::New();
     this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->SetInputImage(this->m_MovingImageMarginalPDF);
     }

  /*
   // Threaded data
   if (this->m_ThreaderJointPDF != NULL)
    {
    delete[] this->m_ThreaderJointPDF;
    }
   this->m_ThreaderJointPDF = new typename JointPDFType::Pointer[this->GetNumberOfThreads() ];

   for (ThreadIdType threadID = 0; threadID < this->GetNumberOfThreads(); threadID++)
     {
     this->m_ThreaderJointPDF[threadID] = JointPDFType::New();
     }
   std::cout << " Init 9 - Done " << std::endl;
  */

}

/** Compute fixed/moving Parzen Window Index
 *  using the fixed/moving parzen window term */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
int
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::ComputeParzenWindowIndex(double parzenWindowTerm )
{
  int parzenWindowIndex = static_cast<int> (parzenWindowTerm);
  unsigned int padding=0;
  // check that the extreme values are in valid bins
  if (parzenWindowIndex < padding)
    {
    parzenWindowIndex = 0;
    }
  else
    {
    const int nindex = static_cast<int> (this->m_NumberOfHistogramBins) - padding - 1;
    if (parzenWindowIndex > nindex)
      {
      parzenWindowIndex = nindex;
      }
    }
  return parzenWindowIndex;
}

/** Get the value and derivative */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::GetValueAndDerivative(MeasureType & value, DerivativeType & derivative)
{
  // Initialize the outputs to zero
  value = NumericTraits< MeasureType >::Zero;

  // Initialize the joint pdf and the fixed and moving image marginal pdfs
  PDFValueType pdfzero=NumericTraits< PDFValueType >::Zero;
  this->m_JointPDF->FillBuffer(pdfzero);
  this->m_FixedImageMarginalPDF->FillBuffer(pdfzero);
  this->m_MovingImageMarginalPDF->FillBuffer(pdfzero);

  /** First, we compute the joint histogram */
  /* Create an iterator over the virtual sub region */
  ImageRegionConstIteratorWithIndex<typename Superclass::VirtualImageType>
    ItV( this->GetVirtualDomainImage(), this->GetVirtualDomainImage()->GetRequestedRegion() );

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
                                              mappedFixedPoint,
                                              pointIsValid,
                                              fixedImageValue,
                                              false /*compute gradient*/,
                                              fixedImageGradients,
                                              0 );
        if( pointIsValid )
          {
          this->TransformAndEvaluateMovingPoint( ItV.GetIndex(),
                                                virtualPoint,
                                                mappedMovingPoint,
                                                pointIsValid,
                                                movingImageValue,
                                                false /*compute gradient*/,
                                                movingImageGradients,
                                                0 );
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
      typename JointPDFType::PointType jointPDFpoint;
      jointPDFpoint[0]=(fixedImageValue-this->m_FixedImageTrueMin)/(this->m_FixedImageTrueMax-this->m_FixedImageTrueMin);
      jointPDFpoint[1]=(movingImageValue-this->m_MovingImageTrueMin)/(this->m_MovingImageTrueMax-this->m_MovingImageTrueMin);
      JointPDFIndexType  jointPDFIndex;
      this->m_JointPDF->TransformPhysicalPointToIndex(jointPDFpoint,jointPDFIndex);
      this->m_JointPDF->SetPixel(jointPDFIndex,this->m_JointPDF->GetPixel(jointPDFIndex)+1);
      std::cout <<" jointPDFpoint " << jointPDFpoint << " ind " << jointPDFIndex << " val " << this->m_JointPDF->GetPixel(jointPDFIndex) << " nbins " << this->m_NumberOfHistogramBins << std::endl;
    //next index
    ++ItV;
  }
  bool smoothjh=true;
  if (smoothjh)
    {
      typedef DiscreteGaussianImageFilter<JointPDFType,JointPDFType> dgtype;
      typename dgtype::Pointer dg=dgtype::New();
      dg->SetInput(this->m_JointPDF);
      dg->SetVariance(1.);
      dg->SetUseImageSpacingOff();
      dg->SetMaximumError(.01f);
      dg->Update();
      this->m_JointPDF=dg->GetOutput();
    }

  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator ( this->m_JointPDF, this->m_JointPDF->GetRequestedRegion() );

  // Compute joint PDF normalization factor (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;
  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
      jointPDFSum += jointPDFIterator.Get();
      ++jointPDFIterator;
    }
  if ( fabs(jointPDFSum) < 1.e-9 )
    {
      std::cout <<" jointPDFSum " << jointPDFSum << std::endl;
      itkExceptionMacro( "Joint PDF summed to zero" );
    }

  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }

  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter( this->m_JointPDF, this->m_JointPDF->GetBufferedRegion() );

  linearIter.SetDirection( 0 );
  linearIter.GoToBegin();
  unsigned int fixedIndex = 0;
  while( !linearIter.IsAtEnd() )
    {
      double sum = 0.0;
      while( !linearIter.IsAtEndOfLine() )
     {
       sum += linearIter.Get();
       ++linearIter;
     }
      MarginalPDFIndexType mind;
      mind[0]=fixedIndex;
      m_FixedImageMarginalPDF->SetPixel(mind,static_cast<PDFValueType>(sum));
      linearIter.NextLine();
      ++fixedIndex;
    }

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;
  while( !linearIter.IsAtEnd() )
    {
      double sum = 0.0;
      while( !linearIter.IsAtEndOfLine() )
     {
       sum += linearIter.Get();
       ++linearIter;
     }
      MarginalPDFIndexType mind;
      mind[0]=movingIndex;
      m_MovingImageMarginalPDF->SetPixel(mind,static_cast<PDFValueType>(sum));
      linearIter.NextLine();
      ++movingIndex;
    }
  // Multithreaded initiate and process sample
  this->GetValueAndDerivativeMultiThreadedInitiate( derivative );

  // Post processing
  this->GetValueAndDerivativeMultiThreadedPostProcess( true /*doAverage*/ );

}


/** Process the sample point*/
template <class TFixedImage,class TMovingImage,class TVirtualImage>
bool
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           virtualPoint,
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageGradientsType &  fixedImageGradients,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageGradientsType & movingImageGradients,
                    MeasureType &                      metricValueReturn,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType                       threadID)
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

  /* Use a pre-allocated jacobian object for efficiency */
  typename Superclass::MovingTransformJacobianType & jacobian =
                            this->m_MovingTransformJacobianPerThread[threadID];

  /** For dense transforms, this returns identity */
  this->m_MovingTransform->GetJacobianWithRespectToParameters(
                                                            mappedMovingPoint,
                                                            jacobian);
  /** the scalingfactor is the MI specific scaling of the image gradient and jacobian terms */
  double scalingfactor=0;
  itk::ContinuousIndex< double, 2 > pdfind;
  pdfind[1]=(fixedImageValue-this->m_FixedImageTrueMin)/(this->m_FixedImageTrueMax-this->m_FixedImageTrueMin);
  pdfind[0]= (movingImageValue-this->m_MovingImageTrueMin)/(this->m_MovingImageTrueMax-this->m_MovingImageTrueMin);
  double jointPDFValue=this->m_ThreaderJointPDFInterpolator[threadID]->EvaluateAtContinuousIndex(pdfind);
  double dJPDF = 0; //this->m_ThreaderJointPDFInterpolator[threadID]->EvaluateDerivativeAtContinuousIndex( pdfind ))[1];

  itk::ContinuousIndex< double, 1 > mind;
  mind[0]=pdfind[0];
  double movingImagePDFValue = m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->EvaluateAtContinuousIndex(mind);
  double dMmPDF = m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->EvaluateDerivativeAtContinuousIndex(mind)[0];

  double term1=0,term2=0,eps=1.e-12; // FIXME need a 'small value' generalization
  if( jointPDFValue > eps &&  (movingImagePDFValue) > 0)
  {
    term1 = dJPDF/jointPDFValue;
    term2 = dMmPDF/movingImagePDFValue;
    scalingfactor =  (term1*(-1.0)+term2);
  }  // end if-block to check non-zero bin contribution
  else scalingfactor = 0;

  for ( unsigned int par = 0;
          par < this->GetNumberOfLocalParameters(); par++ )
  {
    double sum = 0.0;
    for ( unsigned int dim = 0; dim < this->MovingImageDimension; dim++ )
      {
        sum += scalingfactor * jacobian(dim, par) * movingImageGradients[dim];
      }
    localDerivativeReturn[par] = sum;
  }
  return true;
}


/** Post processing
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::GetValueAndDerivativeMultiThreadedPostProcess(bool doAverage)
{
  typedef itk::ImageIteratorWithIndex<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType it( this->m_JointPDF,this->m_JointPDF->GetRequestedRegion());
  // MY POST PROCESS
  double jointPDFSum = 0.0;
  for (unsigned int t = 0; t < this->GetNumberOfThreads() - 1; t++ )
    {
    it.GoToBegin();
    while(!it.IsAtEnd())
      {
     it.Set(it.Get()+this->m_ThreaderJointPDF[t]->GetPixel(it.GetIndex()));
     jointPDFSum+=it.Get();
     // FIXME also need to update the marginals and other variables
      }
    }
}
*/

/* Compute PDF Derivatives
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::ComputePDFDerivatives(int & fixedImageParzenWindowIndex,
                        int & movingImageParzenWindowIndex,
                        const double & cubicBSplineDerivativeValue,
                        const MovingImageGradientsType & movingImageGradients,
                        ThreadIdType threadID)
{
  // Update the bins for the current point
  JointPDFDerivativesValueType *derivPtr;

  derivPtr = this->m_JointPDFDerivatives->GetBufferPointer()
                 + ( fixedImageParzenWindowIndex  * this->m_JointPDFDerivatives->GetOffsetTable()[2] )
                 + ( movingImageParzenWindowIndex * this->m_JointPDFDerivatives->GetOffsetTable()[1] );

  if ( !this->m_TransformIsBSpline )
    {
    // Generic version which works for all transforms.


    // Compute the transform Jacobian. Should pre-compute
    typedef typename TransformType::JacobianType JacobianType;
    MovingTransformPointer movingTransform = this->m_MovingImageTransform;
    movingTransform->GetJacobianWithRespectToParameters(
      this->fixedImageParzenWindowIndex, jacobian);
    for ( unsigned int mu = 0; mu < this->m_NumberOfParameters; mu++ )
      {
      double innerProduct = 0.0;
      for ( unsigned int dim = 0; dim < Superclass::FixedImageDimension; dim++ )
        {
        innerProduct += jacobian[dim][mu] * movingImageGradients[dim];
        }

      // Equation 24 of Thevenaz and Unser.
      const double derivativeContribution = innerProduct * cubicBSplineDerivativeValue;

      if ( this->m_UseExplicitPDFDerivatives )
        {
        *( derivPtr ) -= derivativeContribution;
        ++derivPtr;
        }
      }
    }
}*/

// Get the value
template <class TFixedImage, class TMovingImage, class TVirtualImage>
typename MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::MeasureType
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::GetValue()
{
  MeasureType value = 0;
  return value;
}


} // end namespace itk

#endif
