/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsImageToImageMetric.h,v $
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
  
  // Initialize the B spline kernels
  this->m_CubicBSplineKernel=NULL;
  this->m_CubicBSplineDerivativeKernel=NULL;
  
  this->m_JointPDFSum=(0.0);
  
  // Threading variables
  this->m_ThreaderJointPDF=(NULL);
  this->m_ThreaderJointPDFStartBin=(NULL);
  this->m_ThreaderJointPDFEndBin=(NULL);
  this->m_ThreaderJointPDFSum=(NULL);
}

/**
 * Destructor
 */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::~MattesMutualInformationImageToImageObjectMetric()
{
  if ( m_FixedImageMarginalPDF != NULL )
    {
    delete[] m_FixedImageMarginalPDF;
    }
  m_FixedImageMarginalPDF = NULL;
  
  if ( m_MovingImageMarginalPDF != NULL )
    {
    delete[] m_MovingImageMarginalPDF;
    }
  m_MovingImageMarginalPDF = NULL;
  
  if ( m_ThreaderJointPDF != NULL )
    {
    delete[] m_ThreaderJointPDF;
    }
  m_ThreaderJointPDF = NULL;
  
  if ( m_ThreaderJointPDFStartBin != NULL )
    {
    delete[] m_ThreaderJointPDFStartBin;
    }
  m_ThreaderJointPDFStartBin = NULL;

  if ( m_ThreaderJointPDFEndBin != NULL )
    {
    delete[] m_ThreaderJointPDFEndBin;
    }
  m_ThreaderJointPDFEndBin = NULL;

  if ( m_ThreaderJointPDFSum != NULL )
    {
    delete[] m_ThreaderJointPDFSum;
    }
  m_ThreaderJointPDFSum = NULL;
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
  if( this->m_JointPDFDerivatives.IsNotNull() )
    {
    os << indent << "JointPDFDerivatives: ";
    os << this->m_JointPDFDerivatives;
    }  
}


/** Initialize the metric */
template <class TFixedImage,class TMovingImage,class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::Initialize() throw (itk::ExceptionObject)
{
  Superclass::Initialize();
  std::cout << " Init 1 " << std::endl;
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
                fi(this->m_FixedImage,this->m_FixedImage->GetBufferedRegion());
  
  std::cout << " Init 2 " << std::endl;
  while( !fi.IsAtEnd() )
    {
    typename TFixedImage::PointType fixedSpacePhysicalPoint;
    this->m_FixedImage->TransformIndexToPhysicalPoint
                                      (fi.GetIndex(), fixedSpacePhysicalPoint);
                                         
    
    if ( this->m_FixedImageMask.IsNull()  /* A null mask implies entire space is to be used.*/
         || this->m_FixedImageMask->IsInside(fixedSpacePhysicalPoint) )
       {
       const typename TFixedImage::PixelType currentValue = fi.Get();
       
       // update the Fixed Image true min accordingly
       if (this->m_FixedImageTrueMin < currentValue)
         {
         this->m_FixedImageTrueMin = currentValue;
         }
       // update the Fixed Image true max accordingly
       if (this->m_FixedImageTrueMax > currentValue)
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
              
  std::cout << " Init 3 " << std::endl;
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
       if (this->m_MovingImageTrueMin < currentValue)
         {
         this->m_MovingImageTrueMin = currentValue;
         }
       // update the Moving Image true max accordingly
       if (this->m_MovingImageTrueMax > currentValue)
         {
         this->m_MovingImageTrueMax = currentValue;
         }       
       }       
      ++mi;
    } 
     
  itkDebugMacro(" FixedImageMin: " << m_FixedImageTrueMin
                                   << " FixedImageMax: " 
                                   << m_FixedImageTrueMax << std::endl);
  itkDebugMacro(" MovingImageMin: " << m_MovingImageTrueMin
                                    << " MovingImageMax: " 
                                    << m_MovingImageTrueMax << std::endl);  
  std::cout << " Init 4 " << std::endl;
    
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
  const int padding = 2;  // this will pad by 2 bins

  m_FixedImageBinSize = ( m_FixedImageTrueMax - m_FixedImageTrueMin )
                        / static_cast< double >( m_NumberOfHistogramBins
                                                 - 2 * padding );
  m_FixedImageNormalizedMin = m_FixedImageTrueMin / m_FixedImageBinSize
                              - static_cast< double >( padding );

  m_MovingImageBinSize = ( m_MovingImageTrueMax - m_MovingImageTrueMin )
                         / static_cast< double >( m_NumberOfHistogramBins
                                                  - 2 * padding );
  m_MovingImageNormalizedMin = m_MovingImageTrueMin / m_MovingImageBinSize
                               - static_cast< double >( padding );

  std::cout << " Init 5 " << std::endl;
                   
  itkDebugMacro("FixedImageNormalizedMin: "  << m_FixedImageNormalizedMin);
  itkDebugMacro("MovingImageNormalizedMin: " << m_MovingImageNormalizedMin);
  itkDebugMacro("FixedImageBinSize: "        << m_FixedImageBinSize);
  itkDebugMacro("MovingImageBinSize; "       << m_MovingImageBinSize);
  
  /**
   * Allocate memory for the marginal PDF and initialize values
   * to zero. The marginal PDFs are stored as std::vector.
   */
  if ( m_FixedImageMarginalPDF != NULL )
    {
    delete[] m_FixedImageMarginalPDF;
    }
  m_FixedImageMarginalPDF = new PDFValueType[m_NumberOfHistogramBins];
  
  if ( m_MovingImageMarginalPDF != NULL )
    {
    delete[] m_MovingImageMarginalPDF;
    }
  m_MovingImageMarginalPDF = new PDFValueType[m_NumberOfHistogramBins];

  std::cout << " Init 6 " << std::endl;

  // Allocate memory for the joint PDF. 
  m_JointPDF = JointPDFType::New(); 
  
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
  m_JointPDF->SetRegions(jointPDFRegion);
  
  //By setting these values, the joint histogram physical locations will correspond to intensity values.
  typename JointPDFType::PointType origin;
  origin[0]=this->m_FixedImageTrueMin;
  origin[1]=this->m_MovingImageTrueMin;
  m_JointPDF->SetOrigin(origin);
  
  typename JointPDFType::SpacingType spacing;
  spacing[0]=this->m_FixedImageBinSize;
  spacing[1]=this->m_MovingImageBinSize;
  m_JointPDF->SetSpacing(spacing);
  
  m_JointPDF->Allocate();  
  m_JointPDFBufferSize = jointPDFSize[0] * jointPDFSize[1] * sizeof( PDFValueType );
  std::cout << " Init 7 " << std::endl;
   
   // Threaded data
   if (this->m_ThreaderJointPDF != NULL)
    {
    delete[] this->m_ThreaderJointPDF;
    }
   this->m_ThreaderJointPDF = new typename JointPDFType::Pointer[this->GetNumberOfThreads() ];
   
   std::cout << " Init 8 " <<  this->GetNumberOfThreads() << std::endl;
   for (ThreadIdType threadID = 0; threadID < this->GetNumberOfThreads(); threadID++)
     {
     m_ThreaderJointPDF[threadID] = JointPDFType::New();
     }
   std::cout << " Init 9 - Done " << std::endl;
}

/** Compute fixed/moving Parzen Window Index 
 *  using the fixed/moving parzen window term */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
int
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::ComputeParzenWindowIndex(double parzenWindowTerm )
{
  int parzenWindowIndex = static_cast<int> (parzenWindowTerm);
  
  // check that the extreme values are in valid bins
  if (parzenWindowIndex < 2)
    {
    parzenWindowIndex = 2;
    }
  else
    {
    const int nindex = static_cast<int> (this->m_NumberOfHistogramBins) - 3;
    if (parzenWindowIndex > nindex)
      {
      parzenWindowIndex = nindex;
      }
    }
  return parzenWindowIndex;  
}

/** Process the sample point*/
template <class TFixedImage,class TMovingImage,class TVirtualImage>
bool
MattesMutualInformationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
::GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           itkNotUsed (virtualPoint),
                    const FixedImagePointType &        itkNotUsed (mappedFixedPoint),
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageDerivativesType &  itkNotUsed(fixedImageDerivatives),
                    const MovingImagePointType &       itkNotUsed (mappedMovingPoint),
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageDerivativesType & movingImageDerivatives,
                    MeasureType &                      itkNotUsed(metricValueResult),
                    DerivativeType &                   itkNotUsed(localDerivativeReturn),
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
  
  // compute the fixed image parzen window index
  double fixedImageParzenWindowTerm = fixedImageValue/this->m_FixedImageBinSize 
                                             - this->m_FixedImageNormalizedMin;
  int fixedImageParzenWindowIndex 
                        = ComputeParzenWindowIndex( fixedImageParzenWindowTerm );
    
  // increment the value in the fixed marginal pdf for the corresponding index
  this->m_FixedImageMarginalPDF[fixedImageParzenWindowIndex]++;
  
  // compute the moving image parzen window index
  double movingImageParzenWindowTerm = movingImageValue/this->m_MovingImageBinSize
                                             - this->m_MovingImageNormalizedMin;
                                                                                          
  int movingImageParzenWindowIndex = ComputeParzenWindowIndex(movingImageParzenWindowTerm);
  
  // go to the affected bin in the joint pdf
  JointPDFValueType * pdfPtr;
  
  if ( threadID > 0 )
    {
    pdfPtr = m_ThreaderJointPDF[threadID - 1]->GetBufferPointer()
             + ( fixedImageParzenWindowIndex
                 * m_ThreaderJointPDF[threadID - 1]
                 ->GetOffsetTable()[1] );
    }
  
  else
    {
    pdfPtr = m_JointPDF->GetBufferPointer()
             + (fixedImageParzenWindowIndex 
                * m_JointPDF->GetOffsetTable()[1]);
    }
              
  pdfPtr+= movingImageParzenWindowIndex;  
  
  // set the value to be the evaluated B-spline for the moving index.
  double movingImageParzenWindowArg = 
          static_cast<double> (movingImageParzenWindowIndex) 
          - movingImageParzenWindowTerm;
          
  *(pdfPtr) = this->m_CubicBSplineKernel->Evaluate (movingImageParzenWindowArg);  
  
  /* 
  TODO: Do this derivative stuff later 
  
  const double cubicBSplineDerivativeValue = this->m_CubicBSplineDerivativeKernel->Evaluate(movingImageParzenWindowArg);
  
  this->ComputePDFDerivatives (fixedImageParzenWindowIndex, 
                               movingImageParzenWindowIndex, 
                               cubicBSplineDerivativeValue,
                               movingImageDerivatives,
                               threadID);  
                               */
  return true;
}


/** Post processing */ 
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
	/** FIXME also need to update the marginals and other variables */
      }
    }  
   
}

/** Get the value and derivative */
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::GetValueAndDerivative(MeasureType & value, DerivativeType & derivative)
{
  // Initialize the outputs to zero
  value = NumericTraits< MeasureType >::Zero;
  
  /* TODO: Do derivative initialization*/
  
  // Initialize the joint pdf and the fixed and moving image marginal pdfs
  PDFValueType pdfzero=NumericTraits< PDFValueType >::Zero;
  this->m_JointPDF->FillBuffer(pdfzero);
  memset(this->m_FixedImageMarginalPDF,0,m_NumberOfHistogramBins * sizeof( PDFValueType ) );  
  memset(this->m_MovingImageMarginalPDF,0,m_NumberOfHistogramBins * sizeof( PDFValueType ) );  
          
  // Multithreaded initiate and process sample
  this->GetValueAndDerivativeMultiThreadedInitiate( derivative );
  
  // Post processing
  this->GetValueAndDerivativeMultiThreadedPostProcess( true /*doAverage*/ );
  
  if ( m_JointPDFSum == 0.0 )
    {
    itkExceptionMacro("Joint PDF summed to zero");
    }

  double fixedPDFSum = 0.0;
  const double normalizationFactor = 1.0 / m_JointPDFSum;  // this is alpha in the references
 
  // Set the moving marginal pdf and normalize it
  JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer();
  for (unsigned int i = 0; i < m_NumberOfHistogramBins; i++ )
    {
    fixedPDFSum += m_FixedImageMarginalPDF[i];
    PDFValueType * movingMarginalPtr = m_MovingImageMarginalPDF;
    for (unsigned int j = 0; j < m_NumberOfHistogramBins; j++ )
      {
      *( pdfPtr ) *= normalizationFactor;
      *( movingMarginalPtr++ ) += *( pdfPtr++ );
      }
    }
  
  // Normalize the fixed image marginal PDF
  if ( fixedPDFSum == 0.0 )
    {
    itkExceptionMacro("Fixed image marginal PDF summed to zero");
    }
  for ( unsigned int bin = 0; bin < m_NumberOfHistogramBins; bin++ )
    {
    m_FixedImageMarginalPDF[bin] /= fixedPDFSum;
    }
  
  /** Double sum over the histogram and compute the metric. */
     
  // Setup pointer to point to the first bin
  JointPDFValueType *jointPDFPtr = m_JointPDF->GetBufferPointer();
  
  // Initialize sum to zero
  double sum = 0.0;
  
  for ( unsigned int fixedIndex = 0;
        fixedIndex < m_NumberOfHistogramBins;
        ++fixedIndex )
    {
    const double fixedImagePDFValue = m_FixedImageMarginalPDF[fixedIndex];
    
    for ( unsigned int movingIndex = 0;
          movingIndex < m_NumberOfHistogramBins;
          ++movingIndex, jointPDFPtr++ )
      {
      const double movingImagePDFValue = m_MovingImageMarginalPDF[movingIndex];
      const double jointPDFValue = *( jointPDFPtr );
      
      // check for non-zero bin contribution
      if ( jointPDFValue > 1e-16 &&  movingImagePDFValue > 1e-16 )
        {
        const double pRatio = vcl_log(jointPDFValue / movingImagePDFValue);

        if ( fixedImagePDFValue > 1e-16 )
          {
          sum += jointPDFValue * ( pRatio - vcl_log(fixedImagePDFValue) );
          }
        
        }
      
      } // inner loop    
    } // end outer loop
  
  // return the negative value
  value = static_cast<MeasureType> (-1.0*sum); 
}

/* Compute PDF Derivatives 
template <class TFixedImage, class TMovingImage, class TVirtualImage>
void
MattesMutualInformationImageToImageObjectMetric<TFixedImage,TMovingImage,TVirtualImage>
::ComputePDFDerivatives(int & fixedImageParzenWindowIndex, 
                        int & movingImageParzenWindowIndex, 
                        const double & cubicBSplineDerivativeValue,
                        const MovingImageDerivativesType & movingImageDerivatives,
                        ThreadIdType threadID)
{
  // Update the bins for the current point
  JointPDFDerivativesValueType *derivPtr;
  
  derivPtr = m_JointPDFDerivatives->GetBufferPointer()
                 + ( fixedImageParzenWindowIndex  * m_JointPDFDerivatives->GetOffsetTable()[2] )
                 + ( movingImageParzenWindowIndex * m_JointPDFDerivatives->GetOffsetTable()[1] );
  
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
        innerProduct += jacobian[dim][mu] * movingImageDerivatives[dim];
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
