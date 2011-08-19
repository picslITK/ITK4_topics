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

#ifndef __itkMattesMutualInformationImageToImageObjectMetric__h
#define __itkMattesMutualInformationImageToImageObjectMetric__h

#include "itkImageToImageObjectMetric.h"
#include "itkImage.h"
#include "itkBSplineDerivativeKernelFunction.h"

namespace itk
{
/** \class MattesMutualInformationImageToImageMetric
 * \brief Computes the mutual information between two images to be
 * registered using the method of Mattes et al.
 *
 * References:
 * [1] "Nonrigid multimodality image registration"
 *      D. Mattes, D. R. Haynor, H. Vesselle, T. Lewellen and W. Eubank
 *      Medical Imaging 2001: Image Processing, 2001, pp. 1609-1620.
 * [2] "PET-CT Image Registration in the Chest Using Free-form Deformations"
 *      D. Mattes, D. R. Haynor, H. Vesselle, T. Lewellen and W. Eubank
 *      IEEE Transactions in Medical Imaging. Vol.22, No.1,
        January 2003. pp.120-128.
 * [3] "Optimization of Mutual Information for MultiResolution Image
 *      Registration"
 *      P. Thevenaz and M. Unser
 *      IEEE Transactions in Image Processing, 9(12) December 2000.
 *
 * \ingroup ITK-RegistrationRefactoring
 */

template<class TFixedImage,class TMovingImage,class TVirtualImage = TFixedImage>
class ITK_EXPORT MattesMutualInformationImageToImageObjectMetric :
public ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
{
public:

  /** Standard class typedefs. */
  typedef MattesMutualInformationImageToImageObjectMetric       Self;
  typedef ImageToImageObjectMetric<TFixedImage, TMovingImage>   Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MattesMutualInformationImageToImageObjectMetric, ImageToImageObjectMetric);

  /** Type used for representing parameter values  */
  typedef typename Superclass::CoordinateRepresentationType
                                                  CoordinateRepresentationType;
  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  /** Superclass typedefs */
  typedef typename Superclass::MeasureType                 MeasureType;
  typedef typename Superclass::DerivativeType              DerivativeType;
  typedef typename Superclass::VirtualPointType            VirtualPointType;
  typedef typename Superclass::FixedImagePointType         FixedImagePointType;
  typedef typename Superclass::FixedImagePixelType         FixedImagePixelType;
  typedef typename Superclass::FixedImageDerivativesType   FixedImageDerivativesType;
  typedef typename Superclass::MovingImagePointType        MovingImagePointType;
  typedef typename Superclass::MovingImagePixelType        MovingImagePixelType;
  typedef typename Superclass::MovingImageDerivativesType  MovingImageDerivativesType;
  /** Value type of the PDF */
  typedef double PDFValueType;
  /** Typedef for the joint PDF and marginal PDF are stored as ITK Images. */
  typedef Image<PDFValueType,1> MarginalPDFType;
  typedef typename MarginalPDFType::IndexType MarginalPDFIndexType;
  typedef Image< PDFValueType, 2 >            JointPDFType;
  itkGetConstReferenceMacro(JointPDF,typename JointPDFType::Pointer);

  itkSetClampMacro( NumberOfHistogramBins, SizeValueType,
                    5, NumericTraits< SizeValueType >::max() );
  itkGetConstReferenceMacro(NumberOfHistogramBins, SizeValueType);

  /** Initialize the metric. Make sure all essential inputs are plugged in. */
  virtual void Initialize() throw (itk::ExceptionObject);

  /** Get both the value and derivative intializes the processing.
   *  For Mattes MI, we just compute the joint histogram / pdf here.
   *  This implementation single-threads the JH computation but it
   *  could be multi-threaded in the future.
   */
  void GetValueAndDerivative(MeasureType & value, DerivativeType & derivative);

  /** Called by GetValueAndDerivative */
  void ComputeJointPDFandMarginalPDFs();

  /** Compute parzen window index. */
  int ComputeParzenWindowIndex(double);

  /** Compute the PDF Derivatives. */
  void ComputePDFDerivatives(int pdfMovingIndex,
                             const MovingImageDerivativesType & movingImageGradientValue,
                             double & cubicBSplineDerivativeValue,
                             ThreadIdType threadID) const;
  /** Get the value */
  MeasureType GetValue();

protected:

  MattesMutualInformationImageToImageObjectMetric();
  virtual ~MattesMutualInformationImageToImageObjectMetric();
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Initiates multi-threading to evaluate the current metric value
   * and derivatives. */
//  virtual void GetValueAndDerivativeMultiThreadedInitiate( DerivativeType &
//							   derivativeReturn ) // use superclass

  /* Provide the worker routine to process each point */
  bool GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           itkNotUsed(virtualPoint),
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageDerivativesType &  fixedImageDerivatives,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageDerivativesType & movingImageDerivatives,
                    MeasureType &                      metricValueResult,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType                       threadID);

  /** Default post-processing after multi-threaded calculation of
   * value and derivative.*/
//  virtual void GetValueAndDerivativeMultiThreadedPostProcess( bool doAverage ); // use superclass

private:

  //purposely not implemented
  MattesMutualInformationImageToImageObjectMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

  /** The fixed image marginal PDF */
  typename MarginalPDFType::Pointer m_FixedImageMarginalPDF;

  /** The moving image marginal PDF. */
  typename MarginalPDFType::Pointer m_MovingImageMarginalPDF;

  /** The joint PDF and PDF derivatives. */
  typename JointPDFType::Pointer            m_JointPDF;

  /** Buffer sizes of the PDF and PDF derivatives */
  SizeValueType m_JointPDFBufferSize;

  /** Joint PDF types */
  typedef JointPDFType::IndexType             JointPDFIndexType;
  typedef JointPDFType::PixelType             JointPDFValueType;
  typedef JointPDFType::RegionType            JointPDFRegionType;
  typedef JointPDFType::SizeType              JointPDFSizeType;

  /** Variables to define the marginal and joint histograms. */
  SizeValueType m_NumberOfHistogramBins;
  double        m_MovingImageNormalizedMin;
  double        m_FixedImageNormalizedMin;
  double        m_FixedImageTrueMin;
  double        m_FixedImageTrueMax;
  double        m_MovingImageTrueMin;
  double        m_MovingImageTrueMax;
  double        m_FixedImageBinSize;
  double        m_MovingImageBinSize;

  /** Typedefs for BSpline kernel and derivative functions. */
  typedef BSplineKernelFunction<3>           CubicBSplineFunctionType;
  typedef BSplineDerivativeKernelFunction<3> CubicBSplineDerivativeFunctionType;

  /** Cubic B-Spline kernels */
  typename CubicBSplineFunctionType::Pointer             m_CubicBSplineKernel;
  typename CubicBSplineDerivativeFunctionType::Pointer   m_CubicBSplineDerivativeKernel;

  double m_JointPDFSum;

  /** For threading
  JointPDFType::Pointer * m_ThreaderJointPDF;
  int *m_ThreaderJointPDFStartBin;
  int *m_ThreaderJointPDFEndBin;
  double *m_ThreaderJointPDFSum;
  */
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMattesMutualInformationImageToImageObjectMetric.hxx"
#endif

#endif




/*

  //This is of one two evaluation methods that the user may call.
  void GetValueAndDerivative( MeasureType & valueReturn,
                              DerivativeType & derivativeReturn)
  {
    //1) Do any pre-processing required for your metric. To help with
    // threading, you can use ImageToData or Array1DToData classes,
    // or derive your own from ObjectToData.

    //2) Call GetValueAndDerivativeMultiThreadedInitiate.
    //This will iterate over virtual image region and call your
    // GetValueAndDerivativeProcessPoint method, see definition in
    // base.
    this->GetValueAndDerivativeMultiThreadedInitiate( derivativeReturn );

    //3) Optionally call GetValueAndDerivativeMultiThreadedPostProcess for
    // default post-processing, which sums up results from each thread,
    // and optionally averages them. It then assigns the results to
    // 'value' and 'derivative', without copying in the case of 'derivative'.
    //Do your own post-processing as needed.
    this->GetValueAndDerivativeMultiThreadedPostProcess( true  );//doAverage

    //4) Return the value result. The derivative result has already been
    // written to derivativeReturn.
    valueReturn = this->GetValueResult();

    //That's it. Easy as 1, 2, 3 (and 4).
  }


*/
