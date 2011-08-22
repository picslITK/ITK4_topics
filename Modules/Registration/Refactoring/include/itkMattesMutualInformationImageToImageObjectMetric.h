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
  typedef typename Superclass::FixedGradientPixelType   FixedImageGradientsType;
  typedef typename Superclass::MovingImagePointType        MovingImagePointType;
  typedef typename Superclass::MovingImagePixelType        MovingImagePixelType;
  typedef typename Superclass::MovingGradientPixelType  MovingImageGradientsType;
  /** Value type of the PDF */
  typedef double PDFValueType;
  /** Typedef for the joint PDF and marginal PDF are stored as ITK Images. */
  typedef Image<PDFValueType,1> MarginalPDFType;
  typedef typename MarginalPDFType::IndexType MarginalPDFIndexType;
  typedef typename MarginalPDFType::PointType MarginalPDFPointType;
  typedef Image< PDFValueType, 2 >            JointPDFType;
  typedef typename JointPDFType::PointType JointPDFPointType;
  itkGetConstReferenceMacro(JointPDF,typename JointPDFType::Pointer);

  // Declare the type for the derivative calculation
  typedef itk::GradientRecursiveGaussianImageFilter< JointPDFType >  JPDFGradientFilterType;
  typedef typename JPDFGradientFilterType::OutputImageType JPDFGradientImageType;
  typedef typename JPDFGradientImageType::Pointer JPDFGradientImagePointer;
  typedef itk::GradientRecursiveGaussianImageFilter< MarginalPDFType >  MarginalGradientFilterType;
  typedef typename MarginalGradientFilterType::OutputImageType MarginalGradientImageType;
  typedef typename MarginalGradientImageType::Pointer MarginalGradientImagePointer;

  itkSetClampMacro( NumberOfHistogramBins, SizeValueType,
                    5, NumericTraits< SizeValueType >::max() );
  itkGetConstReferenceMacro(NumberOfHistogramBins, SizeValueType);
  /** pdf interpolator */
//  typedef BSplineInterpolateImageFunction<JointPDFType,double> JointPDFInterpolatorType;
  typedef LinearInterpolateImageFunction<JointPDFType,double> JointPDFInterpolatorType;
  typedef typename JointPDFInterpolatorType::Pointer JointPDFInterpolatorPointer;
//  typedef BSplineInterpolateImageFunction<MarginalPDFType,double> MarginalPDFInterpolatorType;
  typedef LinearInterpolateImageFunction<MarginalPDFType,double> MarginalPDFInterpolatorType;
  typedef typename MarginalPDFInterpolatorType::Pointer MarginalPDFInterpolatorPointer;

  /** Initialize the metric. Make sure all essential inputs are plugged in. */
  virtual void Initialize() throw (itk::ExceptionObject);

  /** Get both the value and derivative intializes the processing.
   *  For Mattes MI, we just compute the joint histogram / pdf here.
   *  This implementation single-threads the JH computation but it
   *  could be multi-threaded in the future.
   */
  void GetValueAndDerivative(MeasureType & value, DerivativeType & derivative);

  /** Get the value */
  MeasureType GetValue();

protected:

  MattesMutualInformationImageToImageObjectMetric();
  virtual ~MattesMutualInformationImageToImageObjectMetric();
  void PrintSelf(std::ostream & os, Indent indent) const;

  bool ComputeJointPDFPoint( FixedImagePixelType fixedImageValue, MovingImagePixelType movingImageValue , JointPDFPointType& jointPDFpoint , unsigned int threadID ) {
    double a=(fixedImageValue-this->m_FixedImageTrueMin)/(this->m_FixedImageTrueMax-this->m_FixedImageTrueMin);
    double b=(movingImageValue-this->m_MovingImageTrueMin)/(this->m_MovingImageTrueMax-this->m_MovingImageTrueMin);
    jointPDFpoint[0]=a;
    jointPDFpoint[1]=b;
    return this->m_ThreaderJointPDFInterpolator[threadID]->IsInsideBuffer(jointPDFpoint );
   }

  inline double ComputeFixedImageMarginalPDFDerivative( MarginalPDFPointType margPDFpoint , unsigned int threadID )
  {
    double offset=0.5*this->m_JointPDFSpacing[0];
    double eps=this->m_JointPDFSpacing[0];
    MarginalPDFPointType  leftpoint=margPDFpoint;
    leftpoint[0]-=offset;
    MarginalPDFPointType  rightpoint=margPDFpoint;
    rightpoint[0]+=offset;
    if (leftpoint[0] < eps ) leftpoint[0]=eps;
    if (rightpoint[0] < eps ) rightpoint[0]=eps;
    if (leftpoint[0] > 1 ) leftpoint[0]=1;
    if (rightpoint[0] > 1  ) rightpoint[0]=1;
    double delta=rightpoint[0]-leftpoint[0];
    if ( delta > 0 ) {
    double deriv=this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->Evaluate(rightpoint)-
                 this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->Evaluate(leftpoint);
    return deriv/delta;
    }
    else return 0;
  }

  inline double ComputeMovingImageMarginalPDFDerivative( MarginalPDFPointType margPDFpoint , unsigned int threadID )
  {
    double offset=0.5*this->m_JointPDFSpacing[0];
    double eps=this->m_JointPDFSpacing[0];
    MarginalPDFPointType  leftpoint=margPDFpoint;
    leftpoint[0]-=offset;
    MarginalPDFPointType  rightpoint=margPDFpoint;
    rightpoint[0]+=offset;
    if (leftpoint[0] < eps ) leftpoint[0]=eps;
    if (rightpoint[0] < eps ) rightpoint[0]=eps;
    if (leftpoint[0] > 1 ) leftpoint[0]=1;
    if (rightpoint[0] > 1  ) rightpoint[0]=1;
    double delta=rightpoint[0]-leftpoint[0];
    if ( delta > 0 ) {
    double deriv=this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->Evaluate(rightpoint)-
                 this->m_ThreaderMovingImageMarginalPDFInterpolator[threadID]->Evaluate(leftpoint);
    return deriv/delta;
    }
    else return 0;
  }

  inline double ComputeJointPDFDerivative( JointPDFPointType jointPDFpoint , unsigned int threadID , unsigned int ind  )
  {
    double offset=0.5*this->m_JointPDFSpacing[ind];
    double eps=this->m_JointPDFSpacing[ind];
    JointPDFPointType  leftpoint=jointPDFpoint;
    leftpoint[ind]-=offset;
    JointPDFPointType  rightpoint=jointPDFpoint;
    rightpoint[ind]+=offset;
    if (leftpoint[ind] < eps ) leftpoint[ind]=eps;
    if (rightpoint[ind] < eps ) rightpoint[ind]=eps;
    if (leftpoint[ind] > 1 ) leftpoint[ind]=1;
    if (rightpoint[ind] > 1 ) rightpoint[ind]=1;
    double delta=rightpoint[ind]-leftpoint[ind];
    double deriv=0;
    if ( delta > 0 ) {
      deriv=this->m_ThreaderJointPDFInterpolator[threadID]->Evaluate(rightpoint)-
            this->m_ThreaderJointPDFInterpolator[threadID]->Evaluate(leftpoint);
      return deriv/delta;
    }
    else return deriv;
  }


  /** Initiates multi-threading to evaluate the current metric value
   * and derivatives. */
//  virtual void GetValueAndDerivativeMultiThreadedInitiate( DerivativeType &
//							   derivativeReturn ) // use superclass

  /* Provide the worker routine to process each point */
  bool GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           itkNotUsed(virtualPoint),
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageGradientsType &  fixedImageGradients,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageGradientsType & movingImageGradients,
                    MeasureType &                      metricValueResult,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType                       threadID);

  /** Default post-processing after multi-threaded calculation of
   * value and derivative.*/
//  virtual void GetValueAndDerivativeMultiThreadedPostProcess( bool doAverage ); // use superclass

  void EnforceJointHistogramBoundaryConditions();

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
  typedef JointPDFType::SpacingType           JointPDFSpacingType;

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

  double m_JointPDFSum;
  JointPDFSpacingType m_JointPDFSpacing;

  /** For threading */
  JointPDFInterpolatorPointer* m_ThreaderJointPDFInterpolator;
  MarginalPDFInterpolatorPointer* m_ThreaderFixedImageMarginalPDFInterpolator;
  MarginalPDFInterpolatorPointer* m_ThreaderMovingImageMarginalPDFInterpolator;
  double m_Log2;
  unsigned int m_Padding;
  JPDFGradientImagePointer m_JPDFGradientImage;
  MarginalGradientImagePointer m_MarginalGradientImage;
  /*
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
