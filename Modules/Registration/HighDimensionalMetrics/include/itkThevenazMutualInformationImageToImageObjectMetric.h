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

#ifndef __itkThevenazMutualInformationImageToImageObjectMetric_h
#define __itkThevenazMutualInformationImageToImageObjectMetric_h

#include "itkImageToImageObjectMetric.h"
#include "itkImage.h"
#include "itkBSplineDerivativeKernelFunction.h"

namespace itk
{
/** \class ThevenazMutualInformationImageToImageMetric
 * \brief Computes the mutual information between two images to be
 * registered using the method of Mattes et al.
 *
 * References:
 * [1] "Optimization of Mutual Information for MultiResolution Image
 *      Registration"
 *      P. Thevenaz and M. Unser
 *      IEEE Transactions in Image Processing, 9(12) December 2000.
 *
 * \ingroup ITKHighDimensionalMetrics
 */

template<class TFixedImage,class TMovingImage,class TVirtualImage = TFixedImage>
class ITK_EXPORT ThevenazMutualInformationImageToImageObjectMetric :
  public ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
{
public:

  /** Standard class typedefs. */
  typedef ThevenazMutualInformationImageToImageObjectMetric     Self;
  typedef ImageToImageObjectMetric<TFixedImage, TMovingImage>   Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ThevenazMutualInformationImageToImageObjectMetric,
                ImageToImageObjectMetric);

  /** Type used for representing parameter values  */
  typedef typename Superclass::CoordinateRepresentationType
                                                  CoordinateRepresentationType;
  /** Type used internally for computations */
  typedef typename Superclass::InternalComputationValueType
                                                  InternalComputationValueType;
  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  /** Superclass typedefs */
  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;
  typedef typename Superclass::VirtualPointType         VirtualPointType;
  typedef typename Superclass::FixedImagePointType      FixedImagePointType;
  typedef typename Superclass::FixedImagePixelType      FixedImagePixelType;
  typedef typename Superclass::FixedGradientPixelType
                                                        FixedImageGradientsType;
  typedef typename Superclass::MovingImagePointType     MovingImagePointType;
  typedef typename Superclass::MovingImagePixelType     MovingImagePixelType;
  typedef typename Superclass::MovingGradientPixelType
                                                        MovingImageGradientsType;
  /** Value type of the PDF */
  typedef InternalComputationValueType                  PDFValueType;

  /** Typedef for the joint PDF and marginal PDF are stored as ITK Images. */
  typedef Image<PDFValueType,1>                 MarginalPDFType;
  typedef typename MarginalPDFType::IndexType   MarginalPDFIndexType;
  typedef typename MarginalPDFType::PointType   MarginalPDFPointType;
  typedef Image< PDFValueType, 2>               JointPDFType;
  typedef typename JointPDFType::IndexType      JointPDFIndexType;
  typedef typename JointPDFType::PointType      JointPDFPointType;
  typedef typename JointPDFType::IndexValueType JointPDFIndexValueType;

  itkGetConstReferenceMacro(JointPDF,typename JointPDFType::Pointer);

  // Declare the type for the derivative calculation
  typedef itk::GradientRecursiveGaussianImageFilter< JointPDFType >
                                                         JPDFGradientFilterType;

  typedef typename JPDFGradientFilterType::OutputImageType JPDFGradientImageType;

  typedef typename JPDFGradientImageType::Pointer JPDFGradientImagePointer;

  typedef itk::GradientRecursiveGaussianImageFilter< MarginalPDFType >
                                                  MarginalGradientFilterType;
  typedef typename MarginalGradientFilterType::OutputImageType
                                                  MarginalGradientImageType;
  typedef typename MarginalGradientImageType::Pointer
                                                  MarginalGradientImagePointer;

  /** pdf interpolator */
  typedef LinearInterpolateImageFunction<JointPDFType,double>
                                                     JointPDFInterpolatorType;
  typedef typename JointPDFInterpolatorType::Pointer JointPDFInterpolatorPointer;
  typedef LinearInterpolateImageFunction<MarginalPDFType,double>
                                                     MarginalPDFInterpolatorType;
  typedef typename MarginalPDFInterpolatorType::Pointer
                                                  MarginalPDFInterpolatorPointer;

  /** Get/Set the number of histogram bins */
  itkSetClampMacro( NumberOfHistogramBins, SizeValueType,
                    5, NumericTraits< SizeValueType >::max() );
  itkGetConstReferenceMacro(NumberOfHistogramBins, SizeValueType );

  /** Get/Set option to smooth the joint pdf after it's updated */
  itkSetMacro(DoSmoothJointPDF, bool);
  itkGetConstReferenceMacro(DoSmoothJointPDF, bool);
  itkBooleanMacro(DoSmoothJointPDF);

  /** Initialize the metric. Make sure all essential inputs are plugged in. */
  virtual void Initialize() throw (itk::ExceptionObject);

  /** Get both the value and derivative intializes the processing.
   *  For Mattes MI, we just compute the joint histogram / pdf here.
   *  This implementation single-threads the JH computation but it
   *  could be multi-threaded in the future.
   *  Results are returned in \c value and \c derivative.
   */
  using Superclass::GetValueAndDerivative;
  void GetValueAndDerivative(MeasureType & value, DerivativeType & derivative);

  /** Get the value */
  using Superclass::GetValue;
  MeasureType GetValue();

protected:

  ThevenazMutualInformationImageToImageObjectMetric();
  virtual ~ThevenazMutualInformationImageToImageObjectMetric();
  void PrintSelf(std::ostream & os, Indent indent) const;

  bool ComputeJointPDFPoint( const FixedImagePixelType fixedImageValue,
                             const MovingImagePixelType movingImageValue,
                             JointPDFPointType& jointPDFpoint,
                             const ThreadIdType threadID );

  inline InternalComputationValueType ComputeFixedImageMarginalPDFDerivative(
                                        const MarginalPDFPointType margPDFpoint,
                                        const ThreadIdType threadID );

  inline InternalComputationValueType ComputeMovingImageMarginalPDFDerivative(
                                        const MarginalPDFPointType margPDFpoint,
                                        const ThreadIdType threadID );

  inline InternalComputationValueType ComputeJointPDFDerivative(
                                          const JointPDFPointType jointPDFpoint,
                                          const ThreadIdType threadID,
                                          const SizeValueType ind );

  bool GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           virtualPoint,
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageGradientsType &    fixedImageGradients,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageGradientsType &   movingImageGradients,
                    MeasureType &                      metricValueResult,
                    DerivativeType &                   localDerivativeReturn,
                    const ThreadIdType                 threadID);

  void EnforceJointHistogramBoundaryConditions();

  /** Update the histograms for use in GetValueAndDerivative */
  void UpdateHistograms();

private:

  //purposely not implemented
  ThevenazMutualInformationImageToImageObjectMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

  /** The fixed image marginal PDF */
  typename MarginalPDFType::Pointer m_FixedImageMarginalPDF;

  /** The moving image marginal PDF. */
  typename MarginalPDFType::Pointer m_MovingImageMarginalPDF;

  /** The joint PDF and PDF derivatives. */
  typename JointPDFType::Pointer            m_JointPDF;

  /** Flag to control smoothing of joint pdf */
  bool                                      m_DoSmoothJointPDF;

  /** Joint PDF types */
  typedef typename JointPDFType::PixelType             JointPDFValueType;
  typedef typename JointPDFType::RegionType            JointPDFRegionType;
  typedef typename JointPDFType::SizeType              JointPDFSizeType;
  typedef typename JointPDFType::SpacingType           JointPDFSpacingType;

  /** Variables to define the marginal and joint histograms. */
  SizeValueType                       m_NumberOfHistogramBins;
  InternalComputationValueType        m_FixedImageTrueMin;
  InternalComputationValueType        m_FixedImageTrueMax;
  InternalComputationValueType        m_MovingImageTrueMin;
  InternalComputationValueType        m_MovingImageTrueMax;
  InternalComputationValueType        m_FixedImageBinSize;
  InternalComputationValueType        m_MovingImageBinSize;

  InternalComputationValueType              m_JointPDFSum;
  JointPDFSpacingType                       m_JointPDFSpacing;

  /** For threading */
  JointPDFInterpolatorPointer*    m_ThreaderJointPDFInterpolator;
  MarginalPDFInterpolatorPointer* m_ThreaderFixedImageMarginalPDFInterpolator;
  MarginalPDFInterpolatorPointer* m_ThreaderMovingImageMarginalPDFInterpolator;

  InternalComputationValueType    m_Log2;
  JointPDFIndexValueType          m_Padding;

  JPDFGradientImagePointer      m_JPDFGradientImage;
  MarginalGradientImagePointer  m_MarginalGradientImage;

  /*
  JointPDFType::Pointer * m_ThreaderJointPDF;
  int *m_ThreaderJointPDFStartBin;
  int *m_ThreaderJointPDFEndBin;
  double *m_ThreaderJointPDFSum;
  */

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkThevenazMutualInformationImageToImageObjectMetric.hxx"
#endif

#endif
