/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsImageToImageObjectMetric.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetric_h
#define __itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetric_h

#include "itkANTSNeighborhoodCorrelationImageToImageObjectMetric.h"

namespace itk
{

/** \class ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric
 * \brief Computes normalized cross correlation using a small neighborhood
 * on random sampled voxels between two images. Inherited from
 * ANTSNeighborhoodCorrelationImageToImageObjectMetric. Unlike using a
 * sliding window over all the image region, voxels are randomly sampled.
 *
 *
 *  Example of usage:
 *
 *  typedef itk::ANTSNeighborhoodCorrelationImageToImageObjectMetric<ImageType, ImageType> MetricType;
 *  typedef MetricType::Pointer MetricTypePointer;
 *  MetricTypePointer metric = MetricType::New();
 *
 *  // initialization
 *  metric->SetNumberOfSampling(1000);
 *  Size<Dimension> neighborhoodRadius;
 *  neighborhoodRadius.Fill(2);
 *  metric->SetRadius(neighborhoodRadius);
 *  metric->SetFixedImage(fixedImage);
 *  metric->SetMovingImage(movingImage);
 *  metric->SetFixedTransform(transformFix);
 *  metric->SetMovingTransform(transformMov);
 *
 *  // getting derivative and metric value
 *  metric->GetValueAndDerivative(valueReturn, derivativeReturn);
 *
 *
 * Please cite this reference for more details:
 *
 * Brian B. Avants, Nicholas J. Tustison, Gang Song, Philip A. Cook,
 * Arno Klein, James C. Gee, A reproducible evaluation of ANTs similarity metric
 * performance in brain image registration, NeuroImage, Volume 54, Issue 3,
 * 1 February 2011, Pages 2033-2044, ISSN 1053-8119,
 * DOI: 10.1016/j.neuroimage.2010.09.025.
 *
 * This Class is templated over the type of the two input objects.
 * This is the base class for a hierarchy of similarity metrics that may, in
 * derived classes, operate on meshes, images, etc.  This class computes a
 * value that measures the similarity between the two objects.
 *
 * \ingroup RegistrationMetrics
 *
 */

template <class TFixedImage,
          class TMovingImage,
          class TVirtualImage = TFixedImage >
class ITK_EXPORT ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric :
public ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
{
public:

  /** Standard class typedefs. */
  typedef ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric                      Self;
  typedef ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
                                                              Superclass;
  typedef typename Superclass::Superclass                     MetricBaseclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric, ImageToImageObjectMetric);

  /** superclass types */
  typedef typename Superclass::MeasureType             MeasureType;
  typedef typename Superclass::DerivativeType          DerivativeType;
  typedef typename Superclass::VirtualPointType        VirtualPointType;
  typedef typename Superclass::VirtualIndexType        VirtualIndexType;
  typedef typename Superclass::FixedImagePointType     FixedImagePointType;
  typedef typename Superclass::FixedImagePixelType     FixedImagePixelType;
  typedef typename Superclass::FixedImageDerivativesType
                                                     FixedImageDerivativesType;

  typedef typename Superclass::MovingImagePointType    MovingImagePointType;
  typedef typename Superclass::MovingImagePixelType    MovingImagePixelType;
  typedef typename Superclass::MovingImageDerivativesType
                                                    MovingImageDerivativesType;

  typedef typename Superclass::MovingTransformType     MovingTransformType;
  typedef typename Superclass::MovingTransformJacobianType
                                                  MovingTransformJacobianType;


  typedef typename Superclass::ThreaderInputObjectType ThreaderInputObjectType;
  typedef typename Superclass::FixedImageType          FixedImageType;
  typedef typename Superclass::MovingImageType         MovingImageType;
  typedef typename Superclass::VirtualImageType        VirtualImageType;
  typedef typename Superclass::FixedOutputPointType    FixedOutputPointType;
  typedef typename Superclass::MovingOutputPointType   MovingOutputPointType;


  /** specifics for sliding window based image to image metric **/
  typedef typename VirtualImageType::RegionType VirtualImageRegionType;
//  typedef typename VirtualImageType::IndexType VirtualIndexType;
  typedef typename VirtualImageType::SizeType RadiusType;


  itkStaticConstMacro(VirtualImageDimension, unsigned int,
              ::itk::GetImageDimension<VirtualImageType>::ImageDimension);


  /** Initialize. Must be called before first call to GetValue or
   *  GetValueAndDerivative, after metric settings are changed. */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

  /** Evaluate and return the value and derivative */
  void GetValueAndDerivative( MeasureType & value, DerivativeType & derivative);

  /** Evaluate and return the metric value */
  MeasureType GetValue()
  { itkExceptionMacro("GetValue not yet implemented."); }


  itkSetMacro(NumberOfSampling, unsigned int);
  itkGetMacro(NumberOfSampling, unsigned int);

  itkSetMacro(Radius, RadiusType);
  itkGetMacro(Radius, RadiusType);

protected:

  /* Worker routine to process each point */
  bool SparseGetValueAndDerivativeProcessPoint(
                    const VirtualImageRegionType &scan_region,
                    const VirtualIndexType &           virtualIndex,
                    const VirtualPointType &           virtualPoint,
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageDerivativesType &  fixedImageDerivatives,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageDerivativesType & movingImageDerivatives,
                    MeasureType &                      metricValueResult,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType                       threadID);

  ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric();
  virtual ~ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  //purposely not implemented
  ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

private:

  /**
   * Reimplement the threader callback function from
   *  ImageToImageObjectMetric
   * to add random sampling voxel positions in each thread.
   * Derivatives are computed only at sampled positions.
   * The derivatives at unsampled positions will be kept intact.
   *
   * This callback function is set inside class constructor.
   *
   * Multi-threader callback used to iterate over image region by thread,
   * and call the derived class' user worker method to calculate
   * value and derivative.
   * If a derived class needs to implement its own callback to replace this,
   * define a static method with a different name, and assign it to the
   * threader in the class' constructor by calling
   * \c m_ValueAndDerivativeThreader->SetThreadedGenerateData( mycallback ) */
  static void SparseSamplingGetValueAndDerivativeMultiThreadedCallback(
                          const ThreaderInputObjectType& virtualImageSubRegion,
                          ThreadIdType threadID,
                          MetricBaseclass * self);

  unsigned int m_NumberOfSampling;

  RadiusType m_Radius;

  typedef typename Superclass::ScanningIteratorType ScanningIteratorType;
  typedef typename Superclass::ScanMemType          ScanMemType;
  typedef typename Superclass::ScanParaType         ScanParaType;
  typedef typename Superclass::SumQueueType         SumQueueType;

  // overload this function for arbitrary starting index
  inline void InitializeScanning(
          const VirtualImageRegionType &scan_region,
          const VirtualIndexType &start_index,
          ScanningIteratorType &scan_it,
          ScanMemType &scan_mem,
          ScanParaType &scan_para,
          const ThreadIdType threadID);

protected:
  //sampling routine
  typedef std::vector< VirtualIndexType > VirtualImageSampleContainer;
  /** Uniformly select a sample set from the fixed image domain. */
  void SampleVirtualImageRegion(
          const VirtualImageRegionType &scan_region,
          VirtualImageSampleContainer & samples);

  void FullSampleVirtualImageRegion(
          const VirtualImageRegionType &scan_region,
          VirtualImageSampleContainer & samples);


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetric.txx"
#endif

#endif
