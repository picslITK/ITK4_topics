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
#ifndef __itkImageToImageNeighborhoodNormalizedCrossCorrelationFunction_h
#define __itkImageToImageNeighborhoodNormalizedCrossCorrelationFunction_h

#include "itkObjectToObjectMetric.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImage.h"
#include "itkCovariantVector.h"
#include <deque>

namespace itk
{

/** \class itkImageToImageNeighborhoodNormalizedCrossCorrelation
 * \brief Computes normalized cross correlation using a small neighborhood
 * for each voxel between two images
 *
 * Use on-the-fly queues to compute sliding windows
 *
 * This Class is templated over the type of the two input objects.
 * This is the base class for a hierarchy of similarity metrics that may, in
 * derived classes, operate on meshes, images, etc.  This class computes a
 * value that measures the similarity between the two objects.
 *
 * \ingroup RegistrationMetrics
 *
 */
template< class TFixedImage,  class TMovingImage >
class ITK_EXPORT ImageToImageNeighborhoodNormalizedCrossCorrelationFunction:
public ObjectToObjectMetric<TFixedImage, TMovingImage>
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageNeighborhoodNormalizedCrossCorrelationFunction       Self;
  typedef ObjectToObjectMetric<TFixedImage, TMovingImage>   Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToImageNeighborhoodNormalizedCrossCorrelationFunction, ObjectToObjectMetric);

  /**  Type of the measure. */
  typedef typename Superclass::MeasureType    MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType DerivativeType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;


  typedef double  InternalComputationValueType;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  // virtual void Initialize(void) throw ( ExceptionObject ) = 0;


  /** specifics for sliding window based image to image metric **/
  typedef TFixedImage   FixedImageType;
  typedef TMovingImage  MovingImageType;
  typedef typename FixedImageType::Pointer FixedImagePointerType;
  typedef typename MovingImageType::Pointer MovingImagePointerType;
  typedef typename FixedImageType::RegionType ImageRegionType;
  typedef typename FixedImageType::SizeType RadiusType;
  typedef typename FixedImageType::IndexType IndexType;


//  void SetFixedImage(FixedImagePointerType &fixedImage) { m_FixedImage = fixedImage; }
//  FixedImagePointerType GetFixedImage() { return m_FixedImage; }

  itkSetMacro(FixedImage, FixedImagePointerType);
  itkGetMacro(FixedImage, FixedImagePointerType);
  itkSetMacro(MovingImage, MovingImagePointerType);
  itkGetMacro(MovingImage, MovingImagePointerType);

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> FixedImageGradientCalculatorType;
  typedef typename FixedImageGradientCalculatorType::Pointer   FixedImageGradientCalculatorPointerType;
  typedef CentralDifferenceImageFunction<MovingImageType> MovingImageGradientCalculatorType;
  typedef typename MovingImageGradientCalculatorType::Pointer   MovingImageGradientCalculatorPointerType;



  itkStaticConstMacro(ImageDimension, unsigned int,
      ::itk::GetImageDimension< FixedImageType>::ImageDimension);

  typedef CovariantVector<double,ImageDimension> CovariantVectorType;

  virtual unsigned int GetNumberOfParameters() const;

  virtual MeasureType GetValue( const ParametersType & parameters ) const;

  virtual void GetDerivative(const ParametersType &,
      DerivativeType & derivative) const;

  virtual void Initialize(void) throw ( itk::ExceptionObject );

  typedef Vector<InternalComputationValueType, ImageDimension> VectorType;

  itkSetMacro(Radius, RadiusType);
  itkGetMacro(Radius, RadiusType);

public:
  // interested values here updated during scanning

  typedef InternalComputationValueType QUEUEREALTYPE;
  typedef std::deque<QUEUEREALTYPE> SumQueueType;
  typedef ConstNeighborhoodIterator<FixedImageType> ScanningIteratorType;

  typedef struct ScanMemType{
      // queus used in the scanning
      SumQueueType Qsuma2;
      SumQueueType Qsumb2;
      SumQueueType Qsuma;
      SumQueueType Qsumb;
      SumQueueType Qsumab;
      SumQueueType Qcount;


      QUEUEREALTYPE Ia;
      QUEUEREALTYPE Ja;
      QUEUEREALTYPE sfm;
      QUEUEREALTYPE sff;
      QUEUEREALTYPE smm;
  } ScanMemType;

  typedef struct ScanParaType{
      // const values during scanning
      ImageRegionType scan_region;
      int number_of_fill_zero; // for each queue
      unsigned int window_length; // number of voxels in the scanning window
      int scan_region_begin_index_dim0;

      FixedImagePointerType I;
      MovingImagePointerType J;
      RadiusType r;

  } ScanParaType;


  // computation routines for normalized cross correlation
public:


  virtual void InitializeGradientCalculator();

  virtual inline void InitializeScanning(
      const ImageRegionType &scan_region,
      ScanningIteratorType &scan_it,
      ScanMemType &scan_mem,
      ScanParaType &scan_para);

  virtual inline void UpdateQueuesAtBeginingOfLine(
      const ScanningIteratorType &scan_it,
      ScanMemType &scan_mem,
      const ScanParaType &scan_para);

  virtual inline void UpdateQueuesToNextScanWindow(
      const ScanningIteratorType &scan_it,
      ScanMemType &scan_mem,
      const ScanParaType &scan_para);

  virtual inline void UpdateQueues(
      const ScanningIteratorType &scan_it,
      ScanMemType &scan_mem,
      const ScanParaType &scan_para);

  virtual inline void ComputeInformationFromQueues(
      const ScanningIteratorType &scan_it,
      ScanMemType &scan_mem,
      const ScanParaType &scan_para);

  virtual void ComputeUpdateBothDirection(
      const ScanningIteratorType &scan_it,
      ScanMemType &scan_mem,
      const ScanParaType &scan_para,
      VectorType &update,
      VectorType &updateInv,
      InternalComputationValueType &local_cc);



protected:
  ImageToImageNeighborhoodNormalizedCrossCorrelationFunction();
  virtual ~ImageToImageNeighborhoodNormalizedCrossCorrelationFunction();

  virtual void PrintSelf(std::ostream & os, Indent indent) const;

  mutable ParametersType      m_Parameters;

  FixedImagePointerType m_FixedImage;
  MovingImagePointerType m_MovingImage;

  FixedImageGradientCalculatorPointerType       m_FixedImageGradientCalculator;
  MovingImageGradientCalculatorPointerType       m_MovingImageGradientCalculator;

  RadiusType m_Radius;

private:
  ImageToImageNeighborhoodNormalizedCrossCorrelationFunction(const Self &); //purposely not implemented
  void operator=(const Self &);     //purposely not implemented

};


// functor for threading using the metric function class
// assuming function has output allocated already
template<class TMetricFunction, class TDeformationField>
struct MetricThreadedHolder{

  typedef MetricThreadedHolder          Self;

  typedef TMetricFunction           MetricType;
  typedef typename MetricType::Pointer  MetricTypePointer;
  typedef TDeformationField             DeformationFieldType;
  typedef typename DeformationFieldType::Pointer DeformationFieldPointerType;
  typedef typename MetricType::MeasureType  MeasureType;
  typedef typename MetricType::InternalComputationValueType InternalComputationValueType;
  typedef typename MetricType::ImageRegionType ImageRegionType;

//  itkSetMacro(MetricFunction, MetricTypePointer);
//  itkGetMacro(MetricFunction, MetricTypePointer);
//
//  itkSetMacro(UpdateField, DeformationFieldPointerType);
//  itkGetMacro(UpdateField, DeformationFieldPointerType);
//
//  itkSetMacro(UpdateFieldInv, DeformationFieldPointerType);
//  itkGetMacro(UpdateFieldInv, DeformationFieldPointerType);


public:
  MetricTypePointer           metric;
  DeformationFieldPointerType updateField;
  DeformationFieldPointerType updateFieldInv;
  std::vector<InternalComputationValueType> measure_per_thread;

private:


public:
  InternalComputationValueType AccumulateMeasuresFromAllThreads() {
    InternalComputationValueType energy = NumericTraits<InternalComputationValueType>::Zero;
    for(unsigned int i=0; i<measure_per_thread.size(); i++) energy += measure_per_thread[i];
    return energy;
  }

  static void ComputeMetricValueInRegionOnTheFlyThreaded(const ImageRegionType &regionForThread, int threadId,  Self *holder){

    typename MetricType::ScanningIteratorType scan_it;
    typename MetricType::ScanParaType scan_para;
    typename MetricType::ScanMemType scan_mem;
    typename MetricType::VectorType deriv, derivInv;

    InternalComputationValueType local_cc;

    typedef itk::ImageRegionIterator<DeformationFieldType> UpdateIteratorType;

    UpdateIteratorType       nU(holder->updateField,  regionForThread);
    UpdateIteratorType       nUinv(holder->updateFieldInv, regionForThread);

    nU.GoToBegin();
    nUinv.GoToBegin();

    holder->metric->InitializeScanning(regionForThread, scan_it, scan_mem, scan_para);

    holder->measure_per_thread[threadId] = NumericTraits<InternalComputationValueType>::Zero;

    scan_it.GoToBegin();
    for(; !scan_it.IsAtEnd(); ++scan_it){
      holder->metric->UpdateQueues(scan_it, scan_mem, scan_para);
      holder->metric->ComputeInformationFromQueues(scan_it, scan_mem, scan_para);
      holder->metric->ComputeUpdateBothDirection(scan_it, scan_mem, scan_para, deriv, derivInv, local_cc);

      nU.Value() += deriv;
      nUinv.Value() += derivInv;;

      holder->measure_per_thread[threadId] += local_cc;

      ++nU;
      ++nUinv;

    }


  }


};





} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToImageNeighborhoodNormalizedCrossCorrelationFunction.txx"
#endif

#endif
