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
#ifndef __itkDemonsImageToImageMetric_h
#define __itkDemonsImageToImageMetric_h

#include "itkObjectToObjectMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkIndex.h"
#include "itkTransform.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMultiThreader.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"


namespace itk
{

template <class TFixedImage,class TMovingImage >
class ITK_EXPORT DemonsImageToImageMetric :
public ObjectToObjectMetric<TFixedImage, TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef DemonsImageToImageMetric                     Self;
  typedef ObjectToObjectMetric<TFixedImage, TMovingImage>   Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  typedef double  InternalComputationValueType;

  /** Image-accessor typedefs */
  typedef TFixedImage   FixedImageType;
  typedef typename FixedImageType::PixelType FixedImagePixelType;
  typedef TMovingImage  MovingImageType;
  typedef typename MovingImageType::PixelType MovingImagePixelType;
  typedef typename FixedImageType::Pointer FixedImagePointer;
  typedef typename MovingImageType::Pointer MovingImagePointer;
  typedef typename FixedImageType::RegionType RegionType;
  typedef typename RegionType::SizeType SizeType;
  typedef typename FixedImageType::SpacingType SpacingType;
  typedef typename FixedImageType::PointType OriginType;
  typedef typename FixedImageType::PointType PointType;
  typedef typename FixedImageType::DirectionType DirectionType;
  typedef typename FixedImageType::SizeType RadiusType;
  typedef typename FixedImageType::IndexType IndexType;
  typedef typename FixedImageType::RegionType ImageRegionType;
  itkStaticConstMacro(FixedImageDimension, unsigned int,
      ::itk::GetImageDimension< FixedImageType>::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int,
      ::itk::GetImageDimension< MovingImageType>::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DemonsImageToImageMetric, ObjectToObjectMetric);

  /** Transform types */
  typedef double CoordinateRepresentationType;
  typedef itk::Transform< CoordinateRepresentationType,
                     itkGetStaticConstMacro(MovingImageDimension),
                     itkGetStaticConstMacro(FixedImageDimension) >
  TransformType;
  typedef typename TransformType::Pointer         TransformPointer;
  typedef typename TransformType::InputPointType  InputPointType;
  typedef typename TransformType::OutputPointType OutputPointType;
  typedef typename TransformType::ParametersType  TransformParametersType;
  typedef typename TransformType::JacobianType    TransformJacobianType;
  /**  Type of the Interpolator Base class */
  typedef LinearInterpolateImageFunction< FixedImageType, CoordinateRepresentationType > FixedInterpolatorType;
  typedef LinearInterpolateImageFunction< MovingImageType, CoordinateRepresentationType > MovingInterpolatorType;
  typedef typename FixedInterpolatorType::Pointer FixedInterpolatorPointer;
  typedef typename MovingInterpolatorType::Pointer MovingInterpolatorPointer;

  /**  Type of the measure. */
  typedef typename Superclass::MeasureType    MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType DerivativeType;
  typedef CovariantVector<double,FixedImageDimension> FixedCovariantVectorType;
  typedef CovariantVector<double,MovingImageDimension> MovingCovariantVectorType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  /** Set/get images */
  itkSetMacro(FixedImage, FixedImagePointer);
  itkGetMacro(FixedImage, FixedImagePointer);
  itkSetMacro(MovingImage, MovingImagePointer);
  itkGetMacro(MovingImage, MovingImagePointer);


  virtual unsigned int GetNumberOfParameters() const { return 0; }

  virtual MeasureType GetValue( const ParametersType & parameters ) const { return 0; }

  virtual void GetDerivative(const ParametersType &,
      DerivativeType & derivative) const {}

  /** This function computes the local voxel-wise contribution of 
   *  the metric to the global integral of the metric/derivative.
   */
  double ComputeLocalContributionToMetricAndDerivative(PointType mappedFixedPoint, PointType mappedMovingPoint) 
  {
    double metricval=0;
    /** Only the voxelwise contribution given the point pairs. */
    
    FixedImagePixelType fpix=this->m_FixedInterpolator->Evaluate(mappedFixedPoint); 
    MovingImagePixelType mpix=this->m_MovingInterpolator->Evaluate(mappedMovingPoint);
    FixedImagePixelType diff = fpix - mpix ;
    
    // Jacobian should be evaluated at the unmapped (fixed image) point.
    const TransformJacobianType & jacobian = this->m_MovingImageTransform->GetJacobian(mappedMovingPoint);
    
    for ( unsigned int par = 0; par < this->m_MovingImageTransform->GetNumberOfParameters(); par++ )  
    {
      double sum = 0.0;
      for (unsigned int c=0; c < this->m_InputImageVectorLength; c++) 
      {
        metricval+=fabs(diff[c])/(double)FixedImageDimension;
        for ( unsigned int dim = 0; dim < MovingImageDimension; dim++ )  
        {
	  sum += 2.0 *diff[c]*jacobian(dim, par);// * movingImageGradientValue[dim];
        }
      }
    }
    
    return metricval;
  }

  virtual void Initialize(void) throw ( itk::ExceptionObject );

  

  /** This function is calls the ComputeMetricAndDerivative() function
   *  over the domain of interest.  
   */
  double ComputeMetricAndDerivative(const ImageRegionType &thread_region) 
  {
    /** For each location in the virtual domain, map to both the fixed and moving space
     *  and compute the values of the voxels in the corresponding locations.  There should 
     *  be a transform between the virtual space and the fixed/moving space s.t. the images 
     *  are interpolated in an unbiased manner. 
     */
    double metric_sum=0;
    unsigned long ct=0;
    ImageRegionConstIteratorWithIndex<FixedImageType> ItV( this->m_VirtualImage,
        thread_region );
    ItV.GoToBegin();
    while( !ItV.IsAtEnd() )
    {
      /** use the fixed and moving transforms to compute the corresponding points.*/
      bool sampleOk = true;
      // convert the index to a point 
      PointType mappedPoint;
      PointType mappedFixedPoint;
      PointType mappedMovingPoint;
      this->m_VirtualImage->TransformIndexToPhysicalPoint(ItV.GetIndex(),mappedPoint);
      mappedFixedPoint = this->m_FixedImageTransform->TransformPoint(mappedPoint);
      mappedMovingPoint = this->m_MovingImageTransform->TransformPoint(mappedPoint);
      if ( !this->m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint) ||  
           !this->m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint) ) 
	sampleOk=false;
      if ( sampleOk )
        {
	  double metricval=this->ComputeLocalContributionToMetricAndDerivative(mappedFixedPoint,mappedMovingPoint);
          metric_sum+=metricval;	  
	  ct++;
        }
      ++ItV;
    }
    if ( ct > 0 ) {
      std::cout << " metric_sum " << metric_sum << " ct " << ct << " thread_region " << thread_region.GetIndex() << " sz " << thread_region.GetSize() << std::endl;
      return metric_sum; 
    }
    else return 0;
  }

  /** Define the virtual reference space. This space defines the resolution, 
   *  at which the registration is performed as well as the physical coordinate
   *  system.  Useful for unbiased registration.  If the user does not set this 
   *  explicitly then it is taken from the fixed image.  One can also use this 
   *  functionality to control multi-threading. 
   */  
  itkSetMacro(VirtualDomainSpacing, SpacingType);
  itkGetMacro(VirtualDomainSpacing, SpacingType);
  itkSetMacro(VirtualDomainSize, SizeType);
  itkGetMacro(VirtualDomainSize, SizeType);
  itkSetMacro(VirtualDomainIndex, IndexType);
  itkGetMacro(VirtualDomainIndex, IndexType);
  itkSetMacro(VirtualDomainOrigin, OriginType);
  itkGetMacro(VirtualDomainOrigin, OriginType);
  itkSetMacro(VirtualDomainDirection, DirectionType);
  itkGetMacro(VirtualDomainDirection, DirectionType);


  /** Connect the Transform. */
  itkSetObjectMacro(FixedImageTransform, TransformType);
  /** Get a pointer to the Transform.  */
  itkGetConstObjectMacro(FixedImageTransform, TransformType);
  /** Connect the Transform. */
  itkSetObjectMacro(MovingImageTransform, TransformType);
  /** Get a pointer to the Transform.  */
  itkGetConstObjectMacro(MovingImageTransform, TransformType);

protected:

  DemonsImageToImageMetric();
  virtual ~DemonsImageToImageMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  //purposely not implemented
  DemonsImageToImageMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);


  FixedImagePointer m_FixedImage;
  TransformPointer m_FixedImageTransform;
  MovingImagePointer m_MovingImage;
  TransformPointer m_MovingImageTransform;
  FixedImagePointer m_VirtualImage;
  MeasureType*     m_ThreaderMSE;
  DerivativeType*  m_ThreaderDerivatives;
  double           m_Normalizer;

  unsigned int m_NumberOfThreads;
  unsigned int m_NumberOfParameters;

  IndexType  m_VirtualDomainIndex;
  SizeType  m_VirtualDomainSize;
  SpacingType  m_VirtualDomainSpacing;
  OriginType  m_VirtualDomainOrigin;
  DirectionType  m_VirtualDomainDirection;

  FixedInterpolatorPointer m_FixedInterpolator;
  MovingInterpolatorPointer m_MovingInterpolator;

  unsigned int m_InputImageVectorLength;

};


// functor for threading using the metric function class
// assuming function has output allocated already
template<class TMetricFunction, class TDeformationField>
struct DemonsMetricThreadedHolder{

  typedef DemonsMetricThreadedHolder          Self;

  typedef TMetricFunction           MetricType;
  typedef typename MetricType::Pointer  MetricTypePointer;
  typedef TDeformationField             DeformationFieldType;
  typedef typename DeformationFieldType::Pointer DeformationFieldPointer;
  typedef typename MetricType::MeasureType  MeasureType;
  typedef typename MetricType::InternalComputationValueType InternalComputationValueType;
  typedef typename MetricType::RegionType ImageRegionType;
  typedef typename MetricType::FixedImageType ImageType;
  typedef typename MetricType::FixedImagePointer FixedImagePointer;
  typedef typename MetricType::MovingImagePointer MovingImagePointer;
  typedef typename MetricType::TransformPointer TransformPointer;

public:
  MetricTypePointer           metric;
  FixedImagePointer fixed_image;
  MovingImagePointer moving_image;
  TransformPointer transformF;
  TransformPointer transformM;
  std::vector<InternalComputationValueType> measure_per_thread;

  InternalComputationValueType AccumulateMeasuresFromAllThreads() {
    InternalComputationValueType energy = NumericTraits<InternalComputationValueType>::Zero;
    for(unsigned int i=0; i<measure_per_thread.size(); i++) energy += measure_per_thread[i];
    return energy;
  }

  static void ComputeMetricValueInRegionOnTheFlyThreaded(const ImageRegionType &regionForThread, int threadId,  Self *holder){

    //    std::cout << regionForThread << std::endl;
    InternalComputationValueType local_metric;
    holder->measure_per_thread[threadId] = NumericTraits<InternalComputationValueType>::Zero;
    /** Compute one iteration of the metric */
    local_metric=holder->metric->ComputeMetricAndDerivative(regionForThread);
    holder->measure_per_thread[threadId] += local_metric;

  }


};




} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDemonsImageToImageMetric.txx"
#endif

#endif
