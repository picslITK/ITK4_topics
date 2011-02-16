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
#include "itkMultiThreader.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
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

  /** The basic operations required for local and global metric  */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

  virtual unsigned int GetNumberOfParameters() const { return 0; }

  virtual MeasureType GetValue( const ParametersType & parameters ) const { return 0; }

  virtual void GetDerivative(const ParametersType &,
      DerivativeType & derivative) const {}

  /** This function computes the local voxel-wise contribution of 
   *  the metric to the global integral of the metric/derivative.
   */
  double ComputeLocalContributionToMetricAndDerivative(PointType mappedFixedPoint, PointType mappedMovingPoint) 
  {
    /** Only the voxelwise contribution given the point pairs. */
    double fpix=this->m_FixedInterpolator->Evaluate(mappedFixedPoint); 
    double mpix=this->m_MovingInterpolator->Evaluate(mappedMovingPoint);
    return fabs((double)fpix-(double)mpix);
  }

  /** This function is calls the ComputeMetricAndDerivative() function
   *  over the domain of interest.  
   */
  double ComputeMetricAndDerivative() 
  {
    /** For each location in the virtual domain, map to both the fixed and moving space
     *  and compute the values of the voxels in the corresponding locations.  There should 
     *  be a transform between the virtual space and the fixed/moving space s.t. the images 
     *  are interpolated in an unbiased manner. 
     */
    double metric_sum=0;
    unsigned long ct=0;
    if ( ! this->m_VirtualImage ) {
      // std::cout <<" allocate-a " << std::endl;
      RegionType region;
      // std::cout <<" allocate-size " << this->GetVirtualDomainSize() << std::endl;
      region.SetSize(this->GetVirtualDomainSize() );
      region.SetIndex(this->GetVirtualDomainIndex() );
      // std::cout <<" allocate-c " << std::endl;
      this->m_VirtualImage = FixedImageType::New();
      // std::cout <<" allocate-d " << std::endl;
      this->m_VirtualImage->SetSpacing( this->GetVirtualDomainSpacing() );
      // std::cout <<" allocate-e " << std::endl;
      this->m_VirtualImage->SetOrigin( this->GetVirtualDomainOrigin() );
      // std::cout <<" allocate-f " << std::endl;
      this->m_VirtualImage->SetDirection( this->GetVirtualDomainDirection() );
      // std::cout <<" allocate-g " << std::endl;
      this->m_VirtualImage->SetRegions( region );
      // std::cout <<" allocate-h " << std::endl;
      this->m_VirtualImage->Allocate();
      // std::cout <<" allocate-i " << std::endl;
      this->m_VirtualImage->FillBuffer( 0 );
      // std::cout <<" allocate-j " << std::endl;
      this->m_FixedInterpolator=FixedInterpolatorType::New();
      // std::cout <<" allocate-k " << std::endl;
      this->m_MovingInterpolator=MovingInterpolatorType::New();
      // std::cout <<" allocate-l " << std::endl;
      this->m_FixedInterpolator->SetInputImage(m_FixedImage);
      // std::cout <<" allocate-m " << std::endl;
      this->m_MovingInterpolator->SetInputImage(m_MovingImage);
      // std::cout <<" allocate-n " << std::endl;
    }
    // std::cout <<" allocate-done " << std::endl;

    ImageRegionIteratorWithIndex<FixedImageType> ItV( this->m_VirtualImage,
        this->m_VirtualImage->GetRequestedRegion() );
    ItV.GoToBegin();
    std::cout << "Iterate from Ind " <<  ItV.GetIndex() << std::endl;
    while( !ItV.IsAtEnd() )
    {
      /** use the fixed and moving transforms to compute the corresponding points.*/
      bool sampleOk = true;
      // convert the index to a point 
      PointType mappedPoint;
      PointType mappedFixedPoint;
      PointType mappedMovingPoint;
      //      std::cout << "Tx-1 " <<std::endl;
      this->m_VirtualImage->TransformIndexToPhysicalPoint(ItV.GetIndex(),mappedPoint);

      //      std::cout << "Tx-2 " <<std::endl;
      // Use generic transform to compute mapped position
      if ( this->m_FixedImageTransform ) {
      mappedFixedPoint = this->m_FixedImageTransform->TransformPoint(mappedPoint);
      mappedMovingPoint = this->m_MovingImageTransform->TransformPoint(mappedPoint);
      }
      //      std::cout << "Tx-3 " <<std::endl;
      if ( !this->m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint) ||  
           !this->m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint) ) 
	sampleOk=false;
      if ( sampleOk )
        {
	  //      std::cout << "Samp-1 " <<std::endl;
	  double metricval=this->ComputeLocalContributionToMetricAndDerivative(mappedFixedPoint,mappedMovingPoint);
	  //      std::cout << "Samp-2 " << metricval << " pt " <<mappedFixedPoint << std::endl;
	  //	  std::cout <<" ct " << ct << " mv " << metricval <<  " ind " << ItV.GetIndex() << " buf-check " << this->m_MovingInterpolator->IsInsideBuffer(mappedPoint) <<  " point " << mappedPoint << std::endl;
          metric_sum+=metricval;	  
	  //      std::cout << "Samp-3 " <<std::endl;
	  ct++;
        }
      ++ItV;
      if ( ItV.IsAtEnd() ){
	--ItV;
        std::cout << " ItV end " << ItV.GetIndex() << std::endl;
        ++ItV;
      }
      //  std::cout << "It-2 " <<std::endl;
    }
    if ( ct > 0 ) {
      std::cout << " metric_sum " << metric_sum << " ct " << ct <<std::endl;
      return metric_sum/(double)(ct);
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
    return energy/(InternalComputationValueType)measure_per_thread.size();
  }

  static void ComputeMetricValueInRegionOnTheFlyThreaded(const ImageRegionType &regionForThread, int threadId,  Self *holder){

    std::cout << regionForThread << std::endl;
    InternalComputationValueType local_metric;
    typedef itk::ConstNeighborhoodIterator<ImageType> IteratorType;

    typedef itk::ImageRegionIterator<DeformationFieldType> UpdateIteratorType;
    std::cout << " allocate metric " << std::endl;
    MetricTypePointer objectMetric = MetricType::New();
    std::cout << " allocate metric done " << std::endl;
    objectMetric->SetFixedImage(holder->fixed_image);
    objectMetric->SetMovingImage(holder->moving_image); 
    std::cout << " set t " << std::endl;
    objectMetric->SetFixedImageTransform(holder->transformF);
    objectMetric->SetMovingImageTransform(holder->transformM);
    std::cout << " set done " << std::endl;
 
    std::cout << " setv size " << regionForThread.GetSize() << " ind " << regionForThread.GetIndex() << std::endl;
    objectMetric->SetVirtualDomainSize(regionForThread.GetSize());
    objectMetric->SetVirtualDomainIndex(regionForThread.GetIndex());
    objectMetric->SetVirtualDomainSpacing(holder->fixed_image->GetSpacing());
    objectMetric->SetVirtualDomainOrigin(holder->fixed_image->GetOrigin());
    objectMetric->SetVirtualDomainDirection(holder->fixed_image->GetDirection());
    holder->measure_per_thread[threadId] = NumericTraits<InternalComputationValueType>::Zero;
    /** Compute one iteration of the metric */
    local_metric=objectMetric->ComputeMetricAndDerivative();
    holder->measure_per_thread[threadId] += local_metric;

  }


};




} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDemonsImageToImageMetric.txx"
#endif

#endif
