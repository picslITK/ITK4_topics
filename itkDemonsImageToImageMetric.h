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

#include "itkMultiThreader.h"

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

  /** Image-accessor typedefs */
  typedef TFixedImage   FixedImageType;
  typedef TMovingImage  MovingImageType;
  typedef typename FixedImageType::Pointer FixedImagePointerType;
  typedef typename MovingImageType::Pointer MovingImagePointerType;
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
  itkSetMacro(FixedImage, FixedImagePointerType);
  itkGetMacro(FixedImage, FixedImagePointerType);
  itkSetMacro(MovingImage, MovingImagePointerType);
  itkGetMacro(MovingImage, MovingImagePointerType);


  /** The basic operations required for local and global metric  */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

  virtual unsigned int GetNumberOfParameters() const { return 0; }

  virtual MeasureType GetValue( const ParametersType & parameters ) const { return 0; }

  virtual void GetDerivative(const ParametersType &,
      DerivativeType & derivative) const {}

  /** This function computes the local voxel-wise contribution of 
   *  the metric to the global integral of the metric/derivative.
   */
  void ComputeLocalContributionToMetricAndDerivative(PointType fixedPoint, PointType movingPoint) 
  {
    /** Only the voxelwise contribution given the point pairs. */
  }

  /** This function is calls the ComputeMetricAndDerivative() function
   *  over the domain of interest.  
   */
  void ComputeMetricAndDerivative() 
  {
    /** For each location in the virtual domain, map to both the fixed and moving space
     *  and compute the values of the voxels in the corresponding locations.  There should 
     *  be a transform between the virtual space and the fixed/moving space s.t. the images 
     *  are interpolated in an unbiased manner. 
     */
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
  itkSetMacro(VirtualDomainOrigin, OriginType);
  itkGetMacro(VirtualDomainOrigin, OriginType);
  itkSetMacro(VirtualDomainDirection, DirectionType);
  itkGetMacro(VirtualDomainDirection, DirectionType);


protected:

  DemonsImageToImageMetric();
  virtual ~DemonsImageToImageMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  //purposely not implemented
  DemonsImageToImageMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

  FixedImagePointerType m_FixedImage;
  MovingImagePointerType m_MovingImage;
  MeasureType*     m_ThreaderMSE;
  DerivativeType*  m_ThreaderDerivatives;
  double           m_Normalizer;

  unsigned int m_NumberOfThreads;
  unsigned int m_NumberOfParameters;

  SizeType  m_VirtualDomainSize;
  SpacingType  m_VirtualDomainSpacing;
  OriginType  m_VirtualDomainOrigin;
  DirectionType  m_VirtualDomainDirection;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDemonsImageToImageMetric.txx"
#endif

#endif
