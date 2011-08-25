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
#ifndef __itkTimeVaryingVelocityFieldIntegrationImageFilter_h
#define __itkTimeVaryingVelocityFieldIntegrationImageFilter_h

#include "itkImageToImageFilter.h"

//#include "itkVectorInterpolateImageFunction.h"

namespace itk
{

/* Forward-declaration to avoid including the header file that would
 * introduce a circular dependency in the Transform module. */
template< class TInputImage, class TCoordRep >
class VectorInterpolateImageFunction;

/**
 * \class TimeVaryingVelocityFieldIntegrationImageFilter
 * \brief Integrate a time-varying velocity field using 4th order Runge-Kutta.
 *
 * Diffeomorphisms are topology-preserving mappings that are useful for
 * describing biologically plausible deformations.  Mathematically, a
 * diffeomorphism, \phi, is generated from a time-varying velocity field, v, as
 * described by the integral equation:
 *
 * \phi(t_b) = \phi(t_a) + \int_t_a^t_b v(\phi(t),t) dt
 *
 * In this class, the input is the time-varying velocity field and an initial
 * diffeomorophism.  The output diffeomorphism is produced using fourth order
 * Runge-Kutta.
 *
 * \warning The output deformation field needs to have dimensionality of 1
 * less than the input time-varying velocity field.
 *
 * \author Nick Tustison
 * \author Brian Avants
 *
 * \ingroup ITKTransform
 *
 * \wiki
 * \wikiexample{}
 * \endwiki
 */
template<class TTimeVaryingVelocityField, class TDeformationField =
 Image<typename TTimeVaryingVelocityField::PixelType,
 GetImageDimension<TTimeVaryingVelocityField>::ImageDimension - 1> >
class ITK_EXPORT TimeVaryingVelocityFieldIntegrationImageFilter :
  public ImageToImageFilter<TTimeVaryingVelocityField, TDeformationField>
{
public:
  typedef TimeVaryingVelocityFieldIntegrationImageFilter  Self;
  typedef ImageToImageFilter
    <TTimeVaryingVelocityField, TDeformationField>        Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( TimeVaryingVelocityFieldIntegrationImageFilter,
    ImageToImageFilter );

  /**
   * Dimensionality of input data is assumed to be one more than the output
   * data the same. */
  itkStaticConstMacro( InputImageDimension, unsigned int,
    TTimeVaryingVelocityField::ImageDimension );

  itkStaticConstMacro( OutputImageDimension, unsigned int,
    TDeformationField::ImageDimension );

  typedef TTimeVaryingVelocityField                 TimeVaryingVelocityFieldType;
  typedef TDeformationField                         DeformationFieldType;
  typedef typename DeformationFieldType::Pointer    DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType  VectorType;
  typedef typename VectorType::RealValueType        RealType;
  typedef typename DeformationFieldType::PointType  PointType;
  typedef typename DeformationFieldType::RegionType OutputRegionType;

  typedef VectorInterpolateImageFunction
    <TimeVaryingVelocityFieldType, RealType>    VelocityFieldInterpolatorType;
  typedef typename VelocityFieldInterpolatorType::Pointer
                                                VelocityFieldInterpolatorPointer;

  typedef VectorInterpolateImageFunction
    <DeformationFieldType, RealType>        DeformationFieldInterpolatorType;
  typedef typename DeformationFieldInterpolatorType::Pointer
                                            DeformationFieldInterpolatorPointer;

  /** Set the time-varying velocity field interpolator.  Default = linear. */
  itkSetObjectMacro( VelocityFieldInterpolator, VelocityFieldInterpolatorType );

  /** Get the time-varying velocity field interpolator.  Default = linear. */
  itkGetObjectMacro( VelocityFieldInterpolator, VelocityFieldInterpolatorType );

  /**
   * Set the deformation field interpolator for the initial diffeomorphism
   * (if set).  Default = linear.
   */
  itkSetObjectMacro( DeformationFieldInterpolator,
    DeformationFieldInterpolatorType );

  /**
   * Get the deformation field interpolator for the initial diffeomorphism
   * (if set).  Default = linear.
   */
  itkGetObjectMacro( DeformationFieldInterpolator,
    DeformationFieldInterpolatorType );

  /**
   * Set the initial diffeomorphism
   */
  itkSetObjectMacro( InitialDiffeomorphism, DeformationFieldType );

  /**
   * Get the initial diffeomorphism
   */
  itkGetObjectMacro( InitialDiffeomorphism, DeformationFieldType );

  /**
   * Set the lower time bound defining the integration domain of the transform.
   * We assume that the total possible time domain is [0,1].
   */
  itkSetClampMacro( LowerTimeBound, RealType, 0, 1 );

  /**
   * Get the lower time bound defining the integration domain of the transform.
   * We assume that the total possible time domain is [0,1].
   */
  itkGetConstMacro( LowerTimeBound, RealType );

  /**
   * Set the upper time bound defining the integration domain of the transform.
   * We assume that the total possible time domain is [0,1].
   */
  itkSetClampMacro( UpperTimeBound, RealType, 0, 1 );

  /**
   * Get the upper time bound defining the integration domain of the transform.
   * We assume that the total possible time domain is [0,1].
   */
  itkGetConstMacro( UpperTimeBound, RealType );

  /**
   * Set the number of integration steps used in the Runge-Kutta solution of the
   * initial value problem.  Default = 10.
   */
  itkSetMacro( NumberOfIntegrationSteps, unsigned int );

  /**
   * Get the number of integration steps used in the Runge-Kutta solution of the
   * initial value problem.  Default = 10.
   */
  itkGetConstMacro( NumberOfIntegrationSteps, unsigned int );

protected:
  TimeVaryingVelocityFieldIntegrationImageFilter();
  ~TimeVaryingVelocityFieldIntegrationImageFilter();

  void PrintSelf( std::ostream & os, Indent indent ) const;

  virtual void GenerateOutputInformation();

  virtual void BeforeThreadedGenerateData();

  virtual void ThreadedGenerateData( const OutputRegionType &, ThreadIdType );

private:
  TimeVaryingVelocityFieldIntegrationImageFilter( const Self & );
                                          //purposely not implemented
  void operator=( const Self & );         //purposely not implemented

  VectorType IntegrateVelocityAtPoint( const PointType );

  VelocityFieldInterpolatorPointer          m_VelocityFieldInterpolator;
  RealType                                  m_LowerTimeBound;
  RealType                                  m_UpperTimeBound;
  unsigned int                              m_NumberOfIntegrationSteps;

  DeformationFieldPointer                   m_InitialDiffeomorphism;
  DeformationFieldInterpolatorPointer       m_DeformationFieldInterpolator;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTimeVaryingVelocityFieldIntegrationImageFilter.hxx"
#endif

#endif
