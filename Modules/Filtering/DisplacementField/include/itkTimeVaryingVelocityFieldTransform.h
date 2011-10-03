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
#ifndef __itkTimeVaryingVelocityFieldTransform_h
#define __itkTimeVaryingVelocityFieldTransform_h

#include "itkDisplacementFieldTransform.h"

#include "itkImageVectorTransformParametersHelper.h"

namespace itk
{

/** \class TimeVaryingVelocityFieldTransform
 * \brief Transform objects based on integration of a time-varying velocity
 * field.
 *
 * Diffeomorphisms are topology-preserving mappings that are useful for
 * describing biologically plausible deformations.  Mathematically, a
 * diffeomorphism, \phi, is generated from a time-varying velocity field, v, as
 * described by the first-order differential equation:
 *
 * v(\phi(x,t), t) = \frac{d\phi(x, t)}{dt}, \phi(x, 0) = x
 *
 * In this class, the input is the time-varying velocity field.  The output
 * diffeomorphism is produced using fourth order Runge-Kutta.
 *
 * \warning The output deformation field needs to have dimensionality of 1
 * less than the input time-varying velocity field. It is assumed that the
 * last dimension of the time-varying velocity field corresponds to Time,
 * and the other dimensions represent Space.
 *
 * \author Nick Tustison
 * \author Brian Avants
 *
 * \ingroup Transforms
 * \ingroup ITKDisplacementField
 *
 * \wiki
 * \wikiexample{}
 * \endwiki
 */
template<class TScalar, unsigned int NDimensions>
class ITK_EXPORT TimeVaryingVelocityFieldTransform :
  public DisplacementFieldTransform<TScalar, NDimensions>
{
public:
  /** Standard class typedefs. */
  typedef TimeVaryingVelocityFieldTransform                 Self;
  typedef DisplacementFieldTransform<TScalar, NDimensions>  Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( TimeVaryingVelocityFieldTransform, DisplacementFieldTransform );

  /** New macro for creation of through a Smart Pointer */
  itkSimpleNewMacro( Self );

  /** Create another transform of the same type. */
  virtual::itk::LightObject::Pointer CreateAnother(void) const;

  /** InverseTransform type. */
  typedef typename Superclass:: InverseTransformBasePointer InverseTransformBasePointer;

  /** Scalar type. */
  typedef typename Superclass::ScalarType          ScalarType;

  /** Type of the input parameters. */
  typedef  typename Superclass::ParametersType          ParametersType;
  typedef  typename Superclass::NumberOfParametersType  NumberOfParametersType;

  /** Jacobian type. */
  typedef typename Superclass::JacobianType        JacobianType;

  /** Standard coordinate point type for this class. */
  typedef typename Superclass::InputPointType      InputPointType;
  typedef typename Superclass::OutputPointType     OutputPointType;

  /** Standard vector type for this class. */
  typedef typename Superclass::InputVectorType     InputVectorType;
  typedef typename Superclass::OutputVectorType    OutputVectorType;

  /** Derivative type */
  typedef typename Superclass::DerivativeType       DerivativeType;

  /** Dimension of the domain spaces. */
  itkStaticConstMacro( Dimension, unsigned int, NDimensions );

  /** Dimension of the time varying velocity field. */
  itkStaticConstMacro( TimeVaryingVelocityFieldDimension, unsigned int, NDimensions+1 );

  /**
   * Define the time-varying velocity field type and corresponding interpolator
   * type.
   */
  typedef typename Superclass::DisplacementFieldType   DisplacementFieldType;
  typedef typename DisplacementFieldType::PixelType    VectorType;
  typedef typename DisplacementFieldType::PointType    PointType;

  typedef Image<OutputVectorType,
    itkGetStaticConstMacro( TimeVaryingVelocityFieldDimension )>
                                     TimeVaryingVelocityFieldType;
  typedef VectorInterpolateImageFunction<TimeVaryingVelocityFieldType,
    ScalarType>                      TimeVaryingVelocityFieldInterpolatorType;
  typedef typename TimeVaryingVelocityFieldInterpolatorType::Pointer
                                     TimeVaryingVelocityFieldInterpolatorPointer;

  /** Define the internal parameter helper used to access the field */
  typedef ImageVectorTransformParametersHelper
    <ScalarType, OutputVectorType::Dimension,
    itkGetStaticConstMacro( Dimension ) + 1>      TransformParametersHelperType;

  /** Get the time-varying deformation field. */
  itkGetObjectMacro( TimeVaryingVelocityField, TimeVaryingVelocityFieldType );

  /** Set the time-varying field.  */
  virtual void SetTimeVaryingVelocityField( TimeVaryingVelocityFieldType * );

  /** Set the interpolator for the time-varying velocity field. */
  itkSetObjectMacro( TimeVaryingVelocityFieldInterpolator,
    TimeVaryingVelocityFieldInterpolatorType );

  /** Get the interpolator for the time-varying velocity field. */
  itkGetConstObjectMacro( TimeVaryingVelocityFieldInterpolator,
    TimeVaryingVelocityFieldInterpolatorType );

  /** Get the modification time of deformation field */
  itkGetConstMacro( TimeVaryingVelocityFieldSetTime, unsigned long );

  /**
   * Set the deformation field. We want to override the base class
   * implementation since we don't want to optimize over the deformation
   * field for this class but rather the time-varying velocity field
   */
  itkSetObjectMacro( DisplacementField, DisplacementFieldType );

  /**
   * Set whether or not the time-varying velocity field should be integrated.
   * Default is true.  However, we don't want to integrated when unnecessary
   * so we allow the user to turn it off.
   */
  itkSetMacro( IntegrateTimeVaryingVelocityField, bool );

  /**
   * Get whether or not the time-varying velocity field should be integrated.
   * Default is true.  However, we don't want to integrated when unnecessary
   * so we allow the user to turn it off.
   */
  itkGetConstMacro( IntegrateTimeVaryingVelocityField, bool );

  /**
   * Set/Get whether or not the time-varying velocity field should be integrated.
   * Default is true.  However, we don't want to integrated when unnecessary
   * so we allow the user to turn it off.
   */
  itkBooleanMacro( IntegrateTimeVaryingVelocityField );

  /**
   * Set the transformation parameters. This sets the time-varying velocity
   * field image directly.
   */
  virtual void SetParameters(const ParametersType & params);

  /** Trigger the computation of the displacement field by integrating
   * the time-varying velocity field. */
  virtual void IntegrateVelocityField();

  /** Set the fixed parameters and update internal transformation. */
  virtual void SetFixedParameters( const ParametersType & )
    {
    itkExceptionMacro( "SetFixedParameters unimplemented." );
    }

  /** Get the Fixed Parameters. */
  virtual const ParametersType & GetFixedParameters() const
    {
    itkExceptionMacro( "GetFixedParameters unimplemented." );
    }

  virtual void UpdateTransformParameters( DerivativeType &,
    ScalarType factor = 1.0 );

  /** Return an inverse of this transform. */
  bool GetInverse( Self *inverse ) const;

  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const;

  /** This transform is not linear. */
  virtual bool IsLinear() const { return false; }

  /** Get the number of local parameters */
  NumberOfParametersType GetNumberOfLocalParameters() const
    {
    return Dimension;
    }

  /** Does the transform have local support */
  virtual bool HasLocalSupport() const
    {
    return true;
    }

  /**
   * Set the lower time bound defining the integration domain of the transform.
   * We assume that the total possible time domain is [0,1]
   */
  itkSetClampMacro( LowerTimeBound, ScalarType, 0, 1 );

  /**
   * Get the lower time bound defining the integration domain of the transform.
   * We assume that the total possible time domain is [0,1]
   */
  itkGetConstMacro( LowerTimeBound, ScalarType );

  /**
   * Set the upper time bound defining the integration domain of the transform.
   * We assume that the total possible time domain is [0,1]
   */
  itkSetClampMacro( UpperTimeBound, ScalarType, 0, 1 );

  /**
   * Get the upper time bound defining the integration domain of the transform.
   * We assume that the total possible time domain is [0,1]
   */
  itkGetConstMacro( UpperTimeBound, ScalarType );

  /**
   * Set the number of integration steps used in the Runge-Kutta solution of the
   * initial value problem.  Default = 10;
   */
  itkSetMacro( NumberOfIntegrationSteps, unsigned int );

  /**
   * Get the number of integration steps used in the Runge-Kutta solution of the
   * initial value problem.  Default = 10;
   */
  itkGetConstMacro( NumberOfIntegrationSteps, unsigned int );

protected:
  TimeVaryingVelocityFieldTransform();
  virtual ~TimeVaryingVelocityFieldTransform();
  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  TimeVaryingVelocityFieldTransform( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** The deformation field and its inverse (if it exists). */
  typename TimeVaryingVelocityFieldType::Pointer    m_TimeVaryingVelocityField;

  TimeVaryingVelocityFieldInterpolatorPointer
                                         m_TimeVaryingVelocityFieldInterpolator;

  ScalarType                                m_LowerTimeBound;
  ScalarType                                m_UpperTimeBound;

  unsigned int                              m_NumberOfIntegrationSteps;
  bool                                      m_IntegrateTimeVaryingVelocityField;
  unsigned long                             m_TimeVaryingVelocityFieldSetTime;
};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkTimeVaryingVelocityFieldTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkTimeVaryingVelocityFieldTransform.hxx"
#endif

#endif // __itkTimeVaryingVelocityFieldTransform_h
