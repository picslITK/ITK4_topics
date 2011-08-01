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
#ifndef __itkObjectToObjectMetric_h
#define __itkObjectToObjectMetric_h

#include "itkSingleValuedCostFunction.h"

namespace itk
{

/** \class ObjectToObjectMetric
 * \brief Computes similarity between regions of two objects.
 *
 * This Class is templated over the type of the two input objects.
 * This is the base class for a hierarchy of similarity metrics that may, in
 * derived classes, operate on meshes, images, etc.  This class computes a
 * value that measures the similarity between the two objects.
 *
 * Transform Optimization
 * This hierarchy currently assumes only the moving transform is 'active',
 * i.e. begin optimized when used in an optimizer. The goal however is
 * to allow for either moving, fixed or both transforms to be optimized
 * within a single metric.
 *
 * Derived classes must provide implementations for:
 *  GetValue
 *  GetValueAndDerivative
 *  Initialize
 *  GetNumberOfParameters
 *  GetNumberOfLocalParameters
 *  GetParameters
 *  HasLocalSupport
 *  UpdateTransformParameters
 *
 * \ingroup ITKRegistrationRefactoring
 */

class ITK_EXPORT ObjectToObjectMetric:
  public SingleValuedCostFunction
{
public:
  /** Standard class typedefs. */
  typedef ObjectToObjectMetric       Self;
  typedef SingleValuedCostFunction   Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ObjectToObjectMetric, SingleValuedCostFunction);

  /** Type used for representing object components  */
  typedef Superclass::ParametersValueType CoordinateRepresentationType;

  /** Type for internal computations */
  typedef double                          InternalComputationValueType;

  /**  Type of the measure. */
  typedef  Superclass::MeasureType        MeasureType;

  /**  Type of the derivative. */
  typedef  Superclass::DerivativeType     DerivativeType;

  /**  Type of the parameters. */
  typedef  Superclass::ParametersType       ParametersType;
  typedef  Superclass::ParametersValueType  ParametersValueType;

  /** Source of the object derivatives (image derivatives, in the case of
   * image to image metrics). Defaults to Moving. */
  typedef enum  { Fixed=0, Moving, Both } DerivativeSourceType;

  /**
   * Set source of derivative.  This variable allows the user to switch
   * between calculating the derivative with respect to the fixed
   * object or moving object.
   */
  itkSetMacro( DerivativeSource, DerivativeSourceType );

  /**
   * Get coordinate system type.
   */
  itkGetConstMacro( DerivativeSource, DerivativeSourceType );

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject ) = 0;

  /** This method returns the value of the cost function */
  virtual MeasureType GetValue() = 0;

  /** This method returns the derivative of the cost function.
   * \c derivative will be sized and allocated as needed by metric.
   * If it's already allocated at proper size, no new allocation is done. */
  virtual void GetDerivative(DerivativeType & derivative)
  {
    MeasureType value;
    this->GetValueAndDerivative(value, derivative);
  }

  /** This method returns the value and derivative of the cost function.
   * \c derivative will be sized and allocated as needed by metric.
   * If it's already proper size, no new allocation is done. */
  virtual void GetValueAndDerivative(MeasureType & value,
                                     DerivativeType & derivative) = 0;

  /** Methods for working with the metric's 'active' transform, e.g. the
   * transform being optimized in the case of registration. Some of these are
   * used in non-metric classes, e.g. optimizers. */
  virtual unsigned int GetNumberOfParameters() const = 0;
  virtual unsigned int GetNumberOfLocalParameters() const = 0;

  /** Get a const reference to the active transform's parameters */
  virtual const ParametersType & GetParameters() const = 0;

  /** Return whether the metric's active transform has local support,
   * i.e. is dense. */
  virtual bool HasLocalSupport() const = 0;

  /** Update the parameters of the metric's active transform.
   * Typically this call is passed through directly to the transform.
   * \c factor is a scalar multiplier for each value in update, and
   * defaults to 1.0 .
   * \c derivative must be the proper size, as retrieved
   * from GetNumberOfParameters. */
  virtual void UpdateTransformParameters( DerivativeType & derivative,
                                          ParametersValueType factor = 1.0) = 0;

protected:
  ObjectToObjectMetric() {}
  virtual ~ObjectToObjectMetric() {}

  void PrintSelf(std::ostream & os, Indent indent) const
  { Superclass::PrintSelf(os, indent); os << indent << "TODO..."; }

  DerivativeSourceType       m_DerivativeSource;

private:
  ObjectToObjectMetric(const Self &); //purposely not implemented
  void operator=(const Self &);     //purposely not implemented

  /** Provide these two methods to satisfy pure virtuals within
   * SingleValuedCostFunction. This is a sign that we probalby shouldn't
   * be deriving this from SingleValuedCostFunction. */
  MeasureType GetValue( const ParametersType& ) const
  { itkExceptionMacro("Not implemented. Use GetValue(void)."); }

  /** This method returns the derivative of the cost function */
  void GetDerivative( const ParametersType &, DerivativeType &) const
  { itkExceptionMacro("Not implemented. Use GetDerivative(DerivativeType&).");}

  void GetValueAndDerivative (const ParametersType &parameters,
                              MeasureType &value,
                              DerivativeType &derivative) const
  { itkExceptionMacro("Not implemented. Use GetValueAndDerivative( "
                      "MeasureType & value, DerivativeType & derivative)."); }

};
} // end namespace itk

#endif
