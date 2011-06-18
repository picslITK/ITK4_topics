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
 * \ingroup RegistrationMetrics
 *
 * \ingroup ITK-Review
 */

template< class TFixedObject,  class TMovingObject >
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

  /**  Type of the measure. */
  typedef typename Superclass::MeasureType    MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType DerivativeType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  /** Type of coordinate system used to calculate values, derivatives */
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

  /** This method returns the derivative of the cost function */
  virtual void GetDerivative(DerivativeType & derivative)
  {
    MeasureType value;
    this->GetValueAndDerivative(value, derivative);
  }

  /** This method returns the value and derivative of the cost function */
  virtual void GetValueAndDerivative(MeasureType & value,
                                     DerivativeType & derivative) = 0;

protected:
  ObjectToObjectMetric();
  virtual ~ObjectToObjectMetric();

  void PrintSelf(std::ostream & os, Indent indent) const;

  /* Necessary ?? */
  mutable ParametersType      m_Parameters;

  DerivativeSourceType       m_DerivativeSource;

private:
  ObjectToObjectMetric(const Self &); //purposely not implemented
  void operator=(const Self &);     //purposely not implemented

  /** Provide these two methods to satisfy pure virtuals within
   * SingleValuedCostFunction. This is a sign that we probalby shouldn't
   * be deriving this from SingleValuedCostFunction. */
  MeasureType GetValue( const ParametersType& ) const
    {
    itkExceptionMacro("Not implemented. Use GetValue(void).");
    };

  /** This method returns the derivative of the cost function */
  void GetDerivative( const ParametersType &, DerivativeType &) const
    {
    itkExceptionMacro("Not implemented. Use GetDerivative(DerivativeType&).");
    };

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkObjectToObjectMetric.txx"
#endif

#endif
