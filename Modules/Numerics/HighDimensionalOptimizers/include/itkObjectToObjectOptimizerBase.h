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
#ifndef __itkObjectToObjectOptimizerBase_h
#define __itkObjectToObjectOptimizerBase_h

#include "itkTransformParameters.h"
#include "itkObjectToObjectMetric.h"
#include "itkIntTypes.h"

namespace itk
{
/** \class ObjectToObjectOptimizerBase
 * \brief Abstract base for object-to-object optimizers.
 *
 * The goal of this optimizer hierarchy is to work with metrics
 * of any type, i.e. working with any kind of object, such as
 * image or point-set.
 *
 * Transform parameters are not manipulated directly. Instead,
 * the optimizer retrieves the metric derivative from the metric,
 * modifies the derivative as required, then passes it back to
 * the metric as an update. The metric then processes it as
 * appropriate, typically by passing it to its transform that is
 * being optimized.
 *
 * \c SetScales allows setting of a per-local-parameter scaling array. If
 * unset, the \c m_Scales array will be initialized to all 1's.
 * Note that when used with transforms with local support, these scales
 * correspond to each _local_ parameter, and not to each parameter. For
 * example, in a DisplacementFieldTransform of dimensionality N, the Scales
 * is size N, with each element corresponding to a dimension within the
 * transform's displacement field, and is applied to each vector in the
 * displacement field.
 *
 * Threading of some optimizer operations may be handled within
 * derived classes, for example in GradientDescentOptimizer.
 *
 * \note Derived classes must override StartOptimization, and then call
 * this base class version to perform common initializations.
 *
 * \ingroup ITKHighDimensionalOptimizers
 */

class ITK_EXPORT ObjectToObjectOptimizerBase : public Object
{
public:
  /** Standard class typedefs. */
  typedef ObjectToObjectOptimizerBase                 Self;
  typedef Object                                      Superclass;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ObjectToObjectOptimizerBase, Object);

  /**  Scale type. */
  typedef TransformParameters< double >             ScalesType;

  /**  Parameters type. */
  typedef TransformParameters< double >             ParametersType;

  /** Metric function type */
  typedef ObjectToObjectMetric                      MetricType;
  typedef MetricType::Pointer                       MetricTypePointer;

  /** Number of parameters type */
  typedef MetricType::NumberOfParametersType        NumberOfParametersType;

  /** Measure type */
  typedef MetricType::MeasureType                   MeasureType;

  /** Internal computation value type */
  typedef MetricType::InternalComputationValueType
                                                InternalComputationValueType;

  /** Accessors for Metric */
  itkGetObjectMacro( Metric, MetricType );
  itkSetObjectMacro( Metric, MetricType );

  /** Accessor for metric value */
  itkGetConstReferenceMacro( Value, MeasureType );

  /** Set current parameters scaling. */
  itkSetMacro( Scales, ScalesType );

  /** Get current parameters scaling. */
  itkGetConstReferenceMacro( Scales, ScalesType );

  /** Set the number of threads to use when threading. */
  virtual void SetNumberOfThreads( ThreadIdType number );

  /** Get a reference to the current position of the optimization.
   * This returns the parameters from the assigned metric, since the optimizer
   * itself does not store a position. */
  const ParametersType & GetCurrentPosition();

  /** Run the optimization.
   * \note Derived classes must override and call this superclass method, then
   * perform any additional initialization before performing optimization. */
  virtual void StartOptimization();

protected:

  /** Default constructor */
  ObjectToObjectOptimizerBase();
  virtual ~ObjectToObjectOptimizerBase();

  MetricTypePointer             m_Metric;
  ThreadIdType                  m_NumberOfThreads;

  /** Metric measure value at a given iteration */
  MeasureType                   m_Value;

  /** Scales. Size is expected to be == metric->GetNumberOfLocalParameters().
   * See the main documentation for more details. */
  ScalesType                    m_Scales;


  /** Flag to avoid unnecessary arithmetic when scales are identity. */
  bool                          m_ScalesAreIdentity;

  virtual void PrintSelf(std::ostream & os, Indent indent) const;

private:

  //purposely not implemented
  ObjectToObjectOptimizerBase( const Self & );
  //purposely not implemented
  void operator=( const Self& );

};

} // end namespace itk

#endif
