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
#include "itkMultiThreader.h"
namespace itk
{
/** \class ObjectToObjectMetric
 * \brief Computes similarity between regions of two objects.
 *
 * This Class is templated over the type of the two input objects.
 * This is the base class for a hierarchy of similarity metrics that may,
 * in derived classes, operate on meshes, images, etc.
 * This class computes a value that measures the similarity
 * between the two objects.
 *
 * We expect the implementation to be multi-threaded.
 *
 *
 * \ingroup RegistrationMetrics
 *
 */

template< class TFixedObject,  class TMovingObject >
class ITK_EXPORT ObjectToObjectMetric:
  public SingleValuedCostFunction
{
public:
  /** Standard class typedefs. */
  typedef ObjectToObjectMetric         Self;
  typedef SingleValuedCostFunction   Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Type used for representing point components  */
  typedef typename Superclass::ParametersValueType CoordinateRepresentationType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ObjectToObjectMetric, SingleValuedCostFunction);

  /**  Type of the measure. */
  typedef typename Superclass::MeasureType MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType DerivativeType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Set/Get number of threads to use for computations. */
  void SetNumberOfThreads(unsigned int numberOfThreads);

  itkGetConstReferenceMacro(NumberOfThreads, unsigned int);

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void)
  throw ( ExceptionObject );

  /** Initialize the components related to supporting multiple threads */
  virtual void MultiThreadingInitialize(void)
  throw ( ExceptionObject );

  typedef MultiThreader MultiThreaderType;

  /** Get the Threader. */
  itkGetConstObjectMacro(Threader, MultiThreaderType);

protected:
  ObjectToObjectMetric();
  virtual ~ObjectToObjectMetric();

  void PrintSelf(std::ostream & os, Indent indent) const;

  unsigned long          m_NumberOfParameters;
  mutable ParametersType m_Parameters;

  struct MultiThreaderParameterType {
    ObjectToObjectMetric *metric;
  };

  MultiThreaderType::Pointer m_Threader;
  MultiThreaderParameterType m_ThreaderParameter;
  bool                       m_WithinThreadPreProcess;
  bool                       m_WithinThreadPostProcess;

  unsigned int m_NumberOfThreads;

private:
  ObjectToObjectMetric(const Self &); //purposely not implemented
  void operator=(const Self &);     //purposely not implemented


};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkObjectToObjectMetric.txx"
#endif

#endif
