/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkObjectToObjectOptimizerBase.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkObjectToObjectOptimizerBase_h
#define __itkObjectToObjectOptimizerBase_h

#include "itkTransformParameters.h"
#include "itkObjectToObjectMetric.h"
#include "itkIntTypes.h"

namespace itk
{
/** \class ObjectToObjectOptimizerBase
 * \brief Abstract base for object-to-object metric optimizers.
 *
 * Threading of some optimizer operations may be handled within
 * derived classes, for example in GradientDescentOptimizer.
 *
 * Derived classes must override StartOptimization, which is called
 * to initialize and run the optimization.
 *
 * \ingroup ITK-Optimizers
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
  itkTypeMacro(Self, Superclass);

  /**  Scale type. */
  typedef TransformParameters< double >             ScalesType;

  /** Metric function type */
  typedef ObjectToObjectMetric                      MetricType;
  typedef MetricType::Pointer                       MetricTypePointer;
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
  void SetScales(const ScalesType & scales)
  {
    if( scales != m_Scales )
    itkDebugMacro("setting scales to " <<  scales);
    m_Scales = scales;
    this->Modified();
  }

  /** Get current parameters scaling. */
  itkGetConstReferenceMacro(Scales, ScalesType);

  /** Set the number of threads to use when threading. */
  virtual void SetNumberOfThreads( ThreadIdType number )
  {
    if( number < 1 )
      {
      itkExceptionMacro("Number of threads must be > 0");
      }
    if( number != this->m_NumberOfThreads )
      {
      this->m_NumberOfThreads = number;
      this->Modified();
      }
  }

  /** Run the optimization */
  virtual void StartOptimization() = 0;

protected:

  /** Default constructor */
  ObjectToObjectOptimizerBase()
  {
    this->m_Metric = NULL;
    this->m_Value = 0;
  }
  virtual ~ObjectToObjectOptimizerBase(){}

  MetricTypePointer             m_Metric;
  ThreadIdType                  m_NumberOfThreads;

  /** Metric measure value at a given iteration */
  MeasureType                   m_Value;

  /** Scales */
  ScalesType                    m_Scales;

private:

  //purposely not implemented
  ObjectToObjectOptimizerBase( const Self & );
  void operator=( const Self& );      //purposely not implemented

};

} // end namespace itk

#endif
