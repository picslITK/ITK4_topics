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
 * A \c ScalarScale value can be set instead of an array of parameter
 * scales. A single scalar scale value is useful for memory-efficient
 * uniform scaling of a dense transform (e.g. DisplacementFieldTransform).
 * Note that \c SetUseScalarScale must be called to enable this option.
 *
 * \c SetScales allows setting of a per-parameter scaling array.
 *
 * Threading of some optimizer operations may be handled within
 * derived classes, for example in GradientDescentOptimizer.
 *
 * Derived classes must override StartOptimization, and then call
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
      {
      itkDebugMacro("setting scales to " <<  scales);
      m_Scales = scales;
      this->Modified();
      }
  }

  /** Get current parameters scaling. */
  itkGetConstReferenceMacro(Scales, ScalesType);

  /** Set Scalar scaling value. Does NOT set the m_UseScalarScaling flag
   * to true. */
  itkSetMacro( ScalarScale, InternalComputationValueType );

  /** Get the scalar scaling value. */
  itkGetConstReferenceMacro( ScalarScale, InternalComputationValueType );

  /** Accessors for UseScalarScale flag */
  itkSetMacro(UseScalarScale, bool);
  itkGetConstReferenceMacro(UseScalarScale, bool);
  itkBooleanMacro(UseScalarScale);

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

  /** Get current position of the optimization. */
  itkGetConstReferenceMacro(CurrentPosition, ParametersType);

  /** Run the optimization.
   * \note Derived classes must override and call this supercall method, then
   * perform any additional initialization before performing optimization. */
  virtual void StartOptimization()
  {
    /* Validate some settings */
    if( this->m_Metric.IsNull() )
      {
      itkExceptionMacro("m_Metric must be set.");
      return;
      }

    /* Initialize or validate scales, if not using a scalar scale value. */
    if( ! m_UseScalarScale )
      {
      if ( m_Scales.Size() == 0 )
        {
        ScalesType scales( this->m_Metric->GetNumberOfParameters() );
        scales.Fill(1.0f);
        SetScales(scales);
        }
      else
        {
        if( this->m_Scales.Size() != this->m_Metric->GetNumberOfParameters() )
          {
          itkExceptionMacro("Size of scales (" << this->m_Scales.Size()
                            << ") must match size of parameters (" <<
                            this->m_Metric->GetNumberOfParameters() << ").");
          }
        }
      }
  }
protected:

  /** Default constructor */
  ObjectToObjectOptimizerBase()
  {
    this->m_Metric = NULL;
    this->m_Value = 0;
    this->m_ScalarScale = 1.0;
    this->m_UseScalarScale = false;
  }
  virtual ~ObjectToObjectOptimizerBase(){}

  MetricTypePointer             m_Metric;
  ThreadIdType                  m_NumberOfThreads;

  /** Metric measure value at a given iteration */
  MeasureType                   m_Value;

  /** Scales. Gets set to size of metric parameters and filled
   * with 1.0 in StartOptimization, if has not yet otherwise
   * been set, and not using a scalar scaling factor. */
  ScalesType                    m_Scales;

  /** A scalar scale value to use instead of an array of scales.
   * Useful for dense transforms. Defaults to 1.0 */
  InternalComputationValueType  m_ScalarScale;

  /** Flag to indicate if a scalar scale value is being used instead
   * of a scales array. Defaults to false. */
  bool                          m_UseScalarScale;

  // Keep m_CurrentPosition as a protected var so that subclasses can
  // have fast access.  This is important when optimizing high-dimensional
  // spaces, e.g. bspline transforms.
  ParametersType                m_CurrentPosition;

private:

  //purposely not implemented
  ObjectToObjectOptimizerBase( const Self & );
  void operator=( const Self& );      //purposely not implemented

};

} // end namespace itk

#endif
