/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGradientDescentObjectOptimizerBase.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGradientDescentObjectOptimizerBase_h
#define __itkGradientDescentObjectOptimizerBase_h

#include "itkPoint.h"
#include "itkIndex.h"
#include "itkObjectToObjectOptimizerBase.h"
#include "itkArray1DToData.h"

namespace itk
{
/** \class GradientDescentObjectOptimizerBase
 *  \brief Base class for gradient descent optimizers.
 *
 *  Gradient modification is threaded in \c ModifyGradient.
 *  Derived classes must override \c ModifyGradientOverSubRange.
 *
 * \ingroup ITKRegistrationRefactoring
 */

class ITK_EXPORT GradientDescentObjectOptimizerBase
  : public ObjectToObjectOptimizerBase
{
public:
  /** Standard class typedefs. */
  typedef GradientDescentObjectOptimizerBase     Self;
  typedef ObjectToObjectOptimizerBase            Superclass;
  typedef SmartPointer< Self >                   Pointer;
  typedef SmartPointer< const Self >             ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GradientDescentObjectOptimizerBase, ObjectToObjectOptimizerBase);

  /** Codes of stopping conditions. */
  typedef enum {
    MaximumNumberOfIterations,
    MetricError,
    UpdateTransformParametersError,
    StepTooSmall,
    QuasiNewtonStepError,
    OtherError
    } StopConditionType;

  /** Metric type over which this class is templated */
  typedef Superclass::MetricType                    MetricType;
  typedef MetricType::Pointer                       MetricTypePointer;

  /** Derivative type */
  typedef MetricType::DerivativeType                DerivativeType;

  /** Measure type */
  typedef Superclass::MeasureType                   MeasureType;

  /** Internal computation type, for maintaining a desired precision */
  typedef Superclass::InternalComputationValueType
                                                InternalComputationValueType;

  /** Threader for gradient update */
  typedef Array1DToData<Self>                       ModifyGradientThreaderType;
  typedef ModifyGradientThreaderType::IndexRangeType  IndexRangeType;

  /** Get the most recent gradient values. */
  itkGetConstReferenceMacro( Gradient, DerivativeType );

  /** Get stop condition enum */
  itkGetConstReferenceMacro(StopCondition, StopConditionType);

  /** Set the number of iterations. */
  itkSetMacro(NumberOfIterations, SizeValueType);

  /** Get the number of iterations. */
  itkGetConstReferenceMacro(NumberOfIterations, SizeValueType);

  /** Get the current iteration number. */
  itkGetConstMacro(CurrentIteration, SizeValueType);

  /** Resume optimization.
   * This runs the optimization loop, and allows continuation
   * of stopped optimization */
  virtual void ResumeOptimization() = 0;

  /** Stop optimization. The object is left in a state so the
   * optimization can be resumed by calling ResumeOptimization. */
  virtual void StopOptimization(void);

  /** Get the reason for termination */
  virtual const std::string GetStopConditionDescription() const;

protected:

  /** Default constructor */
  GradientDescentObjectOptimizerBase();
  virtual ~GradientDescentObjectOptimizerBase(){}

  /** Modify the gradient in place, to advance the optimization.
   * This call performs a threaded modification for transforms with
   * local support (assumed to be dense). Otherwise the modification
   * is performed w/out threading.
   * The work is done in ModifyGradientOverSubRange in both cases.
   * At completion, m_Gradient can be used to update the transform
   * parameters. Derived classes may hold additional results in
   * other member variables.
   */
  virtual void ModifyGradient();

  /** Derived classes define this worker method to modify the gradient.
   * Modifications must be performed over the index range defined in
   * \c subrange.
   * Called from ModifyGradient(), either directly or via threaded
   * operation. */
  virtual void ModifyGradientOverSubRange( const IndexRangeType& subrange ) = 0;

  /* Common variables for optimization control and reporting */
  bool               m_Stop;
  StopConditionType  m_StopCondition;
  std::ostringstream m_StopConditionDescription;
  SizeValueType      m_NumberOfIterations;
  SizeValueType      m_CurrentIteration;

  /** Current gradient */
  DerivativeType     m_Gradient;

private:

  //purposely not implemented
  GradientDescentObjectOptimizerBase( const Self & );
  void operator=( const Self& );      //purposely not implemented

  /** Callback used during threaded gradient modification.
   * Gets assigned to the modify-gradient threader's
   * ThreadedGenerateData. Must be static so it can be used as a callback.
   * It simply calls ModifyGradientOverSubRange, which is where
   * derived classes should put their class-specific modification code.
   * An instance of this optimizer class is referenced through
   * \c holder, which is passed in via the threader's user data. */
  static void ModifyGradientThreaded(
                                  const IndexRangeType& rangeForThread,
                                  ThreadIdType threadId,
                                  Self *inHolder );

  /** Threader for grandient modification */
  ModifyGradientThreaderType::Pointer      m_ModifyGradientThreader;

};

} // end namespace itk

//#ifndef ITK_MANUAL_INSTANTIATION
//# include "itkGradientDescentObjectOptimizerBase.hxx"
//#endif

#endif
