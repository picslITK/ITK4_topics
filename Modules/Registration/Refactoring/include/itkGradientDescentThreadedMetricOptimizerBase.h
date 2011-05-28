/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGradientDescentThreadedMetricOptimizerBase.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGradientDescentThreadedMetricOptimizerBase_h
#define __itkGradientDescentThreadedMetricOptimizerBase_h

#include "itkPoint.h"
#include "itkIndex.h"
#include "itkObjectToObjectThreadedMetricOptimizerBase.h"
#include "itkArray1DToData.h"

namespace itk
{

// functor for threading using the metric function class
// assuming function has output allocated already
template<class TMetricFunction, class TMetricThreader>
class ITK_EXPORT GradientDescentThreadedMetricOptimizerBase
  : public ObjectToObjectThreadedMetricOptimizerBase< TMetricFunction,
                                                      TMetricThreader >
{
public:
  /** Standard class typedefs. */
  typedef GradientDescentThreadedMetricOptimizerBase     Self;
  typedef ObjectToObjectThreadedMetricOptimizerBase < TMetricFunction,
                                                      TMetricThreader >
                                                        Superclass;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, Superclass);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Codes of stopping conditions. */
  typedef enum {
    MaximumNumberOfIterations,
    MetricError,
    OtherError
    } StopConditionType;

  /** Derivative type */
  typedef typename TMetricFunction::DerivativeType  DerivativeType;

  /** Metric type over which this class is templated */
  itkSuperclassTraitMacro( MetricType );
  itkSuperclassTraitMacro( MetricTypePointer );
  /** Threader type */
  itkSuperclassTraitMacro( MetricThreaderType );
  itkSuperclassTraitMacro( MetricThreaderTypePointer );
  /** Measure type */
  itkSuperclassTraitMacro( MeasureType );
  /** Image region type */
  itkSuperclassTraitMacro( ImageRegionType );
  /** Fixed image type */
  itkSuperclassTraitMacro( FixedImageType );
  /** Fixed image pointer */
  itkSuperclassTraitMacro( FixedImagePointer );
  /** Moving image pointer */
  itkSuperclassTraitMacro( MovingImagePointer );
  /** Tranform pointer */
  itkSuperclassTraitMacro( TransformPointer );
  /** Internal computation type, for maintaining a desired precision */
  itkSuperclassTraitMacro( InternalComputationValueType );

  /** Threader for grandient update */
  typedef Array1DToData                 ModifyGradientThreaderType;
  typedef ModifyGradientThreaderType::IndexRangeType  IndexRangeType;

  /** Global derivative accessor */
  const DerivativeType & GetGlobalDerivative(void)
  { return this->m_GlobalDerivative; }

  /** Get stop condition enum */
  itkGetConstReferenceMacro(StopCondition, StopConditionType);

  /** Methods to configure the cost function. */
  itkGetConstReferenceMacro(Maximize, bool);
  itkSetMacro(Maximize, bool);
  itkBooleanMacro(Maximize);
  bool GetMinimize() const
  { return !m_Maximize; }
  void SetMinimize(bool v)
  { this->SetMaximize(!v); }
  void MinimizeOn()
  { this->MaximizeOff(); }
  void MinimizeOff()
  { this->MaximizeOn(); }

  /** Set the number of iterations. */
  itkSetMacro(NumberOfIterations, SizeValueType);

  /** Get the number of iterations. */
  itkGetConstReferenceMacro(NumberOfIterations, SizeValueType);

  /** Get the current iteration number. */
  itkGetConstMacro(CurrentIteration, SizeValueType);

  /** Start and run the optimization */
  virtual void StartOptimization() = 0;

  /** Resume optimization.
   * This runs the optimization loop, and allows continuation
   * of stopped optimization */
  virtual void ResumeOptimization() = 0;

  /** Stop optimization. The object is left in a state so the
   * optimization can be resumed by calling ResumeOptimization. */
  virtual void StopOptimization(void);

protected:

  /** Update the metric's value and derivative for the current iteration.
   * This wraps the threading process by calling
   * BeforeMetricThreadedGenerateData, GenerateData and
   * AfterMetricThreadedGenerateData.
   * Results are stored in member variables. */
  virtual void UpdateMetricValueAndDerivative();

  /** Initialize threading memory before starting optimization.
   * This is a one-time initialization, i.e. not performed at every
   * iteration.
   * Must be called from a derived class during its initialization. */
  void InitializeForThreading();

  /** Perform any per-iteration operations before the metric is
   * updated. This method is not threaded by the optimizer, but
   * can contain calls to threaded methods or filters. */
  virtual void BeforeMetricThreadedGenerateData();

  /** Callback to perform metric update. Gets assigned to the threader's
   * ThreadedGenerateData. Must be static so it can be used as a callback.
   * An instance of this optimizer class is referenced through
   * \c holder, which is passed in via the threader's user data. */
  static void ComputeMetricValueInRegionThreaded(
                                  const ImageRegionType & regionForThread,
                                  int threadId,
                                  void *inHolder );

  /** Perform per-iteration operations after the metric values have been
   * updated via threading. In particular this will gather the results
   * from the threads. */
  virtual void AfterMetricThreadedGenerateData();

  /** Modify the gradient in place, to advance the optimization.
   * This call performs a threaded modification for transforms with
   * local support (assumed to be dense). Otherwise the modification
   * is performed w/out threading.
   * The work is done in ModifyGradientOverSubRange in both cases.
   * At completion, m_Gradient can be used to update the transform
   * parameters. Other members may hold additional results in derived classes.
   */
  virtual void ModifyGradient()
  {
    /* Perform the modification either with or without threading */
    if( this->m_Metric->GetMovingImageTransform()->HasLocalSupport() )
      {
      /* This ends up calling ModifyGradientThreaded frome each thread */
      this->m_ModifyGradientThreader->GenerateData();
      }
    else
      {
      IndexRangeType fullrange;
      fullrange[0] = 0;
      fullrange[1] = this->m_Gradient.Size()-1; //range is inclusive
      this->m_ModifyGradientOverSubRange( fullrange );
      }
  }
  virtual void ModifyGradientOverSubRange( IndexRangeType& subrange ) = 0;

  /** Callback used during threaded gradient modification.
   * Gets assigned to the modify-gradient threader's
   * ThreadedGenerateData. Must be static so it can be used as a callback.
   * It simply calls ModifyGradientOverSubRange, which is where
   * derived classes should put their class-specific modification code.
   * An instance of this optimizer class is referenced through
   * \c holder, which is passed in via the threader's user data. */
  static void ModifyGradientThreaded(
                                  const IndexRangeType& rangeForThread,
                                  int threadId,
                                  void *inHolder );

  /** Advance one step following the gradient direction
   * This method verifies if a change in direction is required
   * and if a reduction in steplength is required. */
//  virtual void AdvanceOneStep(void);

  /** Advance one step along the corrected gradient taking into
   * account the steplength represented by factor.
   * This method is invoked by AdvanceOneStep.
   * \sa AdvanceOneStep */
//  virtual void StepAlongGradient(void) = 0;

  /** Cleanup after optimization is complete.
   * Should be called at end of StartOptimization by derived class. */
  virtual void CleanupFromThreading(void);

  /** Get the reason for termination */
  virtual const std::string GetStopConditionDescription() const;

  /* Common variables for optimization control and reporting */
  bool               m_Stop;
  StopConditionType  m_StopCondition;
  std::ostringstream m_StopConditionDescription;
  SizeValueType      m_NumberOfIterations;
  SizeValueType      m_CurrentIteration;
  bool               m_Maximize;

  /** Default constructor */
  GradientDescentThreadedMetricOptimizerBase();
  virtual ~GradientDescentThreadedMetricOptimizerBase(){}

  /** Threader for grandient modification */
  ModifyGradientThreaderType      m_ModifyGradientThreader;

private:

  DerivativeType                  m_GlobalDerivative;
  std::vector< DerivativeType >   m_DerivativesPerThread;

  //purposely not implemented
  GradientDescentThreadedMetricOptimizerBase( const Self & );
  void operator=( const Self& );      //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkGradientDescentThreadedMetricOptimizerBase.txx"
#endif

#endif
