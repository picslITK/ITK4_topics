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

  /** Global derivative accessor */
  const DerivativeType & GetGlobalDerivative(void)
  { return this->m_GlobalDerivative; }

public:

  void StartOptimization() = 0;

protected:

  /** Initialize threading memory before starting optimization.
   * This is a one-time initialization, i.e. not performed at every
   * iteration.
   * Must be called from a derived class during its initialization. */
  void InitializeThreadingMemory();

  /** Perform any per-iteration operations before the metric is
   * updated. This method is not threaded by the optimizer, but
   * can contain calls to threaded methods or filters. */
  virtual void BeforeMetricThreadedGenerateData();

  /** Perform any per-iteration operations after the metric has been
   * updated via threading. In particular this will gather the results
   * from the threads. */
  void AfterMetricThreadedGenerateData();

  /** Advance one step following the gradient direction
   * This method verifies if a change in direction is required
   * and if a reduction in steplength is required. */
//  virtual void AdvanceOneStep(void);

  /** Advance one step along the corrected gradient taking into
   * account the steplength represented by factor.
   * This method is invoked by AdvanceOneStep.
   * \sa AdvanceOneStep */
//  virtual void StepAlongGradient(void) = 0;

  /** Cleanup after optimization is complete. */
  virtual void CleanupAfterOptimization(void);

  /** Callback to perform metric update. Gets assigned to the threader's
   * ThreadedGenerateData. Must be static so it can be used as a callback.
   * An instance of this optimizer class is referenced through
   * \c holder, which is passed in via the threader's user data. */
  static void ComputeMetricValueInRegionThreaded(
                                  const ImageRegionType & regionForThread,
                                  int threadId,
                                  void *inHolder );

  /** Default constructor */
  GradientDescentThreadedMetricOptimizerBase();

  virtual ~GradientDescentThreadedMetricOptimizerBase(){}

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
