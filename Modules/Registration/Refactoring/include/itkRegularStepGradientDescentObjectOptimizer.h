/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRegularStepGradientDescentObjectOptimizer.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRegularStepGradientDescentObjectOptimizer_h
#define __itkRegularStepGradientDescentObjectOptimizer_h

#include "itkPoint.h"
#include "itkIndex.h"
#include "itkGradientDescentObjectOptimizerBase.h"

/*
 * NOTE: Not working or uptodate.
 *
 * See GradientDescentObjectOptimizer.
 */

namespace itk
{

template<class TMetricFunction>
class ITK_EXPORT RegularStepGradientDescentObjectOptimizer
  : public GradientDescentObjectOptimizerBase<TMetricFunction>
{
public:
  /** Standard class typedefs. */
  typedef RegularStepGradientDescentObjectOptimizer     Self;
  typedef GradientDescentObjectOptimizerBase <TMetricFunction>
                                                                Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(RegularStepGradientDescentObjectOptimizer,
               GradientDescentObjectOptimizerBase);

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
  {
    return this->m_GlobalDerivative;
  }

public:

  void StartOptimization()
  {

    this->Initialize();

    bool done = false;

    do
      {
      /* Check # of iterations */

    /* NOTE - a method for threaded metric computation can go into the GradDescent base */
      /* Compute value/derivative, one step.
       * Use threader. */
      /* per-iteration preparation for threading */
      this->BeforeMetricThreadedGenerateData();
      this->m_MetricThreader->GenerateData();
      /* Collect metric results from the threads. */
      this->AfterMetricThreadedGenerateData();

      /* Advance one step */

          /* update transforms */
          /* this is similar to StepAlongGradient */
          /* TODO: needs to be threaded */
          //virtual opt func to do updates, by calling Transform::Update, maybe
          //with scaling first, or other?

      done = true;
      }while( ! done );

    this->Cleanup();
  }

protected:

  /** Accumulate the metric values computed by each thread */
  virtual InternalComputationValueType AccumulateMeasuresFromAllThreads()
  {
    InternalComputationValueType energy =
      NumericTraits<InternalComputationValueType>::Zero;
    for(unsigned int i=0; i< this->m_MeasurePerThread.size(); i++)
      {
      energy += this->m_MeasurePerThread[i];
      }
    return energy;
  }

  /** Initialize and verify settings before starting optimization */
  void Initialize()
  {

    if( this->m_Metric.IsNull() )
      {
      itkExceptionMacro("m_Metric must be set.");
      return;
      }

    // TODO: call m_Metric->Initialize(), or check that
    // it's been init'ed somehow? Do all metrics have Initialize(), or
    // at least all new ones to be used with this optimzier hierarchy?

/*
    const unsigned int dimensions =
      this->m_Metric->GetNumberOfParameters();
    // Verify parameter settings
    if ( m_RelaxationFactor < 0.0 )
      {
      itkExceptionMacro(<< "Relaxation factor must be positive. "
                           "Current value is " << m_RelaxationFactor);
      return;
      }

    if ( m_RelaxationFactor >= 1.0 )
      {
      itkExceptionMacro(<< "Relaxation factor must less than 1.0. "
                           "Current value is " << m_RelaxationFactor);
      return;
      }

    // Make sure the scales have been set properly
    if ( m_Scales.size() != dimensions )
      {
      itkExceptionMacro(<< "The size of Scales is "
                        << m_Scales.size()
                        << ", but the NumberOfParameters for the Metric is "
                        << dimensions
                        << ".");
      }
*/
    //Allocate and initialize memory for holding derivative results.
    this->m_DerivativesPerThread.resize( this->m_NumberOfThreads );
    this->m_MeasurePerThread.resize( this->m_NumberOfThreads );
    unsigned long globalDerivativeSize =
      this->m_Metric->GetMovingTransform()->GetNumberOfParameters();
    std::cout << "  Initialize: deriv size  "
              << globalDerivativeSize << std::endl;
    this->m_GlobalDerivative.SetSize( globalDerivativeSize );
    /* For transforms with local support, e.g. deformation field,
     * use a single derivative container that's updated by region
     * in multiple threads. */
    if ( this->m_Metric->GetMovingTransform()->HasLocalSupport() )
      {
      std::cout << " Initialize: tx has local support\n";
      for (int i=0; i<this->m_NumberOfThreads; i++)
        {
        this->m_DerivativesPerThread[i].SetData(
                                      this->m_GlobalDerivative.data_block(),
                                      this->m_GlobalDerivative.Size(),
                                      false );
        }
      }
    else
      {
      std::cout << " Initialize: tx does NOT have local support\n";
      /* Global transforms get a separate derivatives container for each thread
       * that holds the result for a particular region. */
      for (int i=0; i<this->m_NumberOfThreads; i++)
        {
        this->m_DerivativesPerThread[i].SetSize( globalDerivativeSize );
        /* Be sure to init to 0 here, because the threader may not use
         * all the threads if the region is better split into fewer
         * subregions. */
        this->m_DerivativesPerThread[i].Fill( 0 );
        }
      }
    std::cout << " end Initialize " << std::endl;
  }

  void BeforeMetricThreadedGenerateData()
  {
  }

  void AfterMetricThreadedGenerateData()
  {
    /* For global transforms, sum the derivatives from each region. */
    if (  ! this->m_Metric->GetMovingTransform()->HasLocalSupport() )
      {
      this->m_GlobalDerivative.Fill(0);
      for (int i=0; i<this->m_NumberOfThreads; i++)
        {
        this->m_GlobalDerivative += this->m_DerivativesPerThread[i];
        }
      }

    /* Get and store the metric value */
    this->m_Value = this->AccumulateMeasuresFromAllThreads();

    std::cout << " end after threaded generate data. global derivative: "
              << this->m_GlobalDerivative << std::endl;
  }

  /** Advance one step following the gradient direction
   * This method verifies if a change in direction is required
   * and if a reduction in steplength is required. */
  virtual void AdvanceOneStep(void);

  /** Advance one step along the corrected gradient taking into
   * account the steplength represented by factor.
   * This method is invoked by AdvanceOneStep.
   * \sa AdvanceOneStep */
  virtual void StepAlongGradient(void){}

  /** Cleanup after optimization is complete. */
  void Cleanup(void)
  {
    //Free some memory
    this->m_GlobalDerivative.SetSize(0);
    if ( ! this->m_Metric->GetMovingTransform()->HasLocalSupport() )
      {
      for (int i=0; i<this->m_NumberOfThreads; i++)
        {
        this->m_DerivativesPerThread[i].SetSize(0);
        }
      }
  }

  /** Callback that gets assigned to the threader's ThreadedGenerateData.
   * Make it static so it can be used as a callback.
   * An instance of this optimizer class is referenced through
   * \c holder, which is passed in via the threader's user data. */
  static void ComputeMetricValueInRegionThreaded(
                                  const ImageRegionType & regionForThread,
                                  int threadId,
                                  void *inHolder )
  {
    //    std::cout << regionForThread << std::endl;
    InternalComputationValueType local_metric;
    Self * holder = static_cast<Self*>(inHolder);
    /** Compute one iteration of the metric */
    local_metric = holder->m_Metric->ComputeMetricAndDerivative(
                                    regionForThread,
                                    holder->m_DerivativesPerThread[threadId] );
    holder->m_MeasurePerThread[threadId] = local_metric;
  }

  /** Default constructor */
  RegularStepGradientDescentObjectOptimizer()
  {
    /* Point the threader to the threading worker callback.
     * The rest of the threader is initialed in Superclass. */
    this->m_MetricThreader->SetThreadedGenerateData(
      Self::ComputeMetricValueInRegionThreaded );
  }

  virtual ~RegularStepGradientDescentObjectOptimizer(){}

private:

  DerivativeType                  m_GlobalDerivative;
  std::vector< DerivativeType >   m_DerivativesPerThread;

  //purposely not implemented
  RegularStepGradientDescentObjectOptimizer( const Self & );
  void operator=( const Self& );      //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkRegularStepGradientDescentObjectOptimizer.txx"
#endif

#endif
