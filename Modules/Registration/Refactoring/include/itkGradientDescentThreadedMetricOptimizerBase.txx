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
#ifndef __itkGradientDescentThreadedMetricOptimizerBase_txx
#define __itkGradientDescentThreadedMetricOptimizerBase_txx

#include "itkGradientDescentThreadedMetricOptimizerBase.h"

namespace itk
{

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::GradientDescentThreadedMetricOptimizerBase()
{
  /* Point the metric threader to the threading worker callback.
   * The rest of the threader is initialed in Superclass. */
  this->m_MetricThreader->SetThreadedGenerateData(
    Self::ComputeMetricValueInRegionThreaded );

  this->m_ModifyGradientThreader->SetThreadedGenerateData(
    Self::ModifyGradientThreaded );

  m_Maximize = false;
  m_NumberOfIterations = 100;
  m_CurrentIteration = 0;
  m_StopCondition = MaximumNumberOfIterations;
  m_StopConditionDescription << this->GetNameOfClass() << ": ";
}

//-------------------------------------------------------------------
template<class TMetricFunction>
const std::string
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::GetStopConditionDescription() const
{
  return m_StopConditionDescription.str();
}

//-------------------------------------------------------------------
template<class TMetricFunction>
const std::string
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::StopOptimization(void)
{
  itkDebugMacro("StopOptimization");

  m_Stop = true;
  this->CleanupFromThreading();
  InvokeEvent( EndEvent() );
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::InitializeForThreading()
{
  if( this->m_Metric.IsNull() )
    {
    itkExceptionMacro("m_Metric must be set.");
    return;
    }

  // TODO: call m_Metric->Initialize(), or check that
  // it's been init'ed somehow? Would all metrics have Initialize(), or
  // at least all new ones to be used with this optimzier hierarchy?

/*

for regular step grad descent:

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

*/
  //Allocate and initialize memory for holding derivative results.
  this->m_DerivativesPerThread.resize( this->m_NumberOfThreads );
  this->m_MeasurePerThread.resize( this->m_NumberOfThreads );
  unsigned long globalDerivativeSize =
    this->m_Metric->GetMovingImageTransform()->GetNumberOfParameters();
  std::cout << "  Initialize: deriv size  "
            << globalDerivativeSize << std::endl;
  this->m_GlobalDerivative.SetSize( globalDerivativeSize );
  /* For transforms with local support, e.g. deformation field,
   * use a single derivative container that's updated by region
   * in multiple threads. */
  if ( this->m_Metric->GetMovingImageTransform()->HasLocalSupport() )
    {
    std::cout << " Initialize: tx HAS local support\n";
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

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::UpdateMetricValueAndDerivative()
{
  this->BeforeMetricThreadedGenerateData();
  /* Calculate new metric value and derivative with threading. */
  this->m_MetricThreader->GenerateData();
  /* Collect metric results from the threads. */
  this->AfterMetricThreadedGenerateData();
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::BeforeMetricThreadedGenerateData()
{
  /* Allow the metric to do any per-iteration initialization it
   * requires. This method is not expected to be thread-safe, and
   * may implement its own threading */
  this->m_Metric->BeforeThreadedGetValueAndDerivative();
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::AfterMetricThreadedGenerateData()
{
  /* For global transforms, sum the derivatives from each region. */
  if (  ! this->m_Metric->GetMovingImageTransform()->HasLocalSupport() )
    {
    this->m_GlobalDerivative.Fill(0);
    for (int i=0; i<this->m_NumberOfThreads; i++)
      {
      this->m_GlobalDerivative += this->m_DerivativesPerThread[i];
      }
    }

  /* Accumulate the metric value from threads and store */
  this->m_Value =
    NumericTraits<InternalComputationValueType>::Zero;
  for(unsigned int i=0; i< this->m_MeasurePerThread.size(); i++)
    {
    this->m_Value += this->m_MeasurePerThread[i];
    }

  std::cout << " end after threaded generate data. global derivative: "
            << this->m_GlobalDerivative << std::endl;
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::CleanupFromThreading()
{
  // Free some memory used during threading. This probably
  // doesn't actually amount to much.
  // We could also free m_GlobalDerivatives by setting its size to 0,
  // but user may want to examine it when optimization is stopped.
  if ( ! this->m_Metric->GetMovingImageTransform()->HasLocalSupport() )
    {
    for (int i=0; i<this->m_NumberOfThreads; i++)
      {
      this->m_DerivativesPerThread[i].SetSize(0);
      }
    }
}

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::ComputeMetricValueInRegionThreaded( const ImageRegionType & regionForThread,
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

//-------------------------------------------------------------------
template<class TMetricFunction>
void
GradientDescentThreadedMetricOptimizerBase<TMetricFunction>
::ModifyGradientThreaded( const IndexRangeType& rangeForThread,
                          int threadId,
                          void *inHolder )
{
  Self * holder = static_cast<Self*>(inHolder);
  holder->ModifyGradientOverSubRange( rangeForThread );
}

}//namespace itk

#endif
