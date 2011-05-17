/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkObjectToObjectThreadedMetricOptimizer.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkObjectToObjectThreadedMetricOptimizer_h
#define __itkObjectToObjectThreadedMetricOptimizer_h

#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkIndex.h"
#include "itkTransform.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMultiThreader.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"


namespace itk
{

// functor for threading using the metric function class
// assuming function has output allocated already
template<class TMetricFunction>
class ITK_EXPORT ObjectToObjectThreadedMetricOptimizer : public Object
{
public:
  /** Standard class typedefs. */
  typedef ObjectToObjectThreadedMetricOptimizer     Self;
  typedef Object                                    Superclass;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, Superclass);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Metric type over which this class is templated */
  typedef TMetricFunction                           MetricType;
  typedef typename MetricType::Pointer              MetricTypePointer;
  /** Measure type */
  typedef typename MetricType::MeasureType          MeasureType;
  /** Derivative type */
  typedef typename MetricType::DerivativeType       DerivativeType;
  /** Image region type */
  typedef typename MetricType::RegionType           ImageRegionType;
  /** Fixed image type */
  typedef typename MetricType::FixedImageType       FixedImageType;
  /** Fixed image pointer */
  typedef typename MetricType::FixedImagePointer    FixedImagePointer;
  /** Moving image pointer */
  typedef typename MetricType::MovingImagePointer   MovingImagePointer;
  /** Tranform pointer */
  typedef typename MetricType::TransformPointer     TransformPointer;
  /** Internal computation type, for maintaining a desired precision */
  typedef typename MetricType::InternalComputationValueType InternalComputationValueType;

  /** Accessors for Metric */
  itkGetObjectMacro( Metric, MetricType );
  itkSetObjectMacro( Metric, MetricType );

  const DerivativeType & GetGlobalDerivative(void)
  {
    return m_GlobalDerivative;
  }

public:
  InternalComputationValueType AccumulateMeasuresFromAllThreads()
  {
    InternalComputationValueType energy =
      NumericTraits<InternalComputationValueType>::Zero;
    for(unsigned int i=0; i<m_MeasurePerThread.size(); i++)
      {
      energy += m_MeasurePerThread[i];
      }
    return energy;
  }

  void BeforeThreadedGenerateData(unsigned int number_of_threads )
  {
    this->m_DerivativesPerThread.resize(number_of_threads);
    this->m_MeasurePerThread.resize(number_of_threads);
    unsigned long globalDerivativeSize =
      this->m_Metric->GetMovingImageTransform()->GetNumberOfParameters();
    std::cout << "  before threaded generate data deriv size  "
              << globalDerivativeSize << std::endl;
    this->m_GlobalDerivative.SetSize( globalDerivativeSize );
    /* For global transforms with local support, e.g. deformation field,
     * use a single global derivative container that's updated by region
     * in multiple threads. */
    if ( this->m_Metric->GetMovingImageTransform()->HasLocalSupport() )
      {
      for (unsigned int i=0; i<number_of_threads; i++)
        {
        this->m_DerivativesPerThread[i].SetData(
                                      this->m_GlobalDerivative.data_block(),
                                      this->m_GlobalDerivative.Size(),
                                      false );
        }
      }
    else
      {
      for (unsigned int i=0; i<number_of_threads; i++)
        {
        this->m_DerivativesPerThread[i].SetSize( globalDerivativeSize );
        this->m_DerivativesPerThread[i].Fill( 0 );
        }
      }
    std::cout << " end before threaded generate data " << std::endl;
  }

  void AfterThreadedGenerateData(unsigned int number_of_threads )
  {
    if (  ! this->m_Metric->GetMovingImageTransform()->HasLocalSupport() )
      {
      for (unsigned int i=0; i<number_of_threads; i++)
        {
        this->m_GlobalDerivative += this->m_DerivativesPerThread[i];
        this->m_DerivativesPerThread[i].SetSize(0); //free memory
        }
      }
    std::cout << " end after threaded generate data "
              << this->m_GlobalDerivative << std::endl;
  }

  static void ComputeMetricValueInRegionThreaded(
                                  const ImageRegionType & regionForThread,
                                  int threadId,
                                  Self *holder )
  {
    //    std::cout << regionForThread << std::endl;
    InternalComputationValueType local_metric;
    holder->m_MeasurePerThread[threadId] =
      NumericTraits<InternalComputationValueType>::Zero;
    /** Compute one iteration of the metric */
    local_metric=holder->m_Metric->ComputeMetricAndDerivative(
                                    regionForThread,
                                    holder->m_DerivativesPerThread[threadId] );
    holder->m_MeasurePerThread[threadId] += local_metric;
  }

private:
  MetricTypePointer       m_Metric;
  DerivativeType          m_GlobalDerivative;

  std::vector< DerivativeType >               m_DerivativesPerThread;
  std::vector<InternalComputationValueType>   m_MeasurePerThread;

};

} // end namespace itk

#endif
