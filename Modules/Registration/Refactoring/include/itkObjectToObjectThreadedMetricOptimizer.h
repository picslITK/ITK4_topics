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

#include "itkObjectToObjectMetric.h"
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
struct ObjectToObjectThreadedMetricOptimizer{

  typedef ObjectToObjectThreadedMetricOptimizer          Self;

  typedef TMetricFunction           MetricType;
  typedef typename MetricType::Pointer  MetricTypePointer;
  typedef typename MetricType::MeasureType  MeasureType;
  typedef typename MetricType::DerivativeType  DerivativeType;
  typedef typename MetricType::InternalComputationValueType InternalComputationValueType;
  typedef typename MetricType::RegionType ImageRegionType;
  typedef typename MetricType::FixedImageType ImageType;
  typedef typename MetricType::FixedImagePointer FixedImagePointer;
  typedef typename MetricType::MovingImagePointer MovingImagePointer;
  typedef typename MetricType::TransformPointer TransformPointer;

public:
  MetricTypePointer           metric;
  std::vector<InternalComputationValueType> measure_per_thread;
  DerivativeType global_derivative;
  std::vector<DerivativeType> derivatives_per_thread;

  InternalComputationValueType AccumulateMeasuresFromAllThreads() {
    InternalComputationValueType energy = NumericTraits<InternalComputationValueType>::Zero;
    for(unsigned int i=0; i<measure_per_thread.size(); i++) energy += measure_per_thread[i];
    return energy;
  }

  void BeforeThreadedGenerateData(unsigned int number_of_threads )
  {
    this->derivatives_per_thread.resize(number_of_threads);
    this->measure_per_thread.resize(number_of_threads);
    unsigned long global_derivative_size=this->metric->GetMovingImageTransform()->GetNumberOfParameters();
    std::cout << "  before threaded generate data deriv size  " << global_derivative_size << std::endl;
    this->global_derivative.SetSize(global_derivative_size);
    if ( this->metric->GetMovingImageTransform()->HasLocalSupport() )
      {
        for (unsigned int i=0; i<number_of_threads; i++) {
          this->derivatives_per_thread[i]=this->global_derivative;
        }
      }
    else
      {
        for (unsigned int i=0; i<number_of_threads; i++) {
          DerivativeType thread_derivative(global_derivative_size);
          thread_derivative.Fill(0);
          this->derivatives_per_thread[i]=thread_derivative;
        }
      }
    std::cout << " end before threaded generate data " << std::endl;
  }

  void AfterThreadedGenerateData(unsigned int number_of_threads )
  {
    if (  ! this->metric->GetMovingImageTransform()->HasLocalSupport() )
      {
        for (unsigned int i=0; i<number_of_threads; i++) {
          this->global_derivative=this->global_derivative+derivatives_per_thread[i];
        }
      }
    std::cout << " end after threaded generate data " <<  this->global_derivative << std::endl;
  }

  static void ComputeMetricValueInRegionThreaded(const ImageRegionType &regionForThread, int threadId,  Self *holder){
    //    std::cout << regionForThread << std::endl;
    InternalComputationValueType local_metric;
    holder->measure_per_thread[threadId] = NumericTraits<InternalComputationValueType>::Zero;
    /** Compute one iteration of the metric */
    local_metric=holder->metric->ComputeMetricAndDerivative(regionForThread, holder->derivatives_per_thread[threadId] );
    holder->measure_per_thread[threadId] += local_metric;
  }


};




} // end namespace itk


#endif
