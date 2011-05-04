/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkObjectToObjectThreadedMetricHolder.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkObjectToObjectThreadedMetricHolder_h
#define __itkObjectToObjectThreadedMetricHolder_h

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
struct ObjectToObjectThreadedMetricHolder{

  typedef ObjectToObjectThreadedMetricHolder          Self;

  typedef TMetricFunction           MetricType;
  typedef typename MetricType::Pointer  MetricTypePointer;
  typedef typename MetricType::MeasureType  MeasureType;
  typedef typename MetricType::InternalComputationValueType InternalComputationValueType;
  typedef typename MetricType::RegionType ImageRegionType;
  typedef typename MetricType::FixedImageType ImageType;
  typedef typename MetricType::FixedImagePointer FixedImagePointer;
  typedef typename MetricType::MovingImagePointer MovingImagePointer;
  typedef typename MetricType::TransformPointer TransformPointer;

public:
  MetricTypePointer           metric;
  std::vector<InternalComputationValueType> measure_per_thread;

  InternalComputationValueType AccumulateMeasuresFromAllThreads() {
    InternalComputationValueType energy = NumericTraits<InternalComputationValueType>::Zero;
    for(unsigned int i=0; i<measure_per_thread.size(); i++) energy += measure_per_thread[i];
    return energy;
  }

  static void ComputeMetricValueInRegionOnTheFlyThreaded(const ImageRegionType &regionForThread, int threadId,  Self *holder){

    //    std::cout << regionForThread << std::endl;
    InternalComputationValueType local_metric;
    holder->measure_per_thread[threadId] = NumericTraits<InternalComputationValueType>::Zero;
    /** Compute one iteration of the metric */
    local_metric=holder->metric->ComputeMetricAndDerivative(regionForThread);
    holder->measure_per_thread[threadId] += local_metric;

  }


};




} // end namespace itk


#endif
