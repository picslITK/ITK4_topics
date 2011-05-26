/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkObjectToObjectThreadedMetricOptimizerBase.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkObjectToObjectThreadedMetricOptimizerBase_h
#define __itkObjectToObjectThreadedMetricOptimizerBase_h

#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"


namespace itk
{

// functor for threading using the metric function class
// assuming function has output allocated already
template<class TMetricFunction, class TMetricThreader>
class ITK_EXPORT ObjectToObjectThreadedMetricOptimizerBase : public Object
{
public:
  /** Standard class typedefs. */
  typedef ObjectToObjectThreadedMetricOptimizerBase   Self;
  typedef Object                                      Superclass;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, Superclass);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /**  Scale type. */
  typedef TransformParameters< double >             ScalesType;

  /** Metric type over which this class is templated */
  typedef TMetricFunction                           MetricType;
  typedef typename MetricType::Pointer              MetricTypePointer;
  /** Measure type */
  typedef typename MetricType::MeasureType          MeasureType;
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
  /** Metric Threader type */
  typedef TMetricThreader                           MetricThreaderType;
  typedef typename MetricThreaderType::Pointer      MetricThreaderTypePointer;

  /** Accessors for Metric */
  itkGetObjectMacro( Metric, MetricType );
  itkSetObjectMacro( Metric, MetricType );

  /** Accerssor for metric value */
  itkGetConstReferenceMacro( Value, InternalComputationValueType );

  /** Accessor for Metric Threader */
  itkGetObjectMacro( MetricThreader, MetricThreaderType );

  /** Set current parameters scaling. */
  void SetScales(const ScalesType & scales)
  {
    itkDebugMacro("setting scales to " <<  scales);
    m_Scales = scales;
    this->Modified();
  }

  /** Get current parameters scaling. */
  itkGetConstReferenceMacro(Scales, ScalesType);

  /** Set the number of threads to use when threading.
   * This is initialized by default to the global default number of
   * threads from itkMultiThreader */
  virtual void SetNumberOfThreads( int number )
  {
    if( number < 1 )
      {
      itkExceptionMacro("Number of threads must be > 0");
      }
    if( number != this->m_NumberOfThreads )
      {
      this->m_NumberOfThreads = number;
      this->m_MetricThreader->SetNumberOfThreads( number );
      this->Modified();
      }
  }

  /** Set the full image region over which to optimize */
  void SetOverallRegion( ImageRegionType & region )
  {
    m_OverallRegion = region;
    m_MetricThreader->SetOverallRegion( region );
    this->Modified();
  }

  /** Run the optimization */
  virtual void StartOptimization() = 0;

protected:

  /** Default constructor */
  ObjectToObjectThreadedMetricOptimizerBase()
  {
    /* Setup Metric threader. The callback must be set
     * from derived class. */
    this->m_MetricThreader = MetricThreaderType::New();
    this->m_MetricThreader->SetHolder( static_cast<void*>(this) );
    this->SetNumberOfThreads( this->m_MetricThreader->
      GetMultiThreader()->GetGlobalDefaultNumberOfThreads() );
    this->m_Metric = NULL;
    this->m_Value = 0;
  }

  MetricTypePointer             m_Metric;
  int                           m_NumberOfThreads;
  MetricThreaderTypePointer     m_MetricThreader;

  /** Metric measure value at a given iteration */
  MeasureType                                 m_Value;
  /** Intermediary metric value storage. */
  std::vector<InternalComputationValueType>   m_MeasurePerThread;

  virtual ~ObjectToObjectThreadedMetricOptimizerBase(){}

private:

  ImageRegionType         m_OverallRegion;
  ScalesType              m_Scales;

  //purposely not implemented
  ObjectToObjectThreadedMetricOptimizerBase( const Self & );
  void operator=( const Self& );      //purposely not implemented

};

} // end namespace itk

#endif
