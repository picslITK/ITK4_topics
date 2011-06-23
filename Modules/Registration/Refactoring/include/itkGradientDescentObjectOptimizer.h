/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGradientDescentObjectOptimizer.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGradientDescentObjectOptimizer_h
#define __itkGradientDescentObjectOptimizer_h

#include "itkPoint.h"
#include "itkIndex.h"
#include "itkGradientDescentObjectOptimizerBase.h"

namespace itk
{

// functor for threading using the metric function class
// assuming function has output allocated already
template<class TMetricFunction>
class ITK_EXPORT GradientDescentObjectOptimizer
  : public GradientDescentObjectOptimizerBase< TMetricFunction >
{
public:
  /** Standard class typedefs. */
  typedef GradientDescentObjectOptimizer     Self;
  typedef GradientDescentObjectOptimizerBase< TMetricFunction >
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
//  itkSuperclassTraitMacro( MetricThreaderType );
//  itkSuperclassTraitMacro( MetricThreaderTypePointer );
need these anymore?
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

  /** Set the learning rate. */
  itkSetMacro(LearningRate, double);

  /** Get the learning rate. */
  itkGetConstReferenceMacro(LearningRate, double);

public:

  /** Start and run the optimization */
  virtual void StartOptimization();

  /** Resume the optimization. Can be called after StopOptimization to
   * resume. The bulk of the optimization work loop is here. */
  virtual void ResumeOptimization();

protected:

  /** Advance one Step following the gradient direction.
   * Includes transform update. */
  virtual void AdvanceOneStep(void);

  /** Modify the gradient over a given index range. */
  virtual void ModifyGradientOverSubRange( IndexRangeType& subrange );

  double m_LearningRate;

  /** Default constructor */
  GradientDescentObjectOptimizer();

  virtual ~GradientDescentObjectOptimizer(){}

private:

  //purposely not implemented
  GradientDescentObjectOptimizer( const Self & );
  void operator=( const Self& );      //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkGradientDescentObjectOptimizer.txx"
#endif

#endif
