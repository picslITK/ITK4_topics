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
/** \class GradientDescentObjectOptimizer
 *
 * Unlike the previous version of GradientDescentOptimizer, this version does
 * not have a "maximize/minimize" option to modify the effect of the metric
 * derivative.
 * The assigned metric is assumed to return a parameter derivative result that
 * "improves" the optimization when added to the current parameters via the
 * metric::UpateTransformParameters method, after the optimizer applies scales
 * and a learning rate.
 */
class ITK_EXPORT GradientDescentObjectOptimizer
  : public GradientDescentObjectOptimizerBase
{
public:
  /** Standard class typedefs. */
  typedef GradientDescentObjectOptimizer       Self;
  typedef GradientDescentObjectOptimizerBase   Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GradientDescentObjectOptimizer,
               GradientDescentObjectOptimizerBase);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Derivative type */
  typedef Superclass::DerivativeType      DerivativeType;

  /** Metric type over which this class is templated */
  typedef Superclass::MeasureType                  MeasureType;
  typedef Superclass::InternalComputationValueType InternalComputationValueType;

  /** Set the learning rate. */
  itkSetMacro(LearningRate, InternalComputationValueType);

  /** Get the learning rate. */
  itkGetConstReferenceMacro(LearningRate, InternalComputationValueType);

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
  virtual void ModifyGradientOverSubRange( const IndexRangeType& subrange );

  InternalComputationValueType  m_LearningRate;

  /** Default constructor */
  GradientDescentObjectOptimizer();

  virtual ~GradientDescentObjectOptimizer(){}

private:

  //purposely not implemented
  GradientDescentObjectOptimizer( const Self & );
  void operator=( const Self& );      //purposely not implemented

};

} // end namespace itk

#endif
