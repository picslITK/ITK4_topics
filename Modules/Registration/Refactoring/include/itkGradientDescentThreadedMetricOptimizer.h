/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGradientDescentThreadedMetricOptimizer.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGradientDescentThreadedMetricOptimizer_h
#define __itkGradientDescentThreadedMetricOptimizer_h

#include "itkPoint.h"
#include "itkIndex.h"
#include "itkGradientDescentThreadedMetricOptimizerBase.h"

namespace itk
{

// functor for threading using the metric function class
// assuming function has output allocated already
template<class TMetricFunction, class TMetricThreader>
class ITK_EXPORT GradientDescentThreadedMetricOptimizer
  : public GradientDescentThreadedMetricOptimizerBase< TMetricFunction,
                                                      TMetricThreader >
{
public:
  /** Standard class typedefs. */
  typedef GradientDescentThreadedMetricOptimizer     Self;
  typedef GradientDescentThreadedMetricOptimizerBase < TMetricFunction,
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

  /** Set the learning rate. */
  itkSetMacro(LearningRate, double);

  /** Get the learning rate. */
  itkGetConstReferenceMacro(LearningRate, double);

public:

  /**
   * Start and run the optimization
   */
  virtual void StartOptimization()
  {
    itkDebugMacro("StartOptimization");

    this->InitializeForThreading();
    m_CurrentIteration   = 0;

    /* Validate some settings */
    if( this->m_Scales.Size() != this->m_Metric->GetNumberOfParameters() )
      {
      m_StopCondition = OtherError;
      m_StopConditionDescription << "Scales size error";
      StopOptimization();
      itkExceptionMacro("Size of scales (" << this->m_Scales.Size()
                        << ") must match size of parameters (" <<
                        this->m_Metric->GetNumberOfParameters() << ").");
      }

    this->ResumeOptimization();
  }

  /**
   * Resume optimization.
   */
  virtual void ResumeOptimization()
  {
    m_StopConditionDescription.str("");
    m_StopConditionDescription << this->GetNameOfClass() << ": ";
    InvokeEvent( StartEvent() );

    m_Stop = false;
    while( ! m_Stop )
      {
      /* Compute value/derivative, using threader. */
      try
        {
        this->UpdateMetricValueAndDerivative();
        }
      catch
        {
        m_StopCondition = MetricError;
        m_StopConditionDescription << "Metric error";
        StopOptimization();

        // Pass exception to caller
        throw err;
        }

      /* Check if optimization has been stopped externally.
       * (Presumably this could happen from a multi-threaded client app?) */
      if ( m_Stop )
        {
        m_StopConditionDescription << "StopOptimization() called";
        break;
        }

      /* Advance one step */
      this->AdvanceOneStep();

      /* update transforms */
      /* this is similar to StepAlongGradient */
      /* TODO: needs to be threaded */
      //virtual opt func to do updates, by calling Transform::Update, maybe
      //with scaling first, or other?


      /* Update and check iteration count */
      m_CurrentIteration++;
      if ( m_CurrentIteration >= m_NumberOfIterations )
        {
        m_StopConditionDescription << "Maximum number of iterations ("
                                   << m_NumberOfIterations
                                   << ") exceeded.";
        m_StopCondition = MaximumNumberOfIterations;
        StopOptimization();
        break;
        }
      }

    this->Cleanup();
  }


protected:

/**
 * Advance one Step following the gradient direction
 */
virtual void AdvanceOneStep(void)
{
  itkDebugMacro("AdvanceOneStep");

  this->ModifyGradient();
  this->UpdateParameters();
  this->InvokeEvent( IterationEvent() );
}

  /*
   * Modify the gradient over a given range
   */
  virtual void ModifyGradientOverSubRange( IndexRangeType& subrange )
  {
    double direction;
    if ( this->m_Maximize )
      {
      direction = 1.0;
      }
    else
      {
      direction = -1.0;
      }

    ScalesType& scales = this->GetScales();

    double product = direction * m_LearningRate;
    for ( unsigned int j = subrange[0]; j <= subrange[1]; j++ )
      {
      m_Gradient[j] = m_Gradient[j] / scales[j] * product;
      }
  }

  /** Advance one step along the corrected gradient taking into
   * account the steplength represented by factor.
   * This method is invoked by AdvanceOneStep.
   * \sa AdvanceOneStep */
  virtual void StepAlongGradient(void)

  double m_LearningRate;

  /** Default constructor */
  GradientDescentThreadedMetricOptimizer()
  {
  m_LearningRate = 1.0;
  }

  virtual ~GradientDescentThreadedMetricOptimizer(){}

private:

  //purposely not implemented
  GradientDescentThreadedMetricOptimizer( const Self & );
  void operator=( const Self& );      //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkGradientDescentThreadedMetricOptimizer.txx"
#endif

#endif
