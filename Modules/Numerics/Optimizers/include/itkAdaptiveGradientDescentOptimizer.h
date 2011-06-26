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
#ifndef __itkAdaptiveGradientDescentOptimizer_h
#define __itkAdaptiveGradientDescentOptimizer_h

#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"

#include "itkIntTypes.h"
#include "itkAutomaticGradientDescentOptimizer.h"
#include "itkOptimizerHelper.h"
#include <string>
namespace itk
{
/** \class AdaptiveGradientDescentOptimizer
 * \brief Implement a gradient descent optimizer with adaptive step sizes.
 *
 * This class estimates the learning rate according to the paper
 * "Adaptive Stochastic Gradient Descent Optimisation for Image Registration",
 * Stefan Klein, et al. http://dx.doi.org/10.1007/s11263-008-0168-y
 *
 * A simplification is made in estimating m_Param_a such that the first step
 * yields a specific voxel shift, delta in the paper.
 *
 * \ingroup Numerics Optimizers
 * \ingroup ITK-Optimizers
 */
class ITK_EXPORT AdaptiveGradientDescentOptimizer:
  public AutomaticGradientDescentOptimizer
{
public:
  /** Standard class typedefs. */
  typedef AdaptiveGradientDescentOptimizer  Self;
  typedef AutomaticGradientDescentOptimizer Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AdaptiveGradientDescentOptimizer, AutomaticGradientDescentOptimizer);

  /** Advance one step following the gradient direction. */
  void AdvanceOneStep(void);

  /** Set/Get the maximum of the sigmoid.
  * Should be >0. Default: 1.0 */
  itkSetMacro(SigmoidMax, double);
  itkGetConstMacro(SigmoidMax, double);

  /** Set/Get the maximum of the sigmoid.
  * Should be <0. Default: -0.8 */
  itkSetMacro(SigmoidMin, double);
  itkGetConstMacro(SigmoidMin, double);

  /** Set/Get the scaling of the sigmoid width. Large values
  * cause a more wide sigmoid. Default: 1e-8. Should be >0. */
  itkSetMacro(SigmoidScale, double);
  itkGetConstMacro(SigmoidScale, double);

protected:

  /** The gradient in the previous step */
  DerivativeType  m_PreviousGradient;

  /** The current time, which serves as input for Compute_a */
  double          m_CurrentTime;

  /** Function to update the current time
  * This function just increments the CurrentTime by 1.
  * Inheriting functions may implement something smarter,
  * for example, dependent on the progress */
  virtual void UpdateCurrentTime( void );

  /** Function to compute the learning rate. */
  double EstimateLearningRate() const;

  AdaptiveGradientDescentOptimizer();
  virtual ~AdaptiveGradientDescentOptimizer() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  AdaptiveGradientDescentOptimizer(const Self &); //purposely not implemented
  void operator=(const Self &);           //purposely not implemented

  /**Parameters for adaptive stochastic gradient descent */
  double                        m_Param_a;
  double                        m_Param_A;
  double                        m_Param_alpha;

  double                        m_InitialTime;
  double                        m_SigmoidMax;
  double                        m_SigmoidMin;
  double                        m_SigmoidScale;

};
} // end namespace itk

#endif
