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
#ifndef __itkEuler2DTransform_h
#define __itkEuler2DTransform_h

#include <iostream>
#include "itkRigid2DTransform.h"

namespace itk
{
/** \class Euler2DTransform
 *
 * \brief Euler2DTransform of a vector space (e.g. space coordinates)
 *
 * This transform applies a rigid transformation is 2D space.
 * The transform is specified as a rotation around arbitrary center
 * and is followed by a translation.
 *
 * This transform is basically is a synonym for Rigid2DTransform.
 *
 * \sa Rigid2DTransform
 *
 * \ingroup Transforms
 * \ingroup ITK-Transform
 */
template< class TScalarType = double >
// Data type for scalars
class ITK_EXPORT Euler2DTransform:
  public Rigid2DTransform< TScalarType >
{
public:
  /** Standard class typedefs. */
  typedef Euler2DTransform                Self;
  typedef Rigid2DTransform< TScalarType > Superclass;
  typedef SmartPointer< Self >            Pointer;
  typedef SmartPointer< const Self >      ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Euler2DTransform, Rigid2DTransform);

  /** Dimension of parameters. */
  itkStaticConstMacro(SpaceDimension, unsigned int, 2);
  itkStaticConstMacro(ParametersDimension, unsigned int, 3);

  /** Scalar type. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Parameters type. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Jacobian type. */
  typedef typename Superclass::JacobianType JacobianType;

  /** Point type. */
  typedef typename Superclass::InputPointType  InputPointType;
  typedef typename Superclass::OutputPointType OutputPointType;

  /** Vector type. */
  typedef typename Superclass::InputVectorType  InputVectorType;
  typedef typename Superclass::OutputVectorType OutputVectorType;

  /** CovariantVector type. */
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;

  /** VnlVector type. */
  typedef typename Superclass::InputVnlVectorType  InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType OutputVnlVectorType;
  typedef typename Superclass::MatrixType          MatrixType;

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost.*/
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  /**
   * This method creates and returns a new Euler2DTransform object
   * which is the inverse of self.
   */
  void CloneInverseTo(Pointer & newinverse) const;

  /** Get an inverse of this transform. */
  bool GetInverse(Self *inverse) const;

  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const;

  /**
   * This method creates and returns a new Euler2DTransform object
   * which has the same parameters as self.
   */
  void CloneTo(Pointer & clone) const;

  /**
   * Update the angle from the underlying matrix. This method
   * is old and is retained for backward compatibility.
   */
  void ComputeAngleFromMatrix()
  { this->ComputeMatrixParameters(); }
protected:
  Euler2DTransform();
  ~Euler2DTransform(){}

  Euler2DTransform(unsigned int outputSpaceDimension, unsigned int parametersDimension);

  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  Euler2DTransform(const Self &); //purposely not implemented
  void operator=(const Self &);   //purposely not implemented
};                                //class Euler2DTransform
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_Euler2DTransform(_, EXPORT, TypeX, TypeY)     \
  namespace itk                                                    \
  {                                                                \
  _( 1 ( class EXPORT Euler2DTransform< ITK_TEMPLATE_1 TypeX > ) ) \
  namespace Templates                                              \
  {                                                                \
  typedef Euler2DTransform< ITK_TEMPLATE_1 TypeX >                 \
  Euler2DTransform##TypeY;                                       \
  }                                                                \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkEuler2DTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkEuler2DTransform.txx"
#endif

#endif /* __itkEuler2DTransform_h */
