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
#ifndef __itkCenteredAffineTransform_h
#define __itkCenteredAffineTransform_h

#include "itkAffineTransform.h"

namespace itk
{
/** \class CenteredAffineTransform
 * \brief Affine transformation with a specified center of rotation.
 *
 * This class implements an Affine transform in which the rotation center
 * can be explicitly selected.
 *
 * \ingroup ITKTransform
 */
template <
  class TScalarType = double,      // Data type for scalars
  unsigned int NDimensions = 3>
// Number of dimensions in the input space
class ITK_EXPORT CenteredAffineTransform : public AffineTransform<TScalarType,
                                                                  NDimensions>
{
public:
  /** Standard typedefs   */
  typedef CenteredAffineTransform                   Self;
  typedef AffineTransform<TScalarType, NDimensions> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(CenteredAffineTransform, AffineTransform);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Dimension of the domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro( ParametersDimension, unsigned int,
                       NDimensions * ( NDimensions + 2 ) );

  /** Types taken from the Superclass */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::JacobianType        JacobianType;
  typedef typename Superclass::ScalarType          ScalarType;
  typedef typename Superclass::InputVectorType     InputVectorType;
  typedef typename Superclass::OutputVectorType    OutputVectorType;
  typedef typename Superclass::InputCovariantVectorType
  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType
  OutputCovariantVectorType;

  typedef typename Superclass::InputVnlVectorType    InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType   OutputVnlVectorType;
  typedef typename Superclass::InputPointType        InputPointType;
  typedef typename Superclass::InputPointValueType   InputPointValueType;
  typedef typename Superclass::OutputVectorValueType OutputVectorValueType;
  typedef typename Superclass::OutputPointType       OutputPointType;
  typedef typename Superclass::MatrixType            MatrixType;
  typedef typename Superclass::MatrixValueType       MatrixValueType;
  typedef typename Superclass::OffsetType            OffsetType;

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost. */
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  /** Set/Get the transformation from a container of parameters.
   * The first (NDimension x NDimension) parameters define the
   * matrix, the next N parameters define the center of rotation
   * and the last N parameters define the translation to be applied
   * after the coordinate system has been restored to the rotation center.
   * Note that the Offset of the superclass is no longer in the
   * parameters array since it is fully dependent on the rotation
   * center and the translation parameters. */
  void SetParameters(const ParametersType & parameters);

  const ParametersType & GetParameters(void) const;

  /** Compute the Jacobian of the transformation
   *
   * This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the transform
   * is invertible at this point. */
  virtual void ComputeJacobianWithRespectToParameters( const InputPointType  & p, JacobianType & jacobian) const;

  /** Get an inverse of this transform. */
  bool GetInverse(Self *inverse) const;

  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const;

protected:
  /** Construct an CenteredAffineTransform object */
  CenteredAffineTransform();

  /** Destroy an CenteredAffineTransform object */
  virtual ~CenteredAffineTransform();
private:
  CenteredAffineTransform(const Self & other);
  const Self & operator=(const Self &);

}; // class CenteredAffineTransform
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_CenteredAffineTransform(_, EXPORT, TypeX, TypeY)                       \
  namespace itk                                                                             \
  {                                                                                         \
  _( 2 ( class EXPORT CenteredAffineTransform<ITK_TEMPLATE_2 TypeX> ) )                   \
  namespace Templates                                                                       \
  {                                                                                         \
  typedef CenteredAffineTransform<ITK_TEMPLATE_2 TypeX> CenteredAffineTransform##TypeY; \
  }                                                                                         \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkCenteredAffineTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkCenteredAffineTransform.hxx"
#endif

#endif /* __itkCenteredAffineTransform_h */
