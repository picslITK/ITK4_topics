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
#ifndef __itkVersorRigid3DTransform_h
#define __itkVersorRigid3DTransform_h

#include <iostream>
#include "itkVersorTransform.h"

namespace itk
{
/** \class VersorRigid3DTransform
 *
 * \brief VersorRigid3DTransform of a vector space (e.g. space coordinates)
 *
 * This transform applies a rotation and translation to the space
 * The parameters for this transform can be set either using individual Set
 * methods or in serialized form using SetParameters() and SetFixedParameters().
 *
 * The serialization of the optimizable parameters is an array of 6 elements.
 * The first 3 elements are the components of the versor representation
 * of 3D rotation. The last 3 parameters defines the translation in each
 * dimension.
 *
 * The serialization of the fixed parameters is an array of 3 elements defining
 * the center of rotation.
 *
 *
 * \ingroup Transforms
 * \ingroup ITK-Transform
 */
template< class TScalarType = double >
//Data type for scalars (float or double)
class ITK_EXPORT VersorRigid3DTransform:
  public VersorTransform< TScalarType >
{
public:
  /** Standard class typedefs. */
  typedef VersorRigid3DTransform         Self;
  typedef VersorTransform< TScalarType > Superclass;
  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VersorRigid3DTransform, VersorTransform);

  /** Dimension of parameters. */
  itkStaticConstMacro(SpaceDimension, unsigned int, 3);
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(ParametersDimension, unsigned int, 6);

  /** Parameters Type   */
  typedef typename Superclass::ParametersType            ParametersType;
  typedef typename Superclass::JacobianType              JacobianType;
  typedef typename Superclass::ScalarType                ScalarType;
  typedef typename Superclass::InputPointType            InputPointType;
  typedef typename Superclass::OutputPointType           OutputPointType;
  typedef typename Superclass::InputVectorType           InputVectorType;
  typedef typename Superclass::OutputVectorType          OutputVectorType;
  typedef typename Superclass::InputVnlVectorType        InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType       OutputVnlVectorType;
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;
  typedef typename Superclass::MatrixType                MatrixType;
  typedef typename Superclass::InverseMatrixType         InverseMatrixType;
  typedef typename Superclass::CenterType                CenterType;
  typedef typename Superclass::OffsetType                OffsetType;
  typedef typename Superclass::TranslationType           TranslationType;

  /** Versor type. */
  typedef typename Superclass::VersorType VersorType;
  typedef typename Superclass::AxisType   AxisType;
  typedef typename Superclass::AngleType  AngleType;

  typedef typename Superclass::AxisValueType        AxisValueType;
  typedef typename Superclass::TranslationValueType TranslationValueType;
  typedef typename Superclass::ParameterValueType   ParameterValueType;

  /** Set the transformation from a container of parameters
   * This is typically used by optimizers.
   * There are 6 parameters. The first three represent the
   * versor, the last three represent the translation. */
  void SetParameters(const ParametersType & parameters);

  virtual const ParametersType & GetParameters(void) const;

  /** This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the
   * transform is invertible at this point. */
  const JacobianType & GetJacobian(const InputPointType  & point) const;

protected:
  VersorRigid3DTransform(unsigned int outputSpaceDim,
                         unsigned int paramDim);
  VersorRigid3DTransform(const MatrixType & matrix,
                         const OutputVectorType & offset);
  VersorRigid3DTransform();
  ~VersorRigid3DTransform(){}

  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  VersorRigid3DTransform(const Self &); //purposely not implemented
  void operator=(const Self &);         //purposely not implemented
};                                      //class VersorRigid3DTransform
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_VersorRigid3DTransform(_, EXPORT, TypeX, TypeY)                      \
  namespace itk                                                                           \
  {                                                                                       \
  _( 1 ( class EXPORT VersorRigid3DTransform< ITK_TEMPLATE_1 TypeX > ) )                  \
  namespace Templates                                                                     \
  {                                                                                       \
  typedef VersorRigid3DTransform< ITK_TEMPLATE_1 TypeX > VersorRigid3DTransform##TypeY; \
  }                                                                                       \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkVersorRigid3DTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkVersorRigid3DTransform.txx"
#endif

#endif /* __itkVersorRigid3DTransform_h */
