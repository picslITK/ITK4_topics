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
#ifndef __itkTranslationTransform_h
#define __itkTranslationTransform_h

#include <iostream>
#include "itkTransform.h"
#include "itkMacro.h"
#include "itkMatrix.h"

namespace itk
{
/** \class TranslationTransform
 * \brief Translation transformation of a vector space (e.g. space coordinates)
 *
 * The same functionality could be obtained by using the Affine tranform,
 * but with a large difference in performace.
 *
 * \ingroup Transforms
 * \ingroup ITKTransform
 *
 * \wiki
 * \wikiexample{SimpleOperations/TranslationTransform,Translate an image}
 * \wikiexample{VectorImages/VectorResampleImageFilter,Translate a vector image}
 * \wikiexample{Registration/ImageRegistrationMethod,A basic global registration of two images}
 * \wikiexample{Registration/MutualInformation,Mutual Information}
 * \endwiki
 */
template<
  class TScalarType = double,          // Data type for scalars (float or
                                       // double)
  unsigned int NDimensions = 3 >
// Number of dimensions
class ITK_EXPORT TranslationTransform:
  public Transform< TScalarType, NDimensions, NDimensions >
{
public:
  /** Standard class typedefs. */
  typedef TranslationTransform                               Self;
  typedef Transform< TScalarType, NDimensions, NDimensions > Superclass;
  typedef SmartPointer< Self >                               Pointer;
  typedef SmartPointer< const Self >                         ConstPointer;

  /** New macro for creation of through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TranslationTransform, Transform);

  /** Dimension of the domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);

  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard parameters container. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Standard Jacobian container. */
  typedef typename Superclass::JacobianType JacobianType;

  /** Standard vector type for this class. */
  typedef Vector< TScalarType, NDimensions > InputVectorType;
  typedef Vector< TScalarType, NDimensions > OutputVectorType;

  /** Standard covariant vector type for this class. */
  typedef CovariantVector< TScalarType, NDimensions > InputCovariantVectorType;
  typedef CovariantVector< TScalarType, NDimensions > OutputCovariantVectorType;

  /** Standard vnl_vector type for this class. */
  typedef vnl_vector_fixed< TScalarType, NDimensions > InputVnlVectorType;
  typedef vnl_vector_fixed< TScalarType, NDimensions > OutputVnlVectorType;

  /** Standard coordinate point type for this class. */
  typedef Point< TScalarType, NDimensions > InputPointType;
  typedef Point< TScalarType, NDimensions > OutputPointType;

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost.*/
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  /** This method returns the value of the offset of the
   * TranslationTransform. */
  const OutputVectorType & GetOffset(void) const
  { return m_Offset; }

  /** This method sets the parameters for the transform
   * value specified by the user. */
  void SetParameters(const ParametersType & parameters);

  /** Get the Transformation Parameters. */
  virtual const ParametersType & GetParameters(void) const;

  /** Set offset of an Translation Transform.
   * This method sets the offset of an TranslationTransform to a
   * value specified by the user. */
  void SetOffset(const OutputVectorType & offset)
  { m_Offset = offset; return; }

  /** Compose with another TranslationTransform. */
  void Compose(const Self *other, bool pre = 0);

  /** Compose affine transformation with a translation.
   * This method modifies self to include a translation of the
   * origin.  The translation is precomposed with self if pre is
   * true, and postcomposed otherwise. */
  void Translate(const OutputVectorType & offset, bool pre = 0);

  /** Transform by an affine transformation.
   * This method applies the affine transform given by self to a
   * given point or vector, returning the transformed point or
   * vector. */
  OutputPointType     TransformPoint(const InputPointType  & point) const;

  OutputVectorType    TransformVector(const InputVectorType & vector) const;

  OutputVnlVectorType TransformVector(const InputVnlVectorType & vector) const;

  OutputCovariantVectorType TransformCovariantVector(
    const InputCovariantVectorType & vector) const;

  /** This method finds the point or vector that maps to a given
   * point or vector under the affine transformation defined by
   * self.  If no such point exists, an exception is thrown. */
  inline InputPointType    BackTransform(const OutputPointType  & point) const;

  inline InputVectorType   BackTransform(const OutputVectorType & vector) const;

  inline InputVnlVectorType BackTransform(const OutputVnlVectorType & vector) const;

  inline InputCovariantVectorType BackTransform(
    const OutputCovariantVectorType & vector) const;

  /** Find inverse of an affine transformation.
   * This method creates and returns a new TranslationTransform object
   * which is the inverse of self.  If self is not invertible,
   * false is returned.  */
  bool GetInverse(Self *inverse) const;

  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const;

  /** Compute the Jacobian Matrix of the transformation at one point */
  virtual const JacobianType & GetJacobian(const InputPointType  & point) const;

  /** Compute the Jacobian Matrix of the transformation at one point */
  virtual void GetJacobianWithRespectToParameters(const InputPointType  & point,
          JacobianType &j) const;

  /** Get the jacobian with respect to position, which simply is an identity
   *  jacobian because the transform is position-invariant.
   *  \jac will be resized as needed, but it will be more efficient if
   *  it is already properly sized. */
  virtual void GetJacobianWithRespectToPosition(const InputPointType  &x,
                                                  JacobianType &jac) const;

  /** Set the parameters to the IdentityTransform */
  void SetIdentity(void);

  /** Return the number of parameters that completely define the Transfom  */
  virtual unsigned int GetNumberOfParameters(void) const
  { return NDimensions; }

  /** Indicates that this transform is linear. That is, given two
   * points P and Q, and scalar coefficients a and b, then
   *
   *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
   */
  virtual bool IsLinear() const { return true; }

  /** Set the fixed parameters and update internal transformation.
   * The Translation Transform does not require fixed parameters,
   * therefore the implementation of this method is a null operation. */
  virtual void SetFixedParameters(const ParametersType &)
  {}

  /** Get the Fixed Parameters. The TranslationTransform does not
   * require Fixed parameters, therefore this method returns an
   * parameters array of size zero. */
  virtual const ParametersType & GetFixedParameters(void) const
  {
    this->m_FixedParameters.SetSize(0);
    return this->m_FixedParameters;
  }

protected:
  TranslationTransform();
  ~TranslationTransform();
  /** Print contents of an TranslationTransform. */
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  TranslationTransform(const Self &); //purposely not implemented
  void operator=(const Self &);       //purposely not implemented

  OutputVectorType m_Offset; // Offset of the transformation
};                           //class TranslationTransform

// Back transform a point
template< class TScalarType, unsigned int NDimensions >
inline
typename TranslationTransform< TScalarType, NDimensions >::InputPointType
TranslationTransform< TScalarType, NDimensions >::BackTransform(const OutputPointType & point) const
{
  return point - m_Offset;
}

// Back transform a vector
template< class TScalarType, unsigned int NDimensions >
inline
typename TranslationTransform< TScalarType, NDimensions >::InputVectorType
TranslationTransform< TScalarType, NDimensions >::BackTransform(const OutputVectorType & vect) const
{
  return vect;
}

// Back transform a vnl_vector
template< class TScalarType, unsigned int NDimensions >
inline
typename TranslationTransform< TScalarType, NDimensions >::InputVnlVectorType
TranslationTransform< TScalarType, NDimensions >::BackTransform(const OutputVnlVectorType & vect) const
{
  return vect;
}

// Back Transform a CovariantVector
template< class TScalarType, unsigned int NDimensions >
inline
typename TranslationTransform< TScalarType, NDimensions >::InputCovariantVectorType
TranslationTransform< TScalarType, NDimensions >::BackTransform(const OutputCovariantVectorType & vect) const
{
  return vect;
}
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_TranslationTransform(_, EXPORT, TypeX, TypeY)                    \
  namespace itk                                                                       \
  {                                                                                   \
  _( 2 ( class EXPORT TranslationTransform< ITK_TEMPLATE_2 TypeX > ) )                \
  namespace Templates                                                                 \
  {                                                                                   \
  typedef TranslationTransform< ITK_TEMPLATE_2 TypeX > TranslationTransform##TypeY; \
  }                                                                                   \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkTranslationTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkTranslationTransform.hxx"
#endif

#endif /* __itkTranslationTransform_h */
