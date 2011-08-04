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
#ifndef __itkPolyAffineTransform_h
#define __itkPolyAffineTransform_h

#include <iostream>

#include "itkMatrix.h"
#include "itkTransform.h"
#include "itkMacro.h"

namespace itk
{
/**
 * Matrix and Offset transformation of a vector space (e.g. space coordinates)
 *
 * This class serves as a base class for transforms that can be expressed
 * as a linear transformation plus a constant offset (e.g., affine, similarity
 * and rigid transforms).   This base class also provides the concept of
 * using a center of rotation and a translation instead of an offset.
 *
 * As derived instances of this class are specializations of an affine
 * transform, any two of these transformations may be composed and the result
 * is an affine transformation.  However, the order is important.
 * Given two affine transformations T1 and T2, we will say that
 * "precomposing T1 with T2" yields the transformation which applies
 * T1 to the source, and then applies T2 to that result to obtain the
 * target.  Conversely, we will say that "postcomposing T1 with T2"
 * yields the transformation which applies T2 to the source, and then
 * applies T1 to that result to obtain the target.  (Whether T1 or T2
 * comes first lexicographically depends on whether you choose to
 * write mappings from right-to-left or vice versa; we avoid the whole
 * problem by referring to the order of application rather than the
 * textual order.)
 *
 * There are three template parameters for this class:
 *
 * ScalarT       The type to be used for scalar numeric values.  Either
 *               float or double.
 *
 * NInputDimensions   The number of dimensions of the input vector space.
 *
 * NOutputDimensions   The number of dimensions of the output vector space.
 *
 * This class provides several methods for setting the matrix and offset
 * defining the transform. To support the registration framework, the
 * transform parameters can also be set as an Array<double> of size
 * (NInputDimension + 1) * NOutputDimension using method SetParameters().
 * The first (NOutputDimension x NInputDimension) parameters defines the
 * matrix in row-major order (where the column index varies the fastest).
 * The last NOutputDimension parameters defines the translation
 * in each dimensions.
 *
 * \ingroup Transforms
 *
 */

template< class TScalarType,
          unsigned int NInputDimensions = 3,
          unsigned int NOutputDimensions = 3
>
// Number of dimensions in the output space
class PolyAffineTransform:
  public Transform< TScalarType, NInputDimensions, NOutputDimensions >
{
public:
  /** Standard typedefs   */
  typedef PolyAffineTransform                   Self;
  typedef Transform< TScalarType,
                     NInputDimensions,
                     NOutputDimensions >        Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(PolyAffineTransform, Transform);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Dimension of the domain space. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, NInputDimensions);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, NOutputDimensions);
  //itkStaticConstMacro( ParametersDimension, unsigned int,
  //                     NOutputDimensions * ( NInputDimensions + 1 ) * NTransforms);

  /** Parameters Type   */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::ParametersValueType ParametersValueType;

  /** Jacobian Type   */
  typedef typename Superclass::JacobianType JacobianType;

  /** Standard scalar type for this class */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard vector type for this class   */
  typedef Vector< TScalarType,
                  itkGetStaticConstMacro(InputSpaceDimension) >  InputVectorType;
  typedef Vector< TScalarType,
                  itkGetStaticConstMacro(OutputSpaceDimension) > OutputVectorType;
  typedef typename OutputVectorType::ValueType OutputVectorValueType;

  /** Standard covariant vector type for this class   */
  typedef CovariantVector< TScalarType,
                           itkGetStaticConstMacro(InputSpaceDimension) >
  InputCovariantVectorType;
  typedef CovariantVector< TScalarType,
                           itkGetStaticConstMacro(OutputSpaceDimension) >
  OutputCovariantVectorType;

  typedef typename Superclass::InputVectorPixelType  InputVectorPixelType;
  typedef typename Superclass::OutputVectorPixelType OutputVectorPixelType;

  /** Standard tensor type for this class */
  typedef typename Superclass::InputTensorType  InputTensorType;
  typedef typename Superclass::OutputTensorType OutputTensorType;

  typedef typename Superclass::InputTensorEigenVectorType  InputTensorEigenVectorType;
  typedef typename Superclass::OutputTensorEigenVectorType OutputTensorEigenVectorType;

  /** Standard vnl_vector type for this class   */
  typedef vnl_vector_fixed< TScalarType,
                            itkGetStaticConstMacro(InputSpaceDimension) >
  InputVnlVectorType;
  typedef vnl_vector_fixed< TScalarType,
                            itkGetStaticConstMacro(OutputSpaceDimension) >
  OutputVnlVectorType;

  /** Standard coordinate point type for this class   */
  typedef Point< TScalarType,
                 itkGetStaticConstMacro(InputSpaceDimension) >
  InputPointType;
  typedef typename InputPointType::ValueType InputPointValueType;
  typedef Point< TScalarType,
                 itkGetStaticConstMacro(OutputSpaceDimension) >
  OutputPointType;
  typedef typename OutputPointType::ValueType OutputPointValueType;

  /** Standard matrix type for this class   */
  typedef Matrix< TScalarType, itkGetStaticConstMacro(OutputSpaceDimension),
                  itkGetStaticConstMacro(InputSpaceDimension) >
  MatrixType;
  typedef typename MatrixType::ValueType MatrixValueType;

  /** Standard inverse matrix type for this class   */
  typedef Matrix< TScalarType, itkGetStaticConstMacro(InputSpaceDimension),
                  itkGetStaticConstMacro(OutputSpaceDimension) >
  InverseMatrixType;

  typedef InputPointType CenterType;

  typedef OutputVectorType               OffsetType;
  typedef typename OffsetType::ValueType OffsetValueType;

  typedef OutputVectorType TranslationType;

  typedef typename TranslationType::ValueType TranslationValueType;

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost. */
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  /** Types for atomic transforms */
  typedef MatrixOffsetTransformBase<TScalarType, NInputDimensions,
                                    NOutputDimensions>  AtomTransformType;
  typedef typename AtomTransformType::Pointer           AtomTransformPointer;
  typedef typename AtomTransformType::ParametersType    AtomParametersType;

  /** Transform queue type */
  typedef std::vector<AtomTransformPointer>             AtomTransformSetType;

  AtomTransformSetType & GetAtomTransformSet()
    {
    return m_AtomTransformSet;
    }

  /** Set the transformation to an Identity
   *
   * This sets the matrix to identity and the Offset to null. */
  virtual void SetIdentity(void);

  /** Set the transformation from a container of parameters.
   * The first (NOutputDimension x NInputDimension) parameters define the
   * matrix and the last NOutputDimension parameters the translation.
   * Offset is updated based on current center. */
  void SetParameters(const ParametersType & parameters);

  /** Get the Transformation Parameters. */
  const ParametersType & GetParameters(void) const;

  /** Set the fixed parameters and update internal transformation. */
  virtual void SetFixedParameters(const ParametersType &);

  /** Get the Fixed Parameters. */
  virtual const ParametersType & GetFixedParameters(void) const;

  /** Transform by an affine transformation
   *
   * This method applies the affine transform given by self to a
   * given point or vector, returning the transformed point or
   * vector.  The TransformPoint method transforms its argument as
   * an affine point, whereas the TransformVector method transforms
   * its argument as a vector. */
  OutputPointType       TransformPoint(const InputPointType & point) const;

  OutputVectorType      TransformVector(const InputVectorType & vector) const;

  OutputVectorType      TransformVector(const InputVectorType & vector,
                                        const InputPointType & itkNotUsed(point) ) const
    { return TransformVector( vector ); }

  OutputVnlVectorType   TransformVector(const InputVnlVectorType & vector) const;

  OutputVnlVectorType   TransformVector(const InputVnlVectorType & vector,
                                        const InputPointType & itkNotUsed(point) ) const
    { return TransformVector( vector ); }

  //OutputVectorPixelType TransformVector(const InputVectorPixelType & vector) const;

  //OutputVectorPixelType TransformVector(const InputVectorPixelType & vector,
  //                                      const InputPointType & itkNotUsed(point) ) const
  //  { return TransformVector( vector ); }


  OutputCovariantVectorType TransformCovariantVector(
      const InputCovariantVectorType & vector) const;

  OutputCovariantVectorType TransformCovariantVector(
      const InputCovariantVectorType & vector,
      const InputPointType & itkNotUsed(point) ) const
    { return TransformCovariantVector( vector ); }

  //OutputVectorPixelType TransformCovariantVector(
  //    const InputVectorPixelType & vector) const;

  //OutputVectorPixelType TransformCovariantVector(
  //    const InputVectorPixelType & vector,
  //    const InputPointType & itkNotUsed(point) ) const
  //  { return TransformCovariantVector( vector ); }

  /** Compute the Jacobian of the transformation
   *
   * This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the transform
   * is invertible at this point. */
  const JacobianType & GetJacobian(const InputPointType & point) const;

  /** get local Jacobian for the given point
   * \c j will sized properly as needed.
   * This is a thread-safe version for GetJacobian(). Otherwise,
   * m_Jacobian could be changed for different values in different threads. */
  void GetJacobianWithRespectToParameters(const InputPointType  &x,
                                          JacobianType &j) const;

  /** Create inverse of an affine transformation
   *
   * This populates the parameters an affine transform such that
   * the transform is the inverse of self. If self is not invertible,
   * an exception is thrown.
   * Note that by default the inverese transform is centered at
   * the origin. If you need to compute the inverse centered at a point, p,
   *
   * \code
   * transform2->SetCenter( p );
   * transform1->GetInverse( transform2 );
   * \endcode
   *
   * transform2 will now contain the inverse of transform1 and will
   * with its center set to p. Flipping the two statements will produce an
   * incorrect transform.
   *
   */
  bool GetInverse(Self *inverse) const;

  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const;

  void AddAtomTransform( AtomTransformType *t  )
  {
    AtomTransformPointer tp = t;
    m_AtomTransformSet.push_back(tp);
    this->Modified();
  }

protected:
  /** Construct an PolyAffineTransform object
   *
   * This method constructs a new PolyAffineTransform object and
   * initializes the matrix and offset parts of the transformation
   * to values specified by the caller.  If the arguments are
   * omitted, then the PolyAffineTransform is initialized to an identity
   * transformation in the appropriate number of dimensions. */
  PolyAffineTransform();

  /** Destroy an PolyAffineTransform object */
  virtual ~PolyAffineTransform();

  /** Print contents of an PolyAffineTransform */
  void PrintSelf(std::ostream & s, Indent indent) const;

private:

  PolyAffineTransform(const Self & other);
  const Self & operator=(const Self &);

  AtomTransformSetType            m_AtomTransformSet;

}; //class PolyAffineTransform
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_PolyAffineTransform(_, EXPORT, TypeX, TypeY)                         \
  namespace itk                                                                                 \
  {                                                                                             \
  _( 3 ( class EXPORT PolyAffineTransform< ITK_TEMPLATE_3 TypeX > ) )                     \
  namespace Templates                                                                           \
  {                                                                                             \
  typedef PolyAffineTransform< ITK_TEMPLATE_3 TypeX > PolyAffineTransform##TypeY; \
  }                                                                                             \
  }

#if ITK_TEMPLATE_EXPLICIT
//template < class TScalarType, unsigned int NInputDimensions, unsigned int
// NOutputDimensions>
//   const unsigned int itk::PolyAffineTransform< TScalarType,
// NInputDimensions, NOutputDimensions >::InputSpaceDimension;
//template < class TScalarType, unsigned int NInputDimensions, unsigned int
// NOutputDimensions>
//   const unsigned int itk::PolyAffineTransform< TScalarType,
// NInputDimensions, NOutputDimensions >::OutputSpaceDimension;
//template < class TScalarType, unsigned int NInputDimensions, unsigned int
// NOutputDimensions>
//   const unsigned int itk::PolyAffineTransform< TScalarType,
// NInputDimensions, NOutputDimensions >::ParametersDimension;
#include "Templates/itkPolyAffineTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkPolyAffineTransform.hxx"
#endif

#endif /* __itkPolyAffineTransform_h */
