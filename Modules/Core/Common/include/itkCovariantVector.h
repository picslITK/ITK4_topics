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
#ifndef __itkCovariantVector_h
#define __itkCovariantVector_h

#include "vnl/vnl_vector_ref.h"
#include "itkIndent.h"
#include "itkVector.h"

namespace itk
{
/** \class CovariantVector
 * \brief A templated class holding a n-Dimensional covariant vector.
 *
 * CovariantVector is a templated class that holds a single vector
 * (i.e., an array of values).  CovariantVector can be used as the data
 * type held at each pixel in an Image or at each vertex of an Mesh.
 * The template parameter T can be any data type that behaves like a
 * primitive (or atomic) data type (int, short, float, complex).
 * The NVectorDimension defines the number of components in the vector array.
 *
 * CovariantVector is not a dynamically extendible array like std::vector. It is
 * intended to be used like a mathematical vector.
 *
 * If you wish a simpler pixel types, you can use Scalar, which represents
 * a single data value at a pixel. There is also the more complex type
 * ScalarCovariantVector, which supports (for a given pixel)
 * a single scalar value plus an array of vector values.
 * (The scalar and vectors can be of different data type.)
 *
 * CovariantVector is the type that should be used for representing normals
 * to surfaces and gradients of functions. AffineTransform transform
 * covariant vectors different than vectors.
 *
 * \ingroup Geometry
 * \ingroup DataRepresentation
 *
 * \sa Image
 * \sa Mesh
 * \sa Point
 * \sa Vector
 * \sa Matrix
 * \ingroup ITK-Common
 * \wikiexample{Math/CovariantVector,Create a covariant vector}
 * \wikiexample{Math/CovariantVectorNorm,Compute the norm of a covariant vector}
 */

template< class T, unsigned int NVectorDimension = 3 >
class ITK_EXPORT CovariantVector:public FixedArray< T, NVectorDimension >
{
public:
  /** Standard class typedefs. */
  typedef CovariantVector                   Self;
  typedef FixedArray< T, NVectorDimension > Superclass;

  /** ValueType can be used to declare a variable that is the same type
   * as a data element held in an CovariantVector.   */
  typedef T                                             ValueType;
  typedef typename NumericTraits< ValueType >::RealType RealValueType;

  /** Component value type */
  typedef T ComponentType;

  /** Dimension of the Space */
  itkStaticConstMacro(Dimension, unsigned int, NVectorDimension);

  /** I am a covariant vector. */
  typedef Self CovariantVectorType;

  /** The Array type from which this CovariantVector is derived. */
  typedef FixedArray< T, NVectorDimension > BaseArray;

  /** Get the dimension (size) of the vector. */
  static unsigned int GetCovariantVectorDimension()
  { return NVectorDimension; }

  /** Set a vnl_vector_ref referencing the same memory block. */
  void SetVnlVector(const vnl_vector< T > &);

  /** Get a vnl_vector_ref referencing the same memory block. */
  vnl_vector_ref< T > GetVnlVector(void);

  /** Get a vnl_vector with a copy of the internal memory block. */
  vnl_vector< T > GetVnlVector(void) const;

  /** Set a vnl_vector_ref referencing the same memory block.
   * \deprecated Use SetVnlVector() instead. */
  void Set_vnl_vector(const vnl_vector< T > &);

  /** Get a vnl_vector_ref referencing the same memory block.
   * \deprecated Use GetVnlVector() instead. */
  vnl_vector_ref< T > Get_vnl_vector(void);

  /** Get a vnl_vector with a copy of the internal memory block.
   * \deprecated Use GetVnlVector() instead. */
  vnl_vector< T > Get_vnl_vector(void) const;

  /** Default constructor and copy constructors. */
  CovariantVector():BaseArray() {}
  CovariantVector(const ValueType & r);

  /** Pass-through constructor for the Array base class. Implicit casting is
   * performed to initialize constructor from any another one of datatype. */
  template< class TVectorValueType >
  CovariantVector(const CovariantVector< TVectorValueType,
                                         NVectorDimension > & r):BaseArray(r) {}
  CovariantVector(const ValueType r[Dimension]):BaseArray(r) {}

  /** Assignment operator with implicit casting from another data type */
  template< class Tt >
  Self & operator=(const Tt & v)
  {
    BaseArray::operator=(v);
    return *this;
  }

  /** Pass-through assignment operator for the Array base class. */
  CovariantVector & operator=(const Self & r);

  CovariantVector & operator=(const ValueType r[NVectorDimension]);

  /** Scalar operator*=.  Scales elements by a scalar. */
  template< class Tt >
  inline const Self & operator*=(const Tt & value)
  {
    for ( unsigned int i = 0; i < NVectorDimension; i++ )
      {
      ( *this )[i] = static_cast< ValueType >( ( *this )[i] * value );
      }
    return *this;
  }

  /** Scalar operator/=.  Scales (divides) elements by a scalar. */
  template< class Tt >
  const Self & operator/=(const Tt & value)
  {
    for ( unsigned int i = 0; i < NVectorDimension; i++ )
      {
      ( *this )[i] = static_cast< ValueType >( ( *this )[i] / value );
      }
    return *this;
  }

  /** CovariantVector operator+=.  Adds a vectors to the current vector. */
  const Self & operator+=(const Self & vec);

  /** CovariantVector operator-=.  Subtracts a vector from a current vector. */
  const Self & operator-=(const Self & vec);

  /** CovariantVector negation.  Negate all the elements of a vector.
   *  Return a new vector */
  Self operator-() const;

  /** CovariantVector addition. Add two vectors. Return a new vector. */
  Self operator+(const Self & vec) const;

  /** CovariantVector subtraction. Subtract two vectors. Return a new vector. */
  Self operator-(const Self & vec) const;

  /** CovariantVector operator*.
   * Performs the inner product of two covariant vectors.
   * \warning This is equivalent to the scalar product only if the reference
   * system has orthogonal axis and equal scales.  */
  ValueType operator *(const Self & vec) const;

  /** operator*.  Performs the scalar product with a vector (contravariant).
   * This scalar product is invariant under affine transformations */
  ValueType operator *(const Vector< T, NVectorDimension > & vec) const;

  /** Scalar operator*. Scale the elements of a vector by a scalar.
   * Return a new vector. */
  inline Self operator*(const ValueType & val) const
  {
    Self result;

    for ( unsigned int i = 0; i < NVectorDimension; i++ )
      {
      result[i] = static_cast< ValueType >( ( *this )[i] * val );
      }
    return result;
  }

  /** Scalar operator/. Scale (divide) the elements of a vector by a scalar.
   * Return a new vector. */
  template< class Tt >
  inline Self operator/(const Tt & val) const
  {
    Self result;

    for ( unsigned int i = 0; i < NVectorDimension; i++ )
      {
      result[i] = static_cast< ValueType >( ( *this )[i] / val );
      }
    return result;
  }

  /** Returns the Euclidean Norm of the vector  */
  RealValueType GetNorm(void) const;

  /** Returns the number of components in this vector type */
  static unsigned int GetNumberOfComponents() { return NVectorDimension; }

  /** Divides the covariant vector componets by the norm */
  void Normalize(void);

  /** Returns vector's Squared Euclidean Norm  */
  RealValueType GetSquaredNorm(void) const;

  /** Copy from another CovariantVector with a different representation type.
   *  Casting is done with C-Like rules  */
  template< typename TCoordRepB >
  void CastFrom(const CovariantVector< TCoordRepB, NVectorDimension > & pa)
  {
    for ( unsigned int i = 0; i < NVectorDimension; i++ )
      {
      ( *this )[i] = static_cast< T >( pa[i] );
      }
  }
};

/** Premultiply Operator for product of a vector and a scalar.
 *  CovariantVector< T, N >  =  T * CovariantVector< T,N > */
template< class T, unsigned int NVectorDimension >
inline
CovariantVector< T, NVectorDimension >
operator*(const T & scalar, const CovariantVector< T, NVectorDimension > & v)
{
  return v * scalar;
}

ITKCommon_EXPORT void CrossProduct(CovariantVector< double, 3 > &,
                                   const Vector< double, 3 > &,
                                   const Vector< double, 3 > &);

ITKCommon_EXPORT void CrossProduct(CovariantVector< float, 3 > &,
                                   const Vector< float, 3 > &,
                                   const Vector< float, 3 > &);

ITKCommon_EXPORT void CrossProduct(CovariantVector< int, 3 >,
                                   const Vector< int, 3 > &,
                                   const Vector< int, 3 > &);
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_CovariantVector(_, EXPORT, TypeX, TypeY)     \
  namespace itk                                                   \
  {                                                               \
  _( 2 ( class EXPORT CovariantVector< ITK_TEMPLATE_2 TypeX > ) ) \
  namespace Templates                                             \
  {                                                               \
  typedef CovariantVector< ITK_TEMPLATE_2 TypeX >                 \
  CovariantVector##TypeY;                                       \
  }                                                               \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkCovariantVector+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkNumericTraitsCovariantVectorPixel.h"
#include "itkCovariantVector.txx"
#endif

#endif
