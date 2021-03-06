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
#ifndef __itkShear2DTransform_h
#define __itkShear2DTransform_h

#include <iostream>
#include "itkTransform.h"
#include "itkMacro.h"
#include "itkMatrix.h"
#include "itkMatrixOffsetTransformBase.h"

namespace itk
{
/** \class ScaleTransform
 * \brief Scale transformation of a vector space (e.g. space coordinates)
 *
 * The same functionality could be obtained by using the Affine tranform,
 * but with a large difference in performace since the affine transform will
 * use a matrix multiplication using a diagonal matrix.
 *
 * \ingroup Transforms
 * \ingroup ITKTransform
 */
template<
class TScalarType = float, // Type for cordinate representation type (float or
// double)
unsigned int NDimensions = 3  >
// Number of dimensions
// class ITK_EXPORT Shear2DTransform:public Transform< TScalarType,
class ITK_EXPORT Shear2DTransform:public MatrixOffsetTransformBase< TScalarType,
NDimensions,
NDimensions >
{
public:
    /** Standard class typedefs.   */
    typedef Shear2DTransform                                     Self;
    // typedef Transform< TScalarType, NDimensions, NDimensions > Superclass;
    typedef MatrixOffsetTransformBase< TScalarType, NDimensions, NDimensions > Superclass;
    typedef SmartPointer< Self >                               Pointer;
    typedef SmartPointer< const Self >                         ConstPointer;

    /** New macro for creation of through a smart pointer. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(Shear2DTransform, Transform);

    /** Dimension of the domain space. */
    itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
    itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);

    /** Scalar type. */
    typedef typename Superclass::ScalarType ScalarType;

    /** Parameters type. */
    typedef typename Superclass::ParametersType ParametersType;

    /** Jacobian type. */
    typedef typename Superclass::JacobianType JacobianType;

    /** Standard vector type for this class. */
    // typedef FixedArray< TScalarType, NDimensions > ScaleType;
    // typedef TScalarType ScalarType;

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

    typedef typename Superclass::MatrixType      MatrixType;

    /** Set parameters.  This method sets the parameters for the transform value
     *  specified by the user. The parameters are organized as scale[i] =
     *  parameter[i]. That means that in 3D the scale parameters for the coordinates
     *  {x,y,z} are {parameter[0], parameter[1], parameter[2]} respectively */
    void SetParameters(const ParametersType & parameters);

    /** Get the parameters that uniquely define the transform This is typically
     * used by optimizers during the process of image registration.  The parameters
     * are organized as {scale X, scale Y, scale Z } = { parameter[0],
     * parameter[1], parameter[2] } respectively */
    const ParametersType & GetParameters(void) const;

    /** Set the fixed parameters and update internal
     * transformation. This transform has no fixed paramaters
     */
    virtual void SetFixedParameters(const ParametersType &) {}

    /** Get the fixed parameters */
    virtual const ParametersType & GetFixedParameters(void) const;

    /** Get the Jacobian matrix. */
    virtual const JacobianType & GetJacobian(const InputPointType & point) const;

    virtual void ComputeJacobianWithRespectToParameters(const InputPointType &point, JacobianType &j) const;

    /** Set the factors of an Scale Transform
     * This method sets the factors of an ScaleTransform to a
     * value specified by the user.
     * This method cannot be done with SetMacro because itk::Array has not an
     * operator== defined. The array of scales correspond in order to the factors
     * to be applied to each one of the coordinaates. For example, in 3D,
     * scale[0] corresponds to X, scale[1] corresponds to Y and scale[2]
     * corresponds to Z. */
    void SetK(const ScalarType & K)
    {  this->m_K = K; this->ComputeMatrix(); this->Modified();}

    virtual void ComputeMatrix(void);

    /** Compose with another ScaleTransform. */
    void Compose(const Self *other, bool pre = false);

    //  /** Compose this transform transformation with another scaling.
    //   * The pre argument is irrelevant here since scale transforms are commutative,
    //   * pre and postcomposition are therefore equivalent. */
    //  void Scale(const ScaleType & scale, bool pre = false);

    /** Transform by a scale transformation
     * This method applies the scale transform given by self to a
     * given point or vector, returning the transformed point or
     * vector. */
    OutputPointType     TransformPoint(const InputPointType  & point) const;

    OutputVectorType    TransformVector(const InputVectorType & vector) const;

    OutputVnlVectorType TransformVector(const InputVnlVectorType & vector) const;

   OutputCovariantVectorType TransformCovariantVector(
        const InputCovariantVectorType & vector) const;

    //  /** Back transform by a scale transformation
    //   * This method finds the point or vector that maps to a given
    //   * point or vector under the scale transformation defined by
    //   * self.  If no such point exists, an exception is thrown. */
    //  inline InputPointType     BackTransform(const OutputPointType  & point) const;
    //
    //  inline InputVectorType    BackTransform(const OutputVectorType & vector) const;
    //
    //  inline InputVnlVectorType BackTransform(const OutputVnlVectorType & vector) const;
    //
    //  inline InputCovariantVectorType BackTransform(
    //    const OutputCovariantVectorType & vector) const;

    /** Find inverse of a scale transformation
     * This method creates and returns a new ScaleTransform object
     * which is the inverse of self.  If self is not invertible,
     * false is returned. */
    bool GetInverse(Self *inverse) const;

    /** Return an inverse of this transform. */
    virtual InverseTransformBasePointer GetInverseTransform() const;

    /** Set the transformation to an Identity
     *
     * This sets all the scales to 1.0 */
    void SetIdentity(void)
    { this->m_K = 0; }
    // { m_Scale.Fill(1.0); }

    /** Set/Get the center used as fixed point for the scaling */
    itkSetMacro(Center, InputPointType);
    itkGetConstReferenceMacro(Center, InputPointType);

    /** Get access to scale values */
    // itkGetConstReferenceMacro(Scale, ScaleType);

    /** Indicates that this transform is linear. That is, given two
     * points P and Q, and scalar coefficients a and b, then
     *
     *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
     */
    virtual bool IsLinear() const { return true; }


protected:

    /** Construct an ScaleTransform object. */
    Shear2DTransform();

    /** Destroy an ScaleTransform object. */
    ~Shear2DTransform();

    /** Print contents of an ScaleTransform */
    void PrintSelf(std::ostream & os, Indent indent) const;

private:
    Shear2DTransform(const Self & other);   //purposely not implemented
    const Self & operator=(const Self &); //purposely not implemented

    // ScaleType m_Scale;    // Scales of the transformation
    ScalarType m_K;

    InputPointType         m_Center; // Scaling center
    mutable ParametersType m_FixedParameters;
};                         //class ScaleTransform

//// Back transform a point
//template< class ScalarType, unsigned int NDimensions >
//inline
//typename Shear2DTransform< ScalarType, NDimensions >::InputPointType
//Shear2DTransform< ScalarType, NDimensions >::BackTransform(const OutputPointType & point) const
//{
//  InputPointType result;
//
//  for ( unsigned int i = 0; i < SpaceDimension; i++ )
//    {
//    result[i] = ( point[i] + m_Center[i] ) / m_Scale[i] - m_Center[i];
//    }
//  return result;
//}
//
//// Back transform a vector
//template< class ScalarType, unsigned int NDimensions >
//inline
//typename Shear2DTransform< ScalarType, NDimensions >::InputVectorType
//Shear2DTransform< ScalarType, NDimensions >::BackTransform(const OutputVectorType & vect) const
//{
//  InputVectorType result;
//
//  for ( unsigned int i = 0; i < SpaceDimension; i++ )
//    {
//    result[i] = vect[i] / m_Scale[i];
//    }
//  return result;
//}
//
//// Back transform a vnl_vector
//template< class ScalarType, unsigned int NDimensions >
//inline
//typename Shear2DTransform< ScalarType, NDimensions >::InputVnlVectorType
//Shear2DTransform< ScalarType, NDimensions >::BackTransform(const OutputVnlVectorType & vect) const
//{
//  InputVnlVectorType result;
//
//  for ( unsigned int i = 0; i < SpaceDimension; i++ )
//    {
//    result[i] = vect[i] / m_Scale[i];
//    }
//  return result;
//}
//
//// Back Transform a CovariantVector
//template< class ScalarType, unsigned int NDimensions >
//inline
//typename Shear2DTransform< ScalarType, NDimensions >::InputCovariantVectorType
//Shear2DTransform< ScalarType, NDimensions >::BackTransform(const OutputCovariantVectorType & vect) const
//{
//  // Covariant Vectors are scaled by the inverse
//  InputCovariantVectorType result;
//
//  for ( unsigned int i = 0; i < SpaceDimension; i++ )
//    {
//    result[i] = vect[i] * m_Scale[i];
//    }
//  return result;
//}
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_Shear2DTransform(_, EXPORT, TypeX, TypeY)              \
        namespace itk                                                           \
        {                                                                       \
    _( 2 ( class EXPORT Shear2DTransform< ITK_TEMPLATE_2 TypeX > ) )          \
    namespace Templates                                                     \
    {                                                                       \
            typedef Shear2DTransform< ITK_TEMPLATE_2 TypeX > Shear2DTransform##TypeY; \
    }                                                                       \
        }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkShear2DTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkShear2DTransform.hxx"
#endif

#endif /* __itkShear2DTransform_h */
