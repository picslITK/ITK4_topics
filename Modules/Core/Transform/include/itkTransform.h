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
#ifndef __itkTransform_h
#define __itkTransform_h

#include "itkTransformBase.h"
#include "itkVector.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkDiffusionTensor3D.h"
#include "itkVariableLengthVector.h"
#include "vnl/vnl_vector_fixed.h"
#include "itkMatrix.h"
#include "itkIndex.h"
#include "itkImage.h"

namespace itk
{
/** \class Transform
 * \brief Transform points and vectors from an input space to an output space.
 *
 * This abstract class defines the generic interface for a geometric
 * transformation from one space to another. The class provides methods
 * for mapping points, vectors and covariant vectors from the input space
 * to the output space.
 *
 * Given that transformations are not necessarily invertible, this basic
 * class does not provide the methods for back transformation. Back transform
 * methods are implemented in derived classes where appropriate.
 *
 * \par Registration Framework Support
 * Typically a Transform class has several methods for setting its
 * parameters. For use in the registration framework, the parameters must
 * also be represented by an array of doubles to allow communication
 * with generic optimizers. The Array of transformation parameters is set using
 * the SetParameters() method.
 *
 * Another requirement of the registration framework is the computation
 * of the transform Jacobian. In general, an ImageToImageMetric requires
 * the knowledge of the Jacobian in order to compute the metric derivatives.
 * The Jacobian is a matrix whose element are the partial derivatives
 * of the output point with respect to the array of parameters that defines
 * the transform.
 *
 * Subclasses must provide implementations for:
 *   virtual OutputPointType           TransformPoint(const InputPointType  &) const
 *   virtual OutputVectorType          TransformVector(const InputVectorType &) const
 *   virtual OutputVnlVectorType       TransformVector(const InputVnlVectorType &) const
 *   virtual OutputCovariantVectorType TransformCovariantVector(const InputCovariantVectorType &) const
 *   virtual void                      SetParameters(const ParametersType &)
 *   virtual void                      SetFixedParameters(const ParametersType &)
 *   virtual void                      ComputeJacobianWithRespectToParameters(
 *                                                             const InputPointType &,
 *                                                             JacobianType &) const
 *   virtual void                      ComputeJacobianWithRespectToPosition(
 *                                                             const InputPointType & x,
 *                                                             JacobianType &jacobian ) const;
 *
 * Since TranformVector and TransformCovariantVector have multiple
 * overloaded methods from the base class, subclasses must specify:
 *  using Superclass::TransformVector;
 *  using Superclass::TransformCovariantVector;
 *
 *
 * \ingroup ITKTransform
 */
template <class TScalarType,
          unsigned int NInputDimensions = 3,
          unsigned int NOutputDimensions = 3>
class ITK_EXPORT Transform : public TransformBase
{
public:
  /** Standard class typedefs. */
  typedef Transform                Self;
  typedef TransformBase            Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(Transform, TransformBase);

  /** Dimension of the domain space. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, NInputDimensions);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, NOutputDimensions);

  /** Get the size of the input space */
  unsigned int GetInputSpaceDimension(void) const
  {
    return NInputDimensions;
  }

  /** Get the size of the output space */
  unsigned int GetOutputSpaceDimension(void) const
  {
    return NOutputDimensions;
  }

  /** Type of the scalar representing coordinate and vector elements. */
  typedef  TScalarType ScalarType;

  /** Type of the input parameters. */
  typedef  typename Superclass::ParametersType      ParametersType;
  typedef  typename Superclass::ParametersValueType ParametersValueType;
  typedef  Array<ParametersValueType>               DerivativeType;

  /** Type of the Jacobian matrix. */
  typedef  Array2D<ParametersValueType> JacobianType;

  /** Standard vector type for this class. */
  typedef Vector<TScalarType, NInputDimensions>  InputVectorType;
  typedef Vector<TScalarType, NOutputDimensions> OutputVectorType;

  /** Standard variable length vector type for this class
   *  this provides an interface for the VectorImage class */
  typedef VariableLengthVector<TScalarType> InputVectorPixelType;
  typedef VariableLengthVector<TScalarType> OutputVectorPixelType;

  /* Standard tensor type for this class */
  typedef DiffusionTensor3D<TScalarType> InputDiffusionTensor3DType;
  typedef DiffusionTensor3D<TScalarType> OutputDiffusionTensor3DType;

  /** Standard covariant vector type for this class */
  typedef CovariantVector<TScalarType, NInputDimensions>
  InputCovariantVectorType;
  typedef CovariantVector<TScalarType, NOutputDimensions>
  OutputCovariantVectorType;

  /** Standard vnl_vector type for this class. */
  typedef vnl_vector_fixed<TScalarType, NInputDimensions> InputVnlVectorType;
  typedef vnl_vector_fixed<TScalarType, NOutputDimensions>
  OutputVnlVectorType;

  /** Standard coordinate point type for this class */
  typedef Point<TScalarType, NInputDimensions>  InputPointType;
  typedef Point<TScalarType, NOutputDimensions> OutputPointType;

  /** Input index type for this class */
  typedef Index< NInputDimensions > InputIndexType;

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost. */
  typedef Transform<
    TScalarType, NOutputDimensions, NInputDimensions> InverseTransformBaseType;

  typedef typename InverseTransformBaseType::Pointer
  InverseTransformBasePointer;

  typedef Matrix<TScalarType,
                 itkGetStaticConstMacro(OutputSpaceDimension),
                 itkGetStaticConstMacro(InputSpaceDimension)>     MatrixType;

  typedef Matrix<double,
                 itkGetStaticConstMacro(OutputSpaceDimension),
                 itkGetStaticConstMacro(OutputSpaceDimension)>
  OutputDirectionMatrix;
  typedef Matrix<double,
                 itkGetStaticConstMacro(InputSpaceDimension),
                 itkGetStaticConstMacro(InputSpaceDimension)>
  InputDirectionMatrix;
  typedef Matrix<double,
                 itkGetStaticConstMacro(OutputSpaceDimension),
                 itkGetStaticConstMacro(InputSpaceDimension)>
  DirectionChangeMatrix;

  typedef Superclass::NumberOfParametersType    NumberOfParametersType;

  /** Typedefs for use in CanUseTransformIndex. */
  typedef Image< TScalarType, NInputDimensions >         ImageType;
  typedef typename ImageType::RegionType                 RegionType;
  typedef typename ImageType::SpacingType                SpacingType;
  typedef typename ImageType::PointType                  OriginType;
  typedef typename ImageType::DirectionType              DirectionType;

#if 0
  // this method is currently undocummented, untested and broken when input and output dimensions are
  // not the same
  void SetDirectionChange( const OutputDirectionMatrix fixedDir, const InputDirectionMatrix  movingDir );
#endif

  void SetDirectionChangeMatrix( const DirectionChangeMatrix & changeDir )
  {
    m_DirectionChange = changeDir; this->Modified();
  }

  DirectionChangeMatrix GetDirectionChangeMatrix( void ) const
  {
    return m_DirectionChange;
  }

  /**  Method to transform a point.
   * \warning This method must be thread-safe. See, e.g., its use
   * in ResampleImageFilter.
   */
  virtual OutputPointType TransformPoint(const InputPointType  &) const = 0;

  /** Method to transform a point given its index, returning a point.
   * Useful only for certain transforms, e.g. DisplacementFieldTransform.
   * This method is provided for speed optimization, to avoid cost of
   * interpolation when a dense transform's domain is aligned with
   * the input domain.
   * The user must first check if TransformIndex can be used for a given
   * region, spacing, direction and origin, by calling \c CanUseTransformIndex.
   * This method throws an exception unless overridden by a dervied class.
   */
  virtual OutputPointType TransformIndex(const InputIndexType &) const
  { itkExceptionMacro("TransformAtIndex not implemented for class "
                      << this->GetNameOfClass() ); }

  /** Check if the \c TransformIndex method can be used, given a
   * region, spacing, direction and origin. If the passed parameters
   * match those of the transform's displacement field image (or other dense
   * parameter source), then return true. Transforms that cannot support
   * TransformIndex at all will return false by default.
   * Must be overridden by derived classes that can support TransformIndex
   * given the proper conditions. */
  bool CanUseTransformIndex( RegionType&, OriginType&,
                             SpacingType&, DirectionType& ) const
  { return false; }

  /**  Method to transform a vector. */
  virtual OutputVectorType  TransformVector(const InputVectorType &) const = 0;

  /** Method to transform a vector at a given location.
   * For global transforms, \c point is ignored and \c TransformVector( vector )
   * is called. Local transforms (e.g. displacement
   * field transform) must override and provide required behavior. */
  virtual OutputVectorType    TransformVector(
    const InputVectorType & vector,
    const InputPointType & itkNotUsed(point) ) const
  {
    return TransformVector( vector );
  }

  /**  Method to transform a vnl_vector. */
  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &)
  const = 0;

  /** Method to transform a vnl_vector, at a point.
   * For global transforms, \c point is ignored and \c TransformVector( vector )
   * is called. Local transforms (e.g. displacement
   * field transform) must override and provide required behavior. */
  virtual OutputVnlVectorType TransformVector(
    const InputVnlVectorType & vector,
    const InputPointType & itkNotUsed(point) ) const
  {
    return TransformVector( vector );
  }

  /** Method to transform a vector stored in a VectorImage.  */
  virtual OutputVectorPixelType TransformVector(
    const InputVectorPixelType & itkNotUsed(vector) ) const
  {
    itkExceptionMacro( "TransformVector( const InputVectorPixelType & ) is "
                       "unimplemented for " << this->GetNameOfClass() );
  }

  /** Method to transform a vector stored in a VectorImage, at a point.
   * For global transforms, \c point is ignored and \c TransformVector( vector )
   * is called. Local transforms (e.g. displacement
   * field transform) must override and provide required behavior. */
  virtual OutputVectorPixelType TransformVector(
    const InputVectorPixelType & vector,
    const InputPointType & itkNotUsed(point) ) const
  {
    return TransformVector( vector );
  }

  /**  Method to transform a CovariantVector. */
  virtual OutputCovariantVectorType TransformCovariantVector(const InputCovariantVectorType &) const = 0;

  /** Method to transform a CovariantVector, using a point. Global transforms
   * can ignore the \c point parameter. Local transforms (e.g. displacement
   * field transform) must override and provide required behavior.
   * By default, \c point is ignored and
   * \c TransformCovariantVector(vector) is called */
  virtual OutputCovariantVectorType TransformCovariantVector(
    const InputCovariantVectorType & vector,
    const InputPointType & itkNotUsed(point) )
  const
  {
    return TransformCovariantVector( vector );
  }

  /**  Method to transform a CovariantVector stored in a VectorImage. */
  virtual OutputVectorPixelType TransformCovariantVector(
    const InputVectorPixelType & itkNotUsed(vector) ) const
  {
    itkExceptionMacro( "TransformCovariantVector(const InputVectorPixelType &)"
                       "is unimplemented for " << this->GetNameOfClass() );
  }

  /** Method to transform a CovariantVector, using a point. Global transforms
   * can ignore the \c point parameter. Local transforms (e.g. displacement
   * field transform) must override and provide required behavior.
   * By default, \c point is ignored and \c TransformCovariantVector(vector) is
   * called */
  virtual OutputVectorPixelType TransformCovariantVector(
    const InputVectorPixelType & vector,
    const InputPointType & itkNotUsed(point) )
  const
  {
    return TransformCovariantVector( vector );
  }

  /** Method to transform a diffusion tensor */
  virtual OutputDiffusionTensor3DType TransformDiffusionTensor(
    const InputDiffusionTensor3DType & itkNotUsed(tensor) )
  const
  {
    itkExceptionMacro(
      "TransformDiffusionTensor( const InputDiffusionTensor3DType & ) is "
      "unimplemented for " << this->GetNameOfClass() );
  }

  /** Method to transform a diffusion tensor at a point. Global transforms
   * can ignore the \c point parameter. Local transforms (e.g. displacement
   * field transform) must override and provide required behavior.
   * By default, \c point is ignored and \c TransformDiffusionTensor(tensor) is
   * called */
  virtual OutputDiffusionTensor3DType TransformDiffusionTensor(
    const InputDiffusionTensor3DType & itkNotUsed(tensor),
    const InputPointType & itkNotUsed(point) ) const
  {
    itkExceptionMacro(
      "TransformDiffusionTensor( const InputDiffusionTensor3DType &, const "
      "InputPointType & ) is unimplemented for " << this->GetNameOfClass() );
  }

  /** Method to transform a diffusion tensor stored in a VectorImage */
  virtual OutputVectorPixelType TransformDiffusionTensor(
    const InputVectorPixelType & itkNotUsed(tensor) ) const
  {
    itkExceptionMacro(
      "TransformDiffusionTensor( const InputVectorPixelType & ) is "
      "unimplemented for " << this->GetNameOfClass() );
  }

  /** Method to transform a diffusion tensor stored in a VectorImage, at
   * a point.  Global transforms
   * can ignore the \c point parameter. Local transforms (e.g. displacement
   * field transform) must override and provide required behavior.
   * By default, \c point is ignored and \c TransformDiffusionTensor(tensor) is
   * called */
  virtual OutputVectorPixelType TransformDiffusionTensor(
    const InputVectorPixelType & itkNotUsed(tensor),
    const InputPointType & itkNotUsed(point) ) const
  {
    itkExceptionMacro(
      "TransformDiffusionTensor( const InputVectorPixelType &, const "
      "InputPointType & ) is unimplemented for " << this->GetNameOfClass() );
  }

  /** Set the transformation parameters and update internal transformation.
   * SetParameters gives the transform the option to set it's
   * parameters by keeping a reference to the parameters, or by
   * copying.  To force the transform to copy it's parameters call
   * SetParametersByValue.
   * \sa SetParametersByValue
   */
  virtual void SetParameters(const ParametersType &) = 0;

  /** Set the transformation parameters and update internal transformation.
   * This method forces the transform to copy the parameters.  The
   * default implementation is to call SetParameters.  This call must
   * be overridden if the transform normally implements SetParameters
   * by keeping a reference to the parameters.
   * \sa SetParameters
   */
  virtual void SetParametersByValue(const ParametersType & p)
  {
    this->SetParameters(p);
  }

  /** Get the Transformation Parameters. */
  virtual const ParametersType & GetParameters(void) const
  {
    return m_Parameters;
  }

  /** Set the fixed parameters and update internal transformation. */
  virtual void SetFixedParameters(const ParametersType &) = 0;

  /** Get the Fixed Parameters. */
  virtual const ParametersType & GetFixedParameters(void) const
  {
    return m_FixedParameters;
  }

#ifdef ITKV3_COMPATIBILITY
  /**
   * This function is only here for ITKv3 backwards compatibility.
   *
   * This is not a thread-safe version for GetJacobian(), because the internal
   * class member variable m_Jacobian could be changed for different values
   * in different threads.
   *
   * All derived classes should move the computations of computing a jacobian
   * from GetJacobian to ComputeJacobianWithRespectToParameters and then
   * use this forwarding function for backwards compatibility.
   */
  virtual const JacobianType & GetJacobian(const InputPointType  & x) const
  {
    this->ComputeJacobianWithRespectToParameters(x, m_SharedLocalJacobian);
    return m_SharedLocalJacobian;
  }

#endif

  /**
   * Compute the Jacobian of the transformation
   *
   * This method computes the Jacobian matrix of the transformation
   * at a given input point. The rank of the Jacobian will also indicate
   * if the transform is invertible at this point.
   *
   * The Jacobian is be expressed as a matrix of partial derivatives of the
   * output point components with respect to the parameters that defined
   * the transform:
   *
   * \f[
   *
  J=\left[ \begin{array}{cccc}
  \frac{\partial x_{1}}{\partial p_{1}} &
  \frac{\partial x_{1}}{\partial p_{2}} &
  \cdots  & \frac{\partial x_{1}}{\partial p_{m}}\\
  \frac{\partial x_{2}}{\partial p_{1}} &
  \frac{\partial x_{2}}{\partial p_{2}} &
  \cdots  & \frac{\partial x_{2}}{\partial p_{m}}\\
  \vdots  & \vdots  & \ddots  & \vdots \\
  \frac{\partial x_{n}}{\partial p_{1}} &
  \frac{\partial x_{n}}{\partial p_{2}} &
  \cdots  & \frac{\partial x_{n}}{\partial p_{m}}
  \end{array}\right]
   *
   * \f]
   *
   *  This is also used for efficient computation of a point-local jacobian
   *  for dense transforms.
   *  \c jacobian is assumed to be thread-local variable, otherwise memory corruption
   *  will most likely occur during multi-threading.
   *  To avoid repeatitive memory allocation, pass in 'jacobian' with its size
   *  already set. */
  virtual void ComputeJacobianWithRespectToParameters(const InputPointType  & p, JacobianType & jacobian) const = 0;

  /** This provides the ability to get a local jacobian value
   *  in a dense/local transform, e.g. DisplacementFieldTransform. For such
   *  transforms it would be unclear what parameters would refer to.
   *  Generally, global transforms should return an indentity jacobian
   *  since there is no change with respect to position. */
  virtual void ComputeJacobianWithRespectToPosition(const InputPointType & x, JacobianType & jacobian ) const = 0;

  /** Update the transform's parameters by the adding values in \c update
   * to current parameter values.
   * We assume \c update is of the same length as Parameters. Throw
   * exception otherwise.
   * \c factor is a scalar multiplier for each value in update.
   * SetParameters is called at the end of this method, to allow transforms
   * to perform any required operations on the update parameters, typically
   * a converion to member variables for use in TransformPoint.
   * Derived classes should override to provide specialized behavior.
   */
  virtual void UpdateTransformParameters( DerivativeType & update, TScalarType factor = 1.0 );

  /** Return the number of local parameters that completely defines the
   *  Transform at an individual voxel.
   *  For transforms with local support, this will enable downstream
   *  computation of the jacobian wrt only the local support region.
   *  For instance, in the case of a displacement field, this will be equal to
   *  the number of image dimensions. If it is an affine transform, this will
   *  be the same as the GetNumberOfParameters().
   */
  virtual NumberOfParametersType GetNumberOfLocalParameters(void) const
  {
    return this->GetNumberOfParameters();
  }

  /** Return the number of parameters that completely define the Transfom  */
  virtual NumberOfParametersType GetNumberOfParameters(void) const
  {
    return this->m_Parameters.Size();
  }

  /** Returns a boolean indicating whether it is possible or not to compute the
   * inverse of this current Transform. If it is possible, then the inverse of
   * the transform is returned in the inverseTransform variable passed by the
   * user.  The inverse is recomputed if this current transform has been
   * modified.
   * This method is intended to be overriden as needed by derived classes.
   *
   */
  bool GetInverse( Self *itkNotUsed(inverseTransform) ) const
  {
    return false;
  }

  /** Return an inverse of this transform. If the inverse has not been
   *  implemented, return NULL. The type of the inverse transform
   *  does not necessarily need to match the type of the forward
   *  transform. This allows one to return a numeric inverse transform
   *  instead.
   */
  virtual InverseTransformBasePointer GetInverseTransform() const
  {
    return NULL;
  }

  /** Generate a platform independant name */
  virtual std::string GetTransformTypeAsString() const;

  /** Indicates if this transform is linear. A transform is defined to be
   * linear if the transform of a linear combination of points is equal to the
   * linear combination (with the same coefficients) of the individual
   * transforms of each point. The transform T will be linear if given two
   * points P and Q, and scalar coefficients a and b, then
   *
   *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
   *
   * By default, we assume this to NOT be the case for most transforms.
   * However, transforms for which this is true will overload and reimplement
   * this method accordingly.
   *
   * \warning This method must be thread-safe. See, e.g., its use
   * in ResampleImageFilter.
   */
  virtual bool IsLinear() const
  {
    return false;
  }

  /** Indicates if this transform is a "global" transform
   *  e.g. an affine transform, or a local one, e.g. a displacement field.
   */
  virtual bool HasLocalSupport() const
  {
    return false;
  }
protected:
  Transform();
  Transform(NumberOfParametersType NumberOfParameters);
  virtual ~Transform()
  {
  }

  mutable ParametersType m_Parameters;
  mutable ParametersType m_FixedParameters;

#ifdef ITKV3_COMPATIBILITY
  // This is only needed to provide the old interface that returns a reference to the Jacobian.
  // It is NOT thread-safe and should be avoided whenever possible.
  mutable JacobianType m_SharedLocalJacobian;
#endif

  mutable DirectionChangeMatrix m_DirectionChange;
private:
  Transform(const Self &);      // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  template <typename TType>
  std::string GetTransformTypeAsString(TType *) const
  {
    std::string rval("other");

    return rval;
  }

  std::string GetTransformTypeAsString(float *) const
  {
    std::string rval("float");

    return rval;
  }

  std::string GetTransformTypeAsString(double *) const
  {
    std::string rval("double");

    return rval;
  }

};
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_Transform(_, EXPORT, TypeX, TypeY)         \
  namespace itk                                                 \
  {                                                             \
  _( 3 ( class EXPORT Transform<ITK_TEMPLATE_3 TypeX> ) )     \
  namespace Templates                                           \
  {                                                             \
  typedef Transform<ITK_TEMPLATE_3 TypeX> Transform##TypeY; \
  }                                                             \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkTransform.hxx"
#endif

#endif
