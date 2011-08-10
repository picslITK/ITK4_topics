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
#ifndef __itkDeformationFieldTransform_h
#define __itkDeformationFieldTransform_h

#include "itkTransform.h"

#include "itkImage.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkImageVectorTransformParametersHelper.h"
#include "itkGaussianOperator.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"

namespace itk
{

/* Forward-declaration to avoid including the header file that would
 * introduce a circular dependency in the Transform module. */
template< class TInputImage, class TCoordRep >
class VectorInterpolateImageFunction;

/** \class DeformationFieldTransform
 * \brief Provides local/dense/high-dimensionaltiy transformation via a
 * a deformation field.
 *
 * The deformation field is defined using an itkImage, and must be set
 * before use by the user, using \c SetDeformationField. The image has
 * the same dimensionality as the input and output spaces, defined by
 * template parameter \c NDimensions, and is an image of vectors of
 * type \c OutputVectorType (i.e. with dimensionality NDimensions as well).
 *
 * An interpolator of type \c VectorInterpolateImageFunction is used with
 * the deformation field image. By default, VectorLinearInterpolateImageFunction
 * is used, and the user can override using SetInterpolator.
 *
 * The deformation field data is stored using the common \c TransformParameters
 * type
 * in conjunction with the \c ImageVectorTransformParametersHelper class. This
 * allows access of the deformation field image as if it were an itkArray,
 * allowing transparent use with other classes. Note that the \c SetParameters
 * method
 * will copy the passed parameters, which can be costly since deformation fields
 * are dense and thus potentially very large.
 *
 * UpdateTransformParameters - after adding the passed update, this method
 * performs smoothing using the GaussianOperator filter. The smoothing can
 * be performed on demand by calling \c SmoothDeformationFieldGauss - see
 * this method for details.
 *
 * Because this is a local transform, methods that have a version that takes
 * a point must be used, such as \c TransformVector, \c TransformCovariantVector,
 * and \c TransformDiffusionTensor. Also, \c GetJacobianWithRespectToParameters simply
 * returns
 * an identity matrix, and \c GetJacobianWithRespectToPosition should be used.
 *
 * \ingroup Transforms
 *
 * \ingroup ITKTransform
 */
template
  <class TScalar, unsigned int NDimensions>
class ITK_EXPORT DeformationFieldTransform :
  public Transform<TScalar, NDimensions, NDimensions>
{
public:
  /** Standard class typedefs. */
  typedef DeformationFieldTransform                         Self;
  typedef Transform<TScalar, NDimensions, NDimensions>      Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( DeformationFieldTransform, Transform );

  /** New macro for creation of through a Smart Pointer */
  itkSimpleNewMacro( Self );

  /** Leave CreateAnother undefined. To fully implement here, it must be
   * sure to copy all members. It may be called from transform-cloning
   * that only copies parameters, so override here to prevent
   * its use without copying full members. */
  virtual::itk::LightObject::Pointer CreateAnother(void) const
    {
    itkExceptionMacro("CreateAnother unimplemented. See source comments.");
    }

  /** InverseTransform type. */
  typedef typename Superclass::InverseTransformBasePointer  InverseTransformBasePointer;

  /** Scalar type. */
  typedef typename Superclass::ScalarType  ScalarType;

  /** Type of the input parameters. */
  typedef  typename Superclass::ParametersType      ParametersType;

  /** Jacobian type. */
  typedef typename Superclass::JacobianType  JacobianType;

  /** Standard coordinate point type for this class. */
  typedef typename Superclass::InputPointType   InputPointType;
  typedef typename Superclass::OutputPointType  OutputPointType;

  /** Standard vector type for this class. */
  typedef typename Superclass::InputVectorType      InputVectorType;
  typedef typename Superclass::OutputVectorType     OutputVectorType;

  typedef typename Superclass::InputVectorPixelType   InputVectorPixelType;
  typedef typename Superclass::OutputVectorPixelType  OutputVectorPixelType;

  /** Standard covariant vector type for this class */
  typedef typename Superclass::InputCovariantVectorType
    InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType
    OutputCovariantVectorType;

  /** Standard vnl_vector type for this class. */
  typedef typename Superclass::InputVnlVectorType   InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType  OutputVnlVectorType;

  /** Standard tensor type for this class */
  typedef typename Superclass::InputDiffusionTensorType      InputDiffusionTensorType;
  typedef typename Superclass::OutputDiffusionTensorType     OutputDiffusionTensorType;

  /** Standard tensor type for this class */
  typedef typename Superclass::InputTensorEigenVectorType      InputTensorEigenVectorType;
  typedef typename Superclass::OutputTensorEigenVectorType     OutputTensorEigenVectorType;

    /** Standard tensor type for this class */
  typedef typename Superclass::InputTensorMatrixType      InputTensorMatrixType;
  typedef typename Superclass::OutputTensorMatrixType     OutputTensorMatrixType;

  /** Derivative type */
  typedef typename Superclass::DerivativeType       DerivativeType;

  /** Dimension of the domain spaces. */
  itkStaticConstMacro( Dimension, unsigned int, NDimensions );

  /** Define the deformation field type and corresponding interpolator type. */
  typedef Image<OutputVectorType,
    itkGetStaticConstMacro( Dimension )> DeformationFieldType;
  typedef VectorInterpolateImageFunction
    <DeformationFieldType, ScalarType> InterpolatorType;

  /** Standard Index type for Deformation Field */
  typedef typename DeformationFieldType::IndexType IndexType;

  /* Define tranform based upon ImageDirections of Deformation Field */
  typedef MatrixOffsetTransformBase< double, NDimensions, NDimensions >
                                                          AffineTransformType;
  typedef typename AffineTransformType::Pointer AffineTransformPointer;
  typedef Superclass                            TransformType;

  /** Define the internal parameter helper used to access the field */
  typedef ImageVectorTransformParametersHelper<
                                          ScalarType,
                                          OutputVectorType::Dimension,
                                          itkGetStaticConstMacro( Dimension ) >
                                                TransformParametersHelperType;

  /** Get/Set the deformation field. */
  itkGetObjectMacro( DeformationField, DeformationFieldType );
  /** Set the deformation field. Create special set accessor to update
   * interpolator and assign deformation field to transform parameters
   * container. */
  virtual void SetDeformationField( DeformationFieldType* field );

  /** Get/Set the inverse deformation field. */
  itkGetObjectMacro( InverseDeformationField, DeformationFieldType );
  itkSetObjectMacro( InverseDeformationField, DeformationFieldType );

  /** Get/Set the interpolator. */
  itkGetObjectMacro( Interpolator, InterpolatorType );
  /* Create out own set accessor that assigns the deformation field */
  virtual void SetInterpolator( InterpolatorType* interpolator );

  /** Get/Set the GaussianOperator variance */
  itkSetMacro( GaussianSmoothSigma, ScalarType );
  itkGetConstReferenceMacro( GaussianSmoothSigma, ScalarType );

  /** Get the modification time of deformation field */
  itkGetConstReferenceMacro( DeformationFieldSetTime, unsigned long );

  /**  Method to transform a point. */
  virtual OutputPointType TransformPoint( const InputPointType& thisPoint )
                                                                        const;

  /**  Method to transform a vector. */
  virtual OutputVectorType TransformVector(const InputVectorType &) const
  { itkExceptionMacro( "TransformVector(Vector) unimplemented, use "
    "TransformVector(Vector,Point)" ); }

  virtual OutputVectorPixelType TransformVector(const InputVectorPixelType &)
                                                                          const
  { itkExceptionMacro( "TransformVector(Vector) unimplemented, use "
    "TransformVector(Vector,Point)" ); }

  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &) const
  { itkExceptionMacro( "TransformVector(Vector) unimplemented, use "
  "TransformVector(Vector,Point)" ); }

  virtual OutputVectorType TransformVector(const InputVectorType &,
                                           const InputPointType & ) const;

  virtual OutputVectorPixelType TransformVector(const InputVectorPixelType &,
                                                const InputPointType & ) const;

  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &,
                                              const InputPointType & ) const;

  /** Method to transform a tensor */
  OutputDiffusionTensorType TransformDiffusionTensor(
                                      const InputDiffusionTensorType & ) const
  { itkExceptionMacro( "TransformDiffusionTensor(Tensor) unimplemented, use "
    "TransformDiffusionTensor(Tensor,Point)" ); }

  OutputVectorPixelType TransformDiffusionTensor(const InputVectorPixelType & )
                                                                          const
  { itkExceptionMacro( "TransformDiffusionTensor(Tensor) unimplemented, use "
    "TransformDiffusionTensor(Tensor,Point)" ); }

  OutputDiffusionTensorType TransformDiffusionTensor(
                                              const InputDiffusionTensorType &,
                                              const InputPointType &) const;

  OutputVectorPixelType TransformDiffusionTensor(const InputVectorPixelType &,
                                                 const InputPointType &) const;

  /**  Method to transform a CovariantVector. */
  virtual OutputCovariantVectorType TransformCovariantVector(
                                        const InputCovariantVectorType &) const
  { itkExceptionMacro( "TransformCovariantVector(CovariantVector) "
    "unimplemented, use TransformCovariantVector(CovariantVector,Point)" ); }

  virtual OutputVectorPixelType TransformCovariantVector(
                                            const InputVectorPixelType &) const
  { itkExceptionMacro( "TransformCovariantVector(CovariantVector) "
    "unimplemented, use TransformCovariantVector(CovariantVector,Point)" ); }

  /** Transform a CovariantVector of type InputCovariantVectorType, at point. */
  virtual OutputCovariantVectorType TransformCovariantVector(
                        const InputCovariantVectorType &,
                        const InputPointType &) const;

  /** Transform a CovariantVector of type InputVectorPixelType, at point. */
  virtual OutputVectorPixelType TransformCovariantVector(
                        const InputVectorPixelType &,
                        const InputPointType & ) const;

  /** Set the transformation parameters. This sets the deformation
   * field image directly. */
  virtual void SetParameters(const ParametersType & params)
    {
    if( &(this->m_Parameters) != &params )
      {
      if( params.Size() != this->m_Parameters.Size() )
        {
        itkExceptionMacro("Input parameters size (" << params.Size()
                          << ") does not match internal size ("
                          << this->m_Parameters.Size() << ").");
        }
      /* copy into existing object */
      this->m_Parameters = params;
      this->Modified();
      }
    }

  /** Set the fixed parameters and update internal transformation. */
  virtual void SetFixedParameters(const ParametersType &)
  { itkExceptionMacro("SetFixedParameters unimplemented."); }

  /** Get the Fixed Parameters. */
  virtual const ParametersType & GetFixedParameters(void) const
  { itkExceptionMacro("GetFixedParameters unimplemented."); }

  /**
   * Compute the jacobian with respect to the parameters.
   * Returns identity.
   */
  virtual JacobianType & GetJacobian( const InputPointType & ) const;

  virtual void GetJacobianWithRespectToParameters(const InputPointType  &x,
                                                  JacobianType &j) const
  { j = this->m_IdentityJacobian; }

  virtual void GetJacobianWithRespectToParameters(const IndexType  &x,
                                                  JacobianType &j) const
  { j = this->m_IdentityJacobian; }

  /**
   * Compute the jacobian with respect to the position, by point.
   * \c j will be resized as needed.
   */
  virtual void GetJacobianWithRespectToPosition(const InputPointType  &x,
                                                JacobianType &j ) const;

  /**
   * Compute the jacobian with respect to the position, by index.
   * \c j will be resized as needed.
   */
  virtual void GetJacobianWithRespectToPosition(const IndexType  &x,
                                                JacobianType &j ) const;

  /**
   * Compute the inverse jacobian of the forward deformation field with
   * respect to the position, by point. Note that this is different than
   * the jacobian of the inverse deformation field. This takes advantage
   * of the ability to compute the inverse jacobian of a deformation field
   * by simply reversing the sign of the forward jacobian.
   * \c useSVD ... TODO */
  virtual void GetInverseJacobianOfForwardFieldWithRespectToPosition(
                                  const InputPointType & point,
                                  JacobianType & jacobian,
                                  bool useSVD = false )
                                                                         const;

  /**
   * Compute the inverse jacobian of the forward deformation field with
   * respect to the position, by index.Note that this is different than
   * the jacobian of the inverse deformation field. This takes advantage
   * of the ability to compute the inverse jacobian of a deformation field
   * by simply reversing the sign of the forward jacobian.
   * \c useSVD ... TODO */
  virtual void GetInverseJacobianOfForwardFieldWithRespectToPosition(
                                  const IndexType & index,
                                  JacobianType & jacobian,
                                  bool useSVD = false )
                                                                        const;

  /** Update the transform's parameters by the values in \c update.
   * We assume \c update is of the same length as Parameters. Throw
   * exception otherwise.
   * \c factor is a scalar multiplier for each value in update.
   * \c SmoothDeformationFieldGaussian is called after the update is
   * added to the field.
   * See base class for more details.
   */
  virtual void UpdateTransformParameters( DerivativeType & update,
                                          ScalarType factor = 1.0 );

  /** Smooth the deformation field in-place.
   * Uses m_GaussSmoothSigma to change the variance for the GaussianOperator,
   * the default is 3. Other parameters MaximumError and MaximumKernelWidth
   * are hardcoded.
   * \warning Not thread safe. Does its own threading. */
  virtual void SmoothDeformationFieldGauss();

  /** Return an inverse of this transform. */
  bool GetInverse( Self *inverse ) const;

  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const;

  /** This transform is not linear. */
  virtual bool IsLinear() const { return false; }

  virtual unsigned int GetNumberOfLocalParameters(void) const
  { return Dimension; }

  virtual bool HasLocalSupport() const { return true; }

protected:
  DeformationFieldTransform();
  virtual ~DeformationFieldTransform();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** The deformation field and its inverse (if it exists). */
  typename DeformationFieldType::Pointer      m_DeformationField;
  typename DeformationFieldType::Pointer      m_InverseDeformationField;

  /** The interpolator. */
  typename InterpolatorType::Pointer          m_Interpolator;

  /** Used in SmoothDeformationFieldGauss to as variance for the
   * GaussianOperator */
  ScalarType                                  m_GaussianSmoothSigma;

  /** Track when the deformation field was last set/assiend, as
   * distinct from when it may have had its contents modified. */
  unsigned long                             m_DeformationFieldSetTime;

private:
  DeformationFieldTransform( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented


  /** Internal method for calculating either forward or inverse jacobian,
   * depending on state of \c doInverseJacobian. Used by
   * public methods \c GetJacobianWithRespectToPosition and
   * \c GetInverseJacobianOfForwardFieldWithRespectToPosition to
   * perform actual work.
   * \c doInverseJacobian ... TODO */
  virtual void GetJacobianWithRespectToPositionInternal(
                                  const IndexType & index,
                                  JacobianType & jacobian,
                                  bool doInverseJacobian)
                                                                        const;

  /** Used to holder temporary deformation field during smoothing.
   * Use member variable to avoid allocation on stack. */
  typename DeformationFieldType::Pointer    m_SmoothGaussTempField;

  /** Track when the temporary deformation field used during smoothing
   * was last modified/initialized. We only want to change it if the
   * main deformation field is also changed, i.e. assigned to a new object */
  unsigned long                             m_SmoothGaussTempFieldModifiedTime;

  /** Type of Gaussian Operator used during smoothing. Define here
   * so we can use a member var during the operation. */
  typedef GaussianOperator<ScalarType,Dimension>      SmoothGaussOperatorType;
  typedef VectorNeighborhoodOperatorImageFilter< DeformationFieldType,
                                                 DeformationFieldType >
                                                      SmoothGaussSmootherType;
  SmoothGaussOperatorType                     m_SmoothGaussOperator;
  typename SmoothGaussSmootherType::Pointer   m_SmoothGaussSmoother;
};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDeformationFieldTransform+-.h"
#endif

#if ITK_MANUAL_INSTANTIATION
# include "itkDeformationFieldTransform.hxx"
#endif

#endif // __itkDeformationFieldTransform_h
