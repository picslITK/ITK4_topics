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
#include "itkVectorInterpolateImageFunction.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkImageVectorTransformParametersHelper.h"

namespace itk
{

/** \class DeformationFieldTransform
 * \brief TODO
 *
 * \ingroup Transforms
 *
 * \ingroup ITK-Review
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
  typedef typename Superclass::InputTensorType      InputTensorType;
  typedef typename Superclass::OutputTensorType     OutputTensorType;

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
  typedef MatrixOffsetTransformBase< double, NDimensions, NDimensions > AffineTransformType;
  typedef typename AffineTransformType::Pointer AffineTransformPointer;

  typedef MatrixOffsetTransformBase< double, NDimensions, NDimensions>  LocalTransformType;
  typedef typename LocalTransformType::Pointer                          LocalTransformPointer;

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

  /**  Method to transform a point. */
  virtual OutputPointType TransformPoint( const InputPointType& thisPoint ) const;

  /**  Method to transform a vector. */
  virtual OutputVectorType TransformVector(const InputVectorType &) const
    { itkExceptionMacro( "TransformVector(Vector) unimplemented, use TransformVector(Vector,Point)" ); }

  virtual OutputVectorPixelType TransformVector(const InputVectorPixelType &) const
    { itkExceptionMacro( "TransformVector(Vector) unimplemented, use TransformVector(Vector,Point)" ); }

  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &) const
    { itkExceptionMacro( "TransformVector(Vector) unimplemented, use TransformVector(Vector,Point)" ); }

  virtual OutputVectorType TransformVector(const InputVectorType &, const InputPointType & ) const;

  virtual OutputVectorPixelType TransformVector(const InputVectorPixelType &, const InputPointType & ) const;

  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &, const InputPointType & ) const;

  /** Method to transform a tensor */
  OutputTensorType TransformTensor(const InputTensorType & ) const
    { itkExceptionMacro( "TransformTensor(Tensor) unimplemented, use TransformTensor(Tensor,Point)" ); }

  OutputVectorPixelType TransformTensor(const InputVectorPixelType & ) const
    { itkExceptionMacro( "TransformTensor(Tensor) unimplemented, use TransformTensor(Tensor,Point)" ); }

  OutputTensorType TransformTensor(const InputTensorType &, const InputPointType &) const;

  OutputVectorPixelType TransformTensor(const InputVectorPixelType &, const InputPointType &) const;

  /**  Method to transform a CovariantVector. */
  virtual OutputCovariantVectorType TransformCovariantVector( const InputCovariantVectorType &) const
      { itkExceptionMacro( "TransformCovariantVector(CovariantVector) unimplemented, use TransformCovariantVector(CovariantVector,Point)" ); }

  virtual OutputVectorPixelType TransformCovariantVector( const InputVectorPixelType &) const
      { itkExceptionMacro( "TransformCovariantVector(CovariantVector) unimplemented, use TransformCovariantVector(CovariantVector,Point)" ); }

  virtual OutputCovariantVectorType TransformCovariantVector( const InputCovariantVectorType &, const InputPointType & ) const;

  virtual OutputVectorPixelType TransformCovariantVector( const InputVectorPixelType &, const InputPointType & ) const;


  /** Set the transformation parameters. This sets the deformation
   * field image directly. */
  virtual void SetParameters(const ParametersType & params)
    {
    if( &(this->m_Parameters) != &params )
      {
      this->m_Parameters = params;
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
                                                  JacobianType &j) const;

  /**
   * Compute the jacobian with respect to the position, by index.
   * \c j will be resized as needed.
   */
  virtual void GetJacobianWithRespectToPosition(const IndexType  &x,
                                                  JacobianType &j) const;

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

private:
  DeformationFieldTransform( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDeformationFieldTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkDeformationFieldTransform.txx"
#endif

#endif // __itkDeformationFieldTransform_h
