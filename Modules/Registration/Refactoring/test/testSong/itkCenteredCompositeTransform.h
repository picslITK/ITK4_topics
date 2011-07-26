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
#ifndef __itkCenteredCompositeTransform_h
#define __itkCenteredCompositeTransform_h


#include "itkCompositeTransform.h"
// #include "itkTransform.h"

#include <deque>

namespace itk
{

/** \class CenteredCompositeTransform
 * \brief This class contains a list of transforms and concatenates them by
 * composition.
 *
 * The fixed parameter is the center of the transform, which will NOT be affected
 * by the fixed parameters from all the component transform in the queue.
 *
 * y = T_2 ( T_1 ( T_0(x-c, p0), p1), p2 ) + c
 *
 * For example, T_0 is a rotation transform, (without translation!)
 *  T_0 (x) = R x, R is a rotation matrix
 * T_1 is a scale transform (without translation!)
 *  T_1 (x) = S x, S is a diagonal matrix
 * T_2 is a translation transform.
 *  T_2 (x) = x + t, t is a vector
 *
 * Thus, the CenteredCompositeTransform of T_0, T_1, T_2 is a centered
 * similarity transform :
 *  y =  S R (x - c) + t + c
 *
 * In practice, we recycle the rotation and scale transform in ITK. These
 * transforms, although have their own fixed parameters, must be reset to
 * all zero, when adding to CenteredCompositeTransform. Otherwise, it will
 * not function properly. Also, the transform should guarantee that
 *  GetParameters() don't include m_Center and thus GetJacobianWithRespectToParameters()
 *
 *
 */
template
<class TScalar = double, unsigned int NDimensions = 3>
class ITK_EXPORT CenteredCompositeTransform :
  public CompositeTransform<TScalar, NDimensions>
{
public:
  /** Standard class typedefs. */
  typedef CenteredCompositeTransform                           Self;
  typedef CompositeTransform<TScalar, NDimensions> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( CenteredCompositeTransform, CompositeTransform );

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

  /** Component transform type **/
  typedef typename Superclass::Superclass                    TransformType;
  typedef typename TransformType::Pointer  TransformTypePointer;
  /** InverseTransform type. */
  typedef typename Superclass::InverseTransformBasePointer
                                        InverseTransformBasePointer;
  /** Scalar type. */
  typedef typename Superclass::ScalarType ScalarType;
  /** Parameters type. */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::ParametersValueType ParametersValueType;
  /** Jacobian type. */
  typedef typename Superclass::JacobianType JacobianType;
  /** Standard coordinate point type for this class. */
  typedef typename Superclass::InputPointType  InputPointType;
  typedef typename Superclass::OutputPointType OutputPointType;
  /** Standard vector type for this class. */
  typedef typename Superclass::InputVectorType  InputVectorType;
  typedef typename Superclass::OutputVectorType OutputVectorType;
  /** Standard covariant vector type for this class */
  typedef typename Superclass::InputCovariantVectorType
  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType
  OutputCovariantVectorType;
  /** Standard vnl_vector type for this class. */
  typedef typename Superclass::InputVnlVectorType  InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType OutputVnlVectorType;
  /** Transform queue type */
  typedef typename Superclass::TransformQueueType TransformQueueType;
  /** Optimization flags queue type */
  typedef typename Superclass::TransformsToOptimizeFlagsType TransformsToOptimizeFlagsType;

  /** Dimension of the domain spaces. */
  itkStaticConstMacro( InputDimension, unsigned int, NDimensions );
  itkStaticConstMacro( OutputDimension, unsigned int, NDimensions );

  /** Functionality for sub transforms */

//  /** Add transforms to the queue, as stack. Only allow one method for simplicity.
//   *  Most-recently added transform is always at back of queue, index N-1.
//   */
//  void AddTransform( TransformType *t  )
//  {
//    this->PushBackTransform( t ); /* Also adds to TransformsToOptimize list */
//  }
//
//  /** See transforms at the front and the back of the queue */
//  const
//  TransformTypePointer GetFrontTransform()
//    {
//    return this->m_TransformQueue.front();
//    }
//  const
//  TransformTypePointer GetBackTransform()
//    {
//    return this->m_TransformQueue.back();
//    }
//
//  const
//  TransformTypePointer GetNthTransform( size_t n ) const
//  {
//    return this->m_TransformQueue[n];
//  }
//
//  /** Active Transform state manipulation */
//
//  void SetNthTransformToOptimize( size_t i, bool state )
//  {
//    this->m_TransformsToOptimizeFlags.at(i)=state;
//    this->Modified();
//  }
//
//  void SetNthTransformToOptimizeOn( size_t i )
//  {
//    this->SetNthTransformToOptimize( i, true );
//  }
//
//  void SetNthTransformToOptimizeOff( size_t i )
//  {
//    this->SetNthTransformToOptimize( i, false );
//  }
//
//  void SetAllTransformsToOptimize( bool state )
//  {
//    this->m_TransformsToOptimizeFlags.assign(
//      this->m_TransformsToOptimizeFlags.size(), state );
//    this->Modified();
//  }
//
//  void SetAllTransformsToOptimizeOn()
//  {
//    this->SetAllTransformsToOptimize( true );
//  }
//
//  void SetAllTransformsToOptimizeOff()
//  {
//    this->SetAllTransformsToOptimize( false );
//  }

//  /* With AddTransform() as the only way to add a transform, we
//   * can have this method to easily allow user to optimize only
//   * the transform added most recenlty. */
//  void SetOnlyMostRecentTransformToOptimizeOn()
//  {
//    this->SetAllTransformsToOptimize( false );
//    this->SetNthTransformToOptimizeOn( this->GetNumberOfTransforms()-1 );
//  }
//
//  /* Get whether the Nth transform is set to be optimzied */
//  /* NOTE: ambiguous function name here - are we getting if the Nth transform
//      is set to be optimized, or the Nth of the transforms that are set to be
//      optimized? */
//  bool GetNthTransformToOptimize( size_t i ) const
//  {
//    return this->m_TransformsToOptimizeFlags.at(i);
//  }
//
//  /** Access transform queue */
//  const TransformQueueType & GetTransformQueue() const
//  { return this->m_TransformQueue; }
//
//
//  /** Access optimize flags */
//  const TransformsToOptimizeFlagsType & GetTransformsToOptimizeFlags() const
//  { return this->m_TransformsToOptimizeFlags; }
//
//  /** Misc. functionality */
//  bool IsTransformQueueEmpty()
//  {
//    return this->m_TransformQueue.empty();
//  }
//
//  size_t GetNumberOfTransforms() const
//  {
//    return this->m_TransformQueue.size();
//  }
//
//  void ClearTransformQueue()
//  {
//    this->m_TransformQueue.clear();
//    this->m_TransformsToOptimizeFlags.clear();
//  }

  /** Return an inverse of this transform. */
  bool GetInverse( Self *inverse ) const;

  virtual InverseTransformBasePointer GetInverseTransform() const;

  /** Compute the position of point in the new space.
  *
  * Transforms are applied starting from the *back* of the
  * queue. That is, in reverse order of which they were added, in order
  * to work properly with ResampleFilter.
  *
  * Imagine a user wants to apply an Affine transform followed by a Deformation
  * Field (DF) transform. He adds the Affine, then the DF. Because the user
  * typically conceptualizes a transformation as being applied from the Moving
  * image to the Fixed image, this makes intuitive sense. But since the
  * ResampleFilter expects to transform from the Fixed image to the Moving
  * image, the transforms are applied in reverse order of addition, i.e. from
  * the back of the queue, and thus, DF then Affine.
  */
  virtual OutputPointType
    TransformPoint( const InputPointType &inputPoint ) const;
  /* Note: why was the 'isInsideTransformRegion' flag used below?
  {
    bool isInside = true;

    return this->TransformPoint( inputPoint, isInside );
  }
  virtual OutputPointType TransformPoint( const InputPointType& thisPoint,
                                          bool &isInsideTransformRegion ) const;
  */
  /**  Method to transform a vector. */
  virtual OutputVectorType TransformVector(const InputVectorType &) const
  {
    itkExceptionMacro( "TransformVector unimplemented" );
  }

  /**  Method to transform a vnl_vector. */
  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &) const
  {
    itkExceptionMacro( "TransformVector unimplemented" );
  }

  /**  Method to transform a CovariantVector. */
  virtual OutputCovariantVectorType
  TransformCovariantVector(const InputCovariantVectorType &) const
  {
    itkExceptionMacro( "TransformCovariantVector unimplemented" );
  }

  virtual bool IsLinear() const;

  /**
   * Compute the jacobian with respect to the parameters.
   */
  virtual const JacobianType & GetJacobian(const InputPointType  &) const;


  virtual void GetJacobianWithRespectToParameters(const InputPointType  &x, JacobianType &j) const;



  /** Get/Set Parameter functions work on the current list of transforms
      that are set to be optimized (active) using the
      'Set[Nth|All]TransformToOptimze' routines.
      The parameter data from each active transform is
      concatenated into a single ParametersType object. */
  virtual const ParametersType & GetParameters(void) const;

  virtual void  SetParameters(const ParametersType & p);

  virtual const ParametersType & GetFixedParameters(void) const;

  virtual void  SetFixedParameters(const ParametersType & fixedParameters);

  virtual unsigned int GetNumberOfParameters(void) const;

  virtual unsigned int GetNumberOfFixedParameters(void) const;

//  virtual void SetCenter(const InputVectorType &center) {
//      m_Center = center;
//  }

  virtual void SetCenter(const InputPointType &center) {
      for(int i=0; i<InputDimension; i++) m_Center[i] = center[i];
  }

  const InputPointType & GetCenter() const {
      InputPointType center;
      for(int i=0; i<InputDimension; i++) center[i] = m_Center[i];
      return m_Center;
  }


protected:
  CenteredCompositeTransform();
  virtual ~CenteredCompositeTransform();
  void PrintSelf( std::ostream& os, Indent indent ) const;

//  void PushFrontTransform( TransformTypePointer t  )
//  {
//    this->m_TransformQueue.push_front( t );
//    /* Add element to list of flags, and set true by default */
//    this->m_TransformsToOptimizeFlags.push_front( true );
//    this->Modified();
//  }
//
//  void PushBackTransform( TransformTypePointer t  )
//  {
//    this->m_TransformQueue.push_back( t );
//    /* Add element to list of flags, and set true by default */
//    this->m_TransformsToOptimizeFlags.push_back( true );
//    this->Modified();
//  }
//
//  /** Transform container object. */
//  mutable TransformQueueType m_TransformQueue;
//
//  /** Get a list of transforms to optimize. Helper function. */
//  TransformQueueType& GetTransformsToOptimizeQueue() const;
//
//  mutable TransformQueueType            m_TransformsToOptimizeQueue;
//  mutable TransformsToOptimizeFlagsType m_TransformsToOptimizeFlags;
private:
  CenteredCompositeTransform( const Self & ); //purposely not implemented
  void operator=( const Self& );      //purposely not implemented

  // not sure if InputPointType suitable for Point+Point
  InputVectorType   m_Center;

  mutable unsigned long m_PreviousTransformsToOptimizeUpdateTime;
};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkCenteredCompositeTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkCenteredCompositeTransform.hxx"
#endif

#endif // __itkCompositeTransform_h
