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

#ifndef __itkDiscreteLevelSetImageBase_h
#define __itkDiscreteLevelSetImageBase_h

#include "itkLevelSetImageBase.h"
#include "itkObjectFactory.h"
#include "itkIndex.h"
#include "itkImageBase.h"

namespace itk
{
/**
 *  \class DiscreteLevelSetImageBase
 *  \brief Abstract class for a level-set function on one Image.
 *
 *  \tparam TOutput OutputType of the level-set function value
 *  \tparam VDimension Dimension of the underlying Image.
 *
 *  \ingroup ITKLevelSetsv4
 */
template< typename TOutput, unsigned int VDimension >
class DiscreteLevelSetImageBase :
  public LevelSetImageBase< Index< VDimension >,
                       VDimension,
                       TOutput >
{
public:
  typedef Index< VDimension >             IndexType;

  typedef DiscreteLevelSetImageBase                           Self;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;
  typedef LevelSetImageBase< IndexType, VDimension, TOutput > Superclass;

  /** Run-time type information */
  itkTypeMacro ( DiscreteLevelSetImageBase, LevelSetImageBase );

  itkStaticConstMacro ( Dimension, unsigned int, Superclass::Dimension );

  typedef typename Superclass::InputType        InputType;
  typedef typename Superclass::OutputType       OutputType;
  typedef typename Superclass::OutputRealType   OutputRealType;
  typedef typename Superclass::GradientType     GradientType;
  typedef typename Superclass::HessianType      HessianType;
  typedef typename Superclass::LevelSetDataType LevelSetDataType;

  /** Returns the gradient of the level set function at a given location iP */
  virtual OutputType  Evaluate( const InputType& iP ) const = 0;

  /** Returns the image gradient of the level set function at a given location iP */
  virtual GradientType EvaluateGradient( const InputType& iP ) const;

  /** Returns the image hessian of the level set function at a given location iP */
  virtual HessianType EvaluateHessian( const InputType& iP ) const;

  /** Returns the image Laplacian of the level set function at a given location iP */
  virtual OutputRealType EvaluateLaplacian( const InputType& iP ) const;

  /** Returns the mean curvature of the level set function at a given location iP */
  virtual OutputRealType EvaluateMeanCurvature( const InputType& iP ) const;

  virtual GradientType EvaluateForwardGradient( const InputType& iP ) const;

  virtual GradientType EvaluateBackwardGradient( const InputType& iP ) const;

  /** Returns the value of the level set function at a given location iP */
  virtual void Evaluate( const InputType& iP, LevelSetDataType& ioData ) const;

  /** Returns the gradient of the level set function at a given location iP
   * as part of the LevelSetDataType */
  virtual void EvaluateGradient( const InputType& iP, LevelSetDataType& ioData ) const;

  /** Returns the Hessian of the level set function at a given location iP
   * as part of the LevelSetDataType */
  virtual void EvaluateHessian( const InputType& iP, LevelSetDataType& ioData ) const;

  /** Returns the Hessian of the level set function at a given location iP
   * as part of the LevelSetDataType */
  virtual void EvaluateMeanCurvature( const InputType& iP, LevelSetDataType& ioData ) const;

  /** Returns the Laplacian of the level set function at a given location iP
   * as part of the LevelSetDataType */
  virtual void EvaluateLaplacian( const InputType& iP, LevelSetDataType& ioData ) const;

  /** Returns the gradient of the level set function at a given location iP
   * as part of the LevelSetDataType */
  virtual void EvaluateForwardGradient( const InputType& iP, LevelSetDataType& ioData ) const;

  /** Returns the gradient of the level set function at a given location iP
   * as part of the LevelSetDataType */
  virtual void EvaluateBackwardGradient( const InputType& iP, LevelSetDataType& ioData ) const;

protected:
  DiscreteLevelSetImageBase();

  virtual ~DiscreteLevelSetImageBase();

  /** Initial the level set pointer */
  virtual void Initialize();

  /** Copy level set information from data object */
  virtual void CopyInformation(const DataObject *data);

  /** Graft data object as level set object */
  virtual void Graft( const DataObject* data );

private:

  DiscreteLevelSetImageBase( const Self& ); // purposely not implemented
  void operator = ( const Self& ); // purposely not implemented
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiscreteLevelSetImageBase.hxx"
#endif

#endif // __itkDiscreteLevelSetImageBase_h
