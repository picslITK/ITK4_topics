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
#ifndef __itkUpdateTransformFunction_h
#define __itkUpdateTransformFunction_h


namespace itk
{

/** \class UpdateTransformFunction
 * \brief Provides functionality for a transform to update its
 * parameters during a call to \c Transform::UpdateTransformParameters.
 *
 * This class is assigned to a transform to perform a particular
 * operation during a call to Transform::UpdateTransformParameters. Typically,
 * \c Transform::UpdateTransformParameters is called during optimization
 * by the optimizer with an array of parameter gradient deltas to apply.
 *
 * This base class provides single-threaded operation, simply adding
 * the update array to the transform's parameters after applying an optional
 * scaling factor.
 *
 * Derived classes should override the \c Update method to provide
 * particular functionality. The pointer of the transform holding this
 * function is passed to \c Update to provide flexibility, at the cost
 * of elegance.
 *
 * This functionality is implemented as a function class for two reasons:
 * 1) allow changing parameter update functionality of a transform at runtime,
 * to avoid having to create a new transform of a different type and copy or
 * swap a parameter state.
 * 2) provide different update functionality to a range of transform types
 * without having to create a new derived transform class for each type of
 * functionality and each transform.
 *
 * \ingroup Transforms
 */
template
  <class TTransform >
class ITK_EXPORT UpdateTransformFunction : public Object
{
public:
  /** Standard class typedefs. */
  typedef UpdateTransformFunction                           Self;
  typedef Object                                            Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( UpdateTransformFunction, Object );

  /** New macro for creation through a Smart Pointer */
  itkNewMacro( Self );

  typedef TTransform                                TransformType;
  typedef typename TransformType::Pointer           TransformPointer;
  typedef typename TransformType::DerivativeType    DerivativeType;
  typedef typename TransformType::ScalarType        ScalarType;

  /** Update method. Derived classes should override this to provide
   * new functionality. */
  virtual void Update( DerivativeType & update,
                       ScalarType factor,
                       TransformType * transform );
protected:

  UpdateTransformFunction(){}
  /** Destroy an AffineTransform object   */
  virtual ~UpdateTransformFunction() {}

private:

  UpdateTransformFunction(const Self & other);
  const Self & operator=(const Self &);

};
} //namespace itk

#if ITK_TEMPLATE_TXX
#include "itkUpdateTransformFunction.txx"
#endif

#endif
