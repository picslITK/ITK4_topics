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
#ifndef __itkSmoothedGaussianDeformationFieldUpdateTransformFunction_h
#define __itkSmoothedGaussianDeformationFieldUpdateTransformFunction_h


namespace itk
{

/** \class SmoothedGaussianDeformationFieldUpdateTransformFunction
 * \brief Provides functionality for a deformation field transform to update
 * its parameters during a call to \c Transform::UpdateTransformParameters,
 * using multi-threaded addition (with optional scaling) of the update,
 * followed by gaussian smoothing of the result.
 *
 * This class is used by default in DeformationFieldTransform. To use
 * a different UpdateTransformFunction, call SetUpdateTransformFunction.
 *
 * \ingroup Transforms
 */
template
  <class TDeformationFieldTransform >
class ITK_EXPORT SmoothedGaussianDeformationFieldUpdateTransformFunction
  : public UpdateTransformFunction< TDeformationFieldTransform >
{
public:
  /** Standard class typedefs. */
  typedef SmoothedGaussianDeformationFieldUpdateTransformFunction Self;
  typedef UpdateTransformFunction< TDeformationFieldTransform >   Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(
          SmoothedGaussianDeformationFieldUpdateTransformFunction,
          UpdateTransformFunction );

  /** New macro for creation through a Smart Pointer */
  itkNewMacro( Self );

  typedef TDeformationFieldTransform                TransformType;
  typedef typename TransformType::Pointer           TransformPointer;
  typedef typename TransformType::DerivativeType    DerivativeType;
  typedef typename TransformType::ScalarType        ScalarType;

  /** Update method. Derived classes should override this to provide
   * new functionality. */
  virtual void Update( DerivativeType & update,
                       ScalarType factor,
                       TransformType * transform );
protected:

  SmoothedGaussianDeformationFieldUpdateTransformFunction(){}
  /** Destroy an AffineTransform object   */
  virtual ~SmoothedGaussianDeformationFieldUpdateTransformFunction() {}

private:

  SmoothedGaussianDeformationFieldUpdateTransformFunction(const Self & other);
  const Self & operator=(const Self &);

};
} //namespace itk

#if ITK_TEMPLATE_TXX
#include "itkSmoothedGaussianDeformationFieldUpdateTransformFunction.txx"
#endif

#endif