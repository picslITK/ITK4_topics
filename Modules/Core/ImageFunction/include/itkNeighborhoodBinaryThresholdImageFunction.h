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
#ifndef __itkNeighborhoodBinaryThresholdImageFunction_h
#define __itkNeighborhoodBinaryThresholdImageFunction_h

#include "itkBinaryThresholdImageFunction.h"

namespace itk
{
/**
 * \class NeighborhoodBinaryThresholdImageFunction
 * \brief Determine whether all the pixels in the specified neighborhood meet a threshold criteria
 *
 * Determine whether all the pixels in the specified neighborhood meet
 * a threshold criteria.
 *
 * If called with a ContinuousIndex or Point, the calculation is performed
 * at the nearest neighbor.
 *
 * This class is templated over the input image type and the coordinate
 * representation type (e.g. float or double).
 *
 * \ingroup ImageFunctions
 * \ingroup ITK-ImageFunction
 */
template< class TInputImage, class TCoordRep = float >
class ITK_EXPORT NeighborhoodBinaryThresholdImageFunction:
  public BinaryThresholdImageFunction< TInputImage, TCoordRep >
{
public:
  /** Standard class typedefs. */
  typedef NeighborhoodBinaryThresholdImageFunction               Self;
  typedef BinaryThresholdImageFunction< TInputImage, TCoordRep > Superclass;
  typedef SmartPointer< Self >                                   Pointer;
  typedef SmartPointer< const Self >                             ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(NeighborhoodBinaryThresholdImageFunction, BinaryThresholdImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** InputImageType typedef support. */
  typedef TInputImage InputImageType;

  /** OutputType typdef support. */
  typedef typename Superclass::OutputType OutputType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** PixelType typedef support. */
  typedef typename Superclass::PixelType PixelType;

  /** Dimension of the underlying image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      InputImageType::ImageDimension);

  /** SizeType of the input image */
  typedef typename InputImageType::SizeType InputSizeType;

  /** Set the radius of the neighborhood used in computation. */
  itkSetMacro(Radius, InputSizeType);

  /** Get the radius of the neighborhood used in computation */
  itkGetConstReferenceMacro(Radius, InputSizeType);

  /** Evalulate the function at specified index */
  virtual bool EvaluateAtIndex(const IndexType & index) const;

  /** Evaluate the function at non-integer positions */
  virtual bool Evaluate(const PointType & point) const
  {
    IndexType index;

    this->ConvertPointToNearestIndex(point, index);
    return this->EvaluateAtIndex(index);
  }

  virtual bool EvaluateAtContinuousIndex(
    const ContinuousIndexType & cindex) const
  {
    IndexType index;

    this->ConvertContinuousIndexToNearestIndex(cindex, index);
    return this->EvaluateAtIndex(index);
  }

protected:
  NeighborhoodBinaryThresholdImageFunction();
  ~NeighborhoodBinaryThresholdImageFunction(){}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  NeighborhoodBinaryThresholdImageFunction(const Self &); //purposely not
                                                          // implemented
  void operator=(const Self &);                           //purposely not

  // implemented

  InputSizeType m_Radius;
};
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_NeighborhoodBinaryThresholdImageFunction(_, EXPORT, TypeX, TypeY)     \
  namespace itk                                                                            \
  {                                                                                        \
  _( 2 ( class EXPORT NeighborhoodBinaryThresholdImageFunction< ITK_TEMPLATE_2 TypeX > ) ) \
  namespace Templates                                                                      \
  {                                                                                        \
  typedef NeighborhoodBinaryThresholdImageFunction< ITK_TEMPLATE_2 TypeX >                 \
  NeighborhoodBinaryThresholdImageFunction##TypeY;                                       \
  }                                                                                        \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkNeighborhoodBinaryThresholdImageFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkNeighborhoodBinaryThresholdImageFunction.txx"
#endif

/*
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNeighborhoodBinaryThresholdImageFunction.txx"
#endif
*/

#endif
