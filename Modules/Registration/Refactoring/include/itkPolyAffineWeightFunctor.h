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
#ifndef __itkPolyAffineWeightFunctor_h
#define __itkPolyAffineWeightFunctor_h

#include "itkImageBoundaryCondition.h"
#include "itkImageBase.h"

namespace itk
{
/** \class PolyAffineWeightFunctor
 * \brief Provides accessor interfaces to Get pixels and is meant to be
 * used on pointers contained within Neighborhoods. A typical user should
 * not need to use this class directly. This class is used by the
 * neighborhood iterators to get pixels from pixel pointers or assign
 * a pixel to an address.
 *
 *
 * \note
 * This work is part of the National Alliance for Medical Image Computing
 * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
 * for Medical Research, Grant U54 EB005149.
  * \ingroup ITKCommon
 */
template< class TScalarType, unsigned int NSpaceDimensions >
class PolyAffineWeightFunctor : public Object
{
public:
  typedef PolyAffineWeightFunctor    Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(PolyAffineWeightFunctor, Object);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Dimension of the domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NSpaceDimensions);

  typedef TScalarType ScalarType;

  /** Standard coordinate point type for this class   */
  typedef Point< ScalarType,
                 itkGetStaticConstMacro(SpaceDimension) >
  PointType;

  typedef typename PointType::ValueType PointValueType;

  inline void SetAnchor(const PointType anchor)
  {
    m_Anchor = anchor;
  }

  inline ScalarType GetWeight(const PointType point) const
  {
    return 1;

    double squaredDistance = m_Anchor.SquaredEuclideanDistanceTo( point );
    return 1.0 / (1 + squaredDistance);

    //double sigma = 100;
    //return vcl_exp( - squaredDistance / sigma / sigma );
  }

  /** Constructor */
  PolyAffineWeightFunctor() { m_Anchor.Fill(0.0); }

  /** Destructor */
  virtual ~PolyAffineWeightFunctor() {}

private:

  PointType m_Anchor;

};
} // end namespace itk

//template< class TImage > const unsigned int
// itk::PolyAffineWeightFunctor<TImage>::ImageDimension;

#endif
