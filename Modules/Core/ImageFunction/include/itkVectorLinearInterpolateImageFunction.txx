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
#ifndef __itkVectorLinearInterpolateImageFunction_txx
#define __itkVectorLinearInterpolateImageFunction_txx

#include "itkVectorLinearInterpolateImageFunction.h"

#include "vnl/vnl_math.h"

namespace itk
{
/**
 * Define the number of neighbors
 */
template< class TInputImage, class TCoordRep >
const unsigned long
VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;

/**
 * Constructor
 */
template< class TInputImage, class TCoordRep >
VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::VectorLinearInterpolateImageFunction()
{}

/**
 * PrintSelf
 */
template< class TInputImage, class TCoordRep >
void
VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
}

/**
 * Evaluate at image index position
 */
template< class TInputImage, class TCoordRep >
typename VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::OutputType
VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(
  const ContinuousIndexType & index) const
{
  unsigned int dim;  // index over dimension

  /**
   * Compute base index = closet index below point
   * Compute distance from point to base index
   */
  IndexType baseIndex;
  double    distance[ImageDimension];

  for ( dim = 0; dim < ImageDimension; dim++ )
    {
    baseIndex[dim] = Math::Floor< IndexValueType >(index[dim]);
    distance[dim] = index[dim] - static_cast< double >( baseIndex[dim] );
    }

  /**
   * Interpolated value is the weighted sum of each of the surrounding
   * neighbors. The weight for each neighbor is the fraction overlap
   * of the neighbor pixel with respect to a pixel centered on point.
   */
  OutputType output;
  output.Fill(0.0);

  typedef typename NumericTraits< PixelType >::ScalarRealType ScalarRealType;
  ScalarRealType totalOverlap = NumericTraits< ScalarRealType >::Zero;

  for ( unsigned int counter = 0; counter < m_Neighbors; counter++ )
    {
    double       overlap = 1.0;    // fraction overlap
    unsigned int upper = counter;  // each bit indicates upper/lower neighbour
    IndexType    neighIndex;

    // get neighbor index and overlap fraction
    for ( dim = 0; dim < ImageDimension; dim++ )
      {
      if ( upper & 1 )
        {
        neighIndex[dim] = baseIndex[dim] + 1;
        // Take care of the case where the pixel is just
        // in the outer upper boundary of the image grid.
        if ( neighIndex[dim] > this->m_EndIndex[dim] )
          {
          neighIndex[dim] = this->m_EndIndex[dim];
          }
        overlap *= distance[dim];
        }
      else
        {
        neighIndex[dim] = baseIndex[dim];
        // Take care of the case where the pixel is just
        // in the outer lower boundary of the image grid.
        if ( neighIndex[dim] < this->m_StartIndex[dim] )
          {
          neighIndex[dim] = this->m_StartIndex[dim];
          }
        overlap *= 1.0 - distance[dim];
        }

      upper >>= 1;
      }

    // get neighbor value only if overlap is not zero
    if ( overlap )
      {
      const PixelType input = this->GetInputImage()->GetPixel(neighIndex);
      for ( unsigned int k = 0; k < PixelType::GetNumberOfComponents(); k++ )
        {
        output[k] += overlap * static_cast< RealType >( input[k] );
        }
      totalOverlap += overlap;
      }

    if ( totalOverlap == 1.0 )
      {
      // finished
      break;
      }
    }

  return ( output );
}
} // end namespace itk

#endif
