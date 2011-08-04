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
#ifndef __itkVectorCentralDifferenceImageFunction_hxx
#define __itkVectorCentralDifferenceImageFunction_hxx

#include "itkVectorCentralDifferenceImageFunction.h"

namespace itk
{
/**
 * Constructor
 */
template< class TInputImage, class TCoordRep >
VectorCentralDifferenceImageFunction< TInputImage, TCoordRep >
::VectorCentralDifferenceImageFunction()
{
  this->m_UseImageDirection = true;
}

/**
 *
 */
template< class TInputImage, class TCoordRep >
void
VectorCentralDifferenceImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "UseImageDirection = " << this->m_UseImageDirection << std::endl;
}

/**
 *
 */
template< class TInputImage, class TCoordRep >
typename VectorCentralDifferenceImageFunction< TInputImage, TCoordRep >::OutputType
VectorCentralDifferenceImageFunction< TInputImage, TCoordRep >
::EvaluateAtIndex(const IndexType & index) const
{

  unsigned int numberOfComponents
    = NumericTraits<PixelType>::GetLength(this->GetInputImage()->GetPixel( index ));

  OutputType derivative;
  derivative.SetSize(numberOfComponents, ImageDimension );
  derivative.Fill(0.0);

  IndexType neighIndex = index;

  const InputImageType *inputImage = this->GetInputImage();

  const typename InputImageType::RegionType region =
    inputImage->GetBufferedRegion();

  const typename InputImageType::SizeType & size   = region.GetSize();
  const typename InputImageType::IndexType & start = region.GetIndex();

  const unsigned int MaxDims = Self::ImageDimension;
  for ( unsigned int dim = 0; dim < MaxDims; dim++ )
    {
    // bounds checking
    if ( index[dim] < start[dim] + 1
         || index[dim] > ( start[dim] + static_cast< OffsetValueType >( size[dim] ) - 2 ) )
      {
      for (unsigned int i=0; i<numberOfComponents; i++)
        {
        derivative(i,dim) = 0.0;
        }
      continue;
      }

    // compute derivative
    for (unsigned int i=0; i<numberOfComponents; i++)
      {
      neighIndex[dim] += 1;
      derivative(i,dim) = inputImage->GetPixel(neighIndex)[i];

      neighIndex[dim] -= 2;
      derivative(i,dim) -= inputImage->GetPixel(neighIndex)[i];

      derivative(i,dim) *= 0.5 / inputImage->GetSpacing()[dim];
      neighIndex[dim] += 1;
      }
    }

  if ( this->m_UseImageDirection )
    {
    //OutputType orientedDerivative;
    //inputImage->TransformLocalVectorToPhysicalVector(derivative, orientedDerivative);
    //return orientedDerivative;
    }

  return ( derivative );
}
} // end namespace itk

#endif
