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
#ifndef __itkNewClassExample_txx
#define __itkNewClassExample_txx

#include "itkNewClassExample.h"

namespace itk
{
template< class TFixedImage, class TMovingImage >
typename NewClassExample< TFixedImage, TMovingImage >
::MeasureType
NewClassExample< TFixedImage, TMovingImage >
::EvaluateMeasure(HistogramType & histogram) const
{
  MeasureType            measure = NumericTraits< MeasureType >::Zero;
  HistogramIteratorType  it = histogram.Begin();
  HistogramIteratorType  end = histogram.End();
  HistogramFrequencyType totalNumberOfSamples =
    NumericTraits< HistogramFrequencyType >::Zero;

  while ( it != end )
    {
    HistogramFrequencyType freq = it.GetFrequency();
    if ( freq > 0 )
      {
      HistogramMeasurementVectorType value = it.GetMeasurementVector();
      measure += ( value[0] - value[1] ) * ( value[0] - value[1] ) * freq;
      totalNumberOfSamples += freq;
      }
    ++it;
    }

  measure /= totalNumberOfSamples;

  return measure;
}

template< class TFixedImage, class TMovingImage >
template< typename TType >
void
NewClassExample< TFixedImage, TMovingImage >
::TestTemplatedMethod( TType input )
{
  TType val = input + NumericTraits< TType >::Zero;
  std::cout << val << std::endl;
}
} // End namespace itk

#endif // itkNewClassExample_txx
