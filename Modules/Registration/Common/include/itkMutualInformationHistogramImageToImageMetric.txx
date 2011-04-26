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
#ifndef __itkMutualInformationHistogramImageToImageMetric_txx
#define __itkMutualInformationHistogramImageToImageMetric_txx

#include "itkMutualInformationHistogramImageToImageMetric.h"
#include "itkHistogram.h"

namespace itk
{
template< class TFixedImage, class TMovingImage, typename TValueType >
typename MutualInformationHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >::MeasureType
MutualInformationHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >
::EvaluateMeasure(HistogramType & histogram) const
{
  MeasureType entropyX = NumericTraits< MeasureType >::Zero;
  MeasureType entropyY = NumericTraits< MeasureType >::Zero;
  MeasureType jointEntropy = NumericTraits< MeasureType >::Zero;

  typedef typename NumericTraits< HistogramFrequencyType >::RealType HistogramFrequencyRealType;

  HistogramFrequencyRealType totalFreq =
    static_cast< HistogramFrequencyRealType >( histogram.GetTotalFrequency() );

  for ( unsigned int i = 0; i < this->GetHistogramSize()[0]; i++ )
    {
    HistogramFrequencyRealType freq =
      static_cast< HistogramFrequencyRealType >( histogram.GetFrequency(i, 0) );
    if ( freq > 0 )
      {
      entropyX += freq * vcl_log(freq);
      }
    }

  entropyX = -entropyX / static_cast< MeasureType >( totalFreq ) + vcl_log(totalFreq);

  for ( unsigned int i = 0; i < this->GetHistogramSize()[1]; i++ )
    {
    HistogramFrequencyRealType freq =
      static_cast< HistogramFrequencyRealType >( histogram.GetFrequency(i, 1) );
    if ( freq > 0 )
      {
      entropyY += freq * vcl_log(freq);
      }
    }

  entropyY = -entropyY / static_cast< MeasureType >( totalFreq ) + vcl_log(totalFreq);

  HistogramIteratorType it = histogram.Begin();
  HistogramIteratorType end = histogram.End();
  while ( it != end )
    {
    HistogramFrequencyRealType freq =
      static_cast< HistogramFrequencyRealType >( it.GetFrequency() );
    if ( freq > 0 )
      {
      jointEntropy += freq * vcl_log(freq);
      }
    ++it;
    }

  jointEntropy = -jointEntropy
                 / static_cast< MeasureType >( totalFreq ) + vcl_log(totalFreq);

  return entropyX + entropyY - jointEntropy;
}
} // End namespace itk

#endif // itkMutualInformationHistogramImageToImageMetric_txx
