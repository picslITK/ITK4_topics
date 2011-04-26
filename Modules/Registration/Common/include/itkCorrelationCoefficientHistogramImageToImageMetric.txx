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
#ifndef __itkCorrelationCoefficientHistogramImageToImageMetric_txx
#define __itkCorrelationCoefficientHistogramImageToImageMetric_txx

#include "itkCorrelationCoefficientHistogramImageToImageMetric.h"

namespace itk
{
template< class TFixedImage, class TMovingImage, typename TValueType >
typename CorrelationCoefficientHistogramImageToImageMetric< TFixedImage,
                                                            TMovingImage,
                                                            TValueType >
                                                            ::MeasureType
CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >
::EvaluateMeasure(HistogramType & histogram) const
{
  const MeasureType varianceX  = this->VarianceX(histogram);
  const MeasureType varianceY  = this->VarianceY(histogram);
  const MeasureType covariance = this->Covariance(histogram);

  return vcl_fabs( covariance / ( vcl_sqrt(varianceX) * vcl_sqrt(varianceY) ) );
}

template< class TFixedImage, class TMovingImage, typename TValueType >
typename CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, \
                                                            TMovingImage,
                                                            TValueType >
                                                            ::MeasureType
CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >
::MeanX(HistogramType & histogram) const
{
  MeasureType meanX = NumericTraits< MeasureType >::Zero;

  for ( unsigned int i = 0; i < this->m_HistogramSize[0]; i++ )
    {
    MeasureType            valX = histogram.GetMeasurement(i, 0);
    HistogramFrequencyType freq = histogram.GetFrequency(i, 0);
    meanX += valX * freq;
    }

  meanX /= histogram.GetTotalFrequency();

  return meanX;
}

template< class TFixedImage, class TMovingImage, typename TValueType >
typename CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, \
                                                            TMovingImage,
                                                            TValueType >
                                                            ::MeasureType
CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >
::MeanY(HistogramType & histogram) const
{
  MeasureType meanY = NumericTraits< MeasureType >::Zero;

  for ( unsigned int i = 0; i < this->m_HistogramSize[1]; i++ )
    {
    MeasureType            valY = histogram.GetMeasurement(i, 1);
    HistogramFrequencyType freq = histogram.GetFrequency(i, 1);
    meanY += valY * freq;
    }

  meanY /= histogram.GetTotalFrequency();

  return meanY;
}

template< class TFixedImage, class TMovingImage, typename TValueType >
typename CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, \
                                                            TMovingImage,
                                                            TValueType >
                                                            ::MeasureType
CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >
::VarianceX(HistogramType & histogram) const
{
  MeasureType varX = NumericTraits< MeasureType >::Zero;

  for ( unsigned int i = 0; i < this->m_HistogramSize[0]; i++ )
    {
    varX += static_cast< double >( histogram.GetFrequency(i, 0) )
            / histogram.GetTotalFrequency()
            * vcl_pow(histogram.GetMeasurement(i, 0), 2);
    }

  return varX - vcl_pow(MeanX(histogram), 2);
}

template< class TFixedImage, class TMovingImage, typename TValueType >
typename CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, \
                                                            TMovingImage,
                                                            TValueType >
                                                            ::MeasureType
CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >
::VarianceY(HistogramType & histogram) const
{
  MeasureType varY = NumericTraits< MeasureType >::Zero;

  for ( unsigned int i = 0; i < this->m_HistogramSize[1]; i++ )
    {
    varY += static_cast< double >( histogram.GetFrequency(i, 1) )
            / histogram.GetTotalFrequency()
            * vcl_pow(histogram.GetMeasurement(i, 1), 2);
    }

  return varY - vcl_pow(MeanY(histogram), 2);
}

template< class TFixedImage, class TMovingImage, typename TValueType >
typename CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, \
                                                            TMovingImage,
                                                            TValueType >
                                                            ::MeasureType
CorrelationCoefficientHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >
::Covariance(HistogramType & histogram) const
{
  MeasureType var = NumericTraits< MeasureType >::Zero;
  MeasureType meanX = MeanX(histogram);
  MeasureType meanY = MeanY(histogram);

  for ( unsigned int j = 0; j < this->m_HistogramSize[1]; j++ )
    {
    for ( unsigned int i = 0; i < this->m_HistogramSize[0]; i++ )
      {
      typename HistogramType::IndexType index;
      index.SetSize(2);
      index[0] = i;
      index[1] = j;

      var += histogram.GetFrequency(index)
             * ( histogram.GetMeasurement(i, 0) - meanX )
             * ( histogram.GetMeasurement(j, 1) - meanY );
      }
    }

  var /= histogram.GetTotalFrequency();

  return var;
}
}

#endif // itkCorrelationCoefficientHistogramImageToImageMetric_txx
