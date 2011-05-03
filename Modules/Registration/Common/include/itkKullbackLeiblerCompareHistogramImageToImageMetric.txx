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
#ifndef __itkKullbackLeiblerCompareHistogramImageToImageMetric_txx
#define __itkKullbackLeiblerCompareHistogramImageToImageMetric_txx

#include "itkKullbackLeiblerCompareHistogramImageToImageMetric.h"
#include "itkHistogram.h"

// Todo: need to access Use_Padding in parent. Make in protected
// need to figure out what to do when "stuff" is not in histogram
// kernel function?

namespace itk
{
template< class TFixedImage, class TMovingImage, typename TValueType >
KullbackLeiblerCompareHistogramImageToImageMetric< TFixedImage,
                                                   TMovingImage,
                                                   TValueType>
::KullbackLeiblerCompareHistogramImageToImageMetric()
{
  m_Epsilon                = 1e-12; // should be smaller than 1/numBins^2
}

template< class TFixedImage, class TMovingImage, typename TValueType >
void
KullbackLeiblerCompareHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >
::Initialize()
throw ( ExceptionObject )
{
  Superclass::Initialize();
}

template< class TFixedImage, class TMovingImage, typename TValueType >
typename KullbackLeiblerCompareHistogramImageToImageMetric< TFixedImage, \
                                                            TMovingImage,
                                                            TValueType>
::MeasureType
KullbackLeiblerCompareHistogramImageToImageMetric< TFixedImage, \
                                                   TMovingImage,
                                                   TValueType >
::EvaluateMeasure(HistogramType & histogram) const
{
  // Two terms.
  // First the term that measures the entropy of the term
  // p(x,y) log p(x,y) - p(x,y) log q(x,y)

  MeasureType KullbackLeibler = NumericTraits< MeasureType >::Zero;

  HistogramIteratorType measured_it   = histogram.Begin();
  HistogramIteratorType measured_end  = histogram.End();

  HistogramIteratorType training_it   = this->GetTrainingHistogram()->Begin();
  HistogramIteratorType training_end  = this->GetTrainingHistogram()->End();

  while ( measured_it != measured_end )
    {
    // Every bin gets epsilon added to it
    double TrainingFreq = training_it.GetFrequency() + m_Epsilon;
    double MeasuredFreq = measured_it.GetFrequency() + m_Epsilon;

    KullbackLeibler += MeasuredFreq * vcl_log(MeasuredFreq / TrainingFreq);

    ++measured_it;
    ++training_it;
    }

  if ( training_it != training_end )
    {
    itkWarningMacro("The Measured and Training Histograms have different number of bins.");
    }

  // Get the total frequency for each histogram.
  HistogramFrequencyType totalTrainingFreq = this->GetTrainingHistogram()->GetTotalFrequency();
  HistogramFrequencyType totalMeasuredFreq = histogram.GetTotalFrequency();

  // The actual number of total frequency is a bit larger
  // than the number of counts because we add m_Epsilon to every bin
  double AdjustedTotalTrainingFreq = totalTrainingFreq
                                     + this->GetHistogramSize()[0] * this->GetHistogramSize()[1] * m_Epsilon;
  double AdjustedTotalMeasuredFreq = totalMeasuredFreq
                                     + this->GetHistogramSize()[0] * this->GetHistogramSize()[1] * m_Epsilon;

  KullbackLeibler = KullbackLeibler / static_cast< MeasureType >( AdjustedTotalMeasuredFreq )
                    - vcl_log(AdjustedTotalMeasuredFreq / AdjustedTotalTrainingFreq);

  return KullbackLeibler;
}

template< class TFixedImage, class TMovingImage, typename TValueType >
void KullbackLeiblerCompareHistogramImageToImageMetric< TFixedImage, TMovingImage, TValueType >::PrintSelf(std::ostream & os,
                                                                                               Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Epsilon: " << m_Epsilon << std::endl;
}
} // End namespace itk

#endif // itkKullbackLeiblerCompareHistogramImageToImageMetric_txx
