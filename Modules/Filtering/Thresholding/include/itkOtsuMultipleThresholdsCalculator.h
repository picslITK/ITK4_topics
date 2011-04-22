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
#ifndef __itkOtsuMultipleThresholdsCalculator_h
#define __itkOtsuMultipleThresholdsCalculator_h

#include "itkHistogramAlgorithmBase.h"
#include "itkHistogram.h"

namespace itk
{
/** \class OtsuMultipleThresholdsCalculator
 * \brief Computes Otsu's thresholds for a histogram.
 *
 * You plug in the target histogram using SetInputHistogram method and
 * specify the number of thresholds you want to be computed. Then call
 * the GenerateData method to run the alogithm.
 *
 * The thresholds are computed so that the between-class variance is
 * maximized.
 *
 * \ingroup Calculators
 * \ingroup ITK-Thresholding
 */

template< class TInputHistogram >
class ITK_EXPORT OtsuMultipleThresholdsCalculator:
  public HistogramAlgorithmBase< TInputHistogram >
{
public:
  /**Standard class typedefs. */
  typedef OtsuMultipleThresholdsCalculator          Self;
  typedef HistogramAlgorithmBase< TInputHistogram > Superclass;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;

  typedef typename TInputHistogram::MeasurementType       MeasurementType;
  typedef typename TInputHistogram::AbsoluteFrequencyType FrequencyType;

  typedef typename NumericTraits< MeasurementType >::RealType MeanType;
  typedef typename NumericTraits< MeasurementType >::RealType VarianceType;

  typedef std::vector< MeanType >      MeanVectorType;
  typedef std::vector< FrequencyType > FrequencyVectorType;

  typedef typename TInputHistogram::InstanceIdentifier InstanceIdentifierType;
  typedef std::vector< InstanceIdentifierType >        InstanceIdentifierVectorType;

  /**Standard Macros */
  itkTypeMacro(OtsuMultipleThresholdsCalculator, HistogramAlgorithmsBase);
  itkNewMacro(Self);

  /** Typedef for the thresholds output */
  typedef std::vector< MeasurementType > OutputType;

  /** Returns the thresholds vector */
  const OutputType & GetOutput();

  /** Set/Get the number of thresholds. */
  itkSetClampMacro( NumberOfThresholds, SizeValueType, 1, NumericTraits< SizeValueType >::max() );
  itkGetConstMacro(NumberOfThresholds, SizeValueType);
protected:
  OtsuMultipleThresholdsCalculator();
  virtual ~OtsuMultipleThresholdsCalculator() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Calculates the thresholds and save them */
  void GenerateData();

  /** Increment the thresholds of one position */
  bool IncrementThresholds(InstanceIdentifierVectorType & thresholdIds,
                           MeanType totalMean,
                           MeanVectorType & classMean,
                           FrequencyVectorType & classFrequency);

private:
  /** Internal thresholds storage */
  SizeValueType m_NumberOfThresholds;
  OutputType    m_Output;
}; // end of class
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOtsuMultipleThresholdsCalculator.txx"
#endif

#endif
