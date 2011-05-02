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
#ifndef __itkHistogramToEntropyImageFilter_h
#define __itkHistogramToEntropyImageFilter_h

#include "itkHistogramToImageFilter.h"

namespace itk
{
/** \class HistogramToEntropyImageFilter
 * \brief The class takes a histogram as an input and gives the entropy
 * image as the output. A pixel, at position I,  in the output image is given by
 *
 * \f[
 * f(I) = -p \log_2 p
 * \f]
 *
 * where
 * \f[
 * p = \frac{q_I}{\sum_{i \in I} q_I}
 * \f]
 *  where  \f$q_I\f$ is the frequency of measurement vector, I.
 *
 * \f$p\f$ is the frequency of a measurement vector by the sum of all frequencies =
 * Probability of the the measurement vector
 *
 * The output image is of type double.
 *
 * This is useful in plotting the joint histograms during registration.
 *
 *  \sa HistogramToImageFilter, HistogramToLogProbabilityImageFilter,
 *  HistogramToIntensityImageFilter, HistogramToProbabilityImageFilter
 *
 * \ingroup ITK-Statistics
 */

namespace Function
{
template< class TInput, class TOutput = double >
class HistogramEntropyFunction
{
public:

  //Probability function = Number of occurances in each bin /
  //   Total Number of occurances.
  //
  // Returns pixels of float..
  typedef  TOutput OutputPixelType;

  HistogramEntropyFunction():
    m_TotalFrequency(1) {}

  ~HistogramEntropyFunction() {}

  inline OutputPixelType operator()(const TInput & A) const
  {
    if ( A )
      {
      const double p = static_cast< OutputPixelType >( A )
                       / static_cast< OutputPixelType >( m_TotalFrequency );
      return static_cast< OutputPixelType >( ( -1 ) * p * vcl_log(p) / vcl_log(2.0) );
      }
    else
      {
      const double p = static_cast< OutputPixelType >( A + 1 )
                       / static_cast< OutputPixelType >( m_TotalFrequency );
      return static_cast< OutputPixelType >( ( -1 ) * p * vcl_log(p) / vcl_log(2.0) );
      }
  }

  void SetTotalFrequency(const SizeValueType n)
  {
    m_TotalFrequency = n;
  }

  SizeValueType GetTotalFrequency() const
  {
    return m_TotalFrequency;
  }

private:
  SizeValueType m_TotalFrequency;
};
}

template< class THistogram, unsigned int NDimension, class TOutputPixel = double >
class ITK_EXPORT HistogramToEntropyImageFilter:
  public HistogramToImageFilter< THistogram, NDimension,
                                 Function::HistogramEntropyFunction< SizeValueType, TOutputPixel > >
{
public:

  /** Standard class typedefs. */
  typedef HistogramToEntropyImageFilter Self;

  /** Standard "Superclass" typedef. */
  typedef HistogramToImageFilter< THistogram, NDimension,
                                  Function::HistogramEntropyFunction< SizeValueType, TOutputPixel > >
  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(HistogramToEntropyImageFilter, HistogramToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
protected:
  HistogramToEntropyImageFilter() {}
  virtual ~HistogramToEntropyImageFilter() {}
private:
  HistogramToEntropyImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                //purposely not implemented
};
} // end namespace itk

#endif
