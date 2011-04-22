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
#ifndef __itkHistogramAlgorithmBase_h
#define __itkHistogramAlgorithmBase_h

#include "itkMacro.h"
#include "itkObjectFactory.h"
#include "itkObject.h"

namespace itk
{
/** \class HistogramAlgorithmBase
 * \brief base class for algorithms operating on histograms
 *
 * You plug in the target sample data using SetInputHistogram method. Then call
 * the GenerateData method to run the alogithm.
 *
 * \ingroup ITK-ImageStatistics
 */

template< class TInputHistogram >
class ITK_EXPORT HistogramAlgorithmBase:public Object
{
public:
  /**Standard class typedefs. */
  typedef HistogramAlgorithmBase     Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /**Standard Macros */
  itkTypeMacro(HistogramAlgorithmBase, Object);

  /** Histogram typedefs alias */
  typedef TInputHistogram InputHistogramType;

  /** Stores the histogram pointer */
  void SetInputHistogram(const TInputHistogram *histogram)
  {
    if ( m_InputHistogram != histogram )
      {
      m_InputHistogram = histogram;
      this->Modified();
      }
  }

  /** Returns the histogram const pointer */
  const TInputHistogram * GetInputHistogram() const
  { return m_InputHistogram.GetPointer(); }

  /** dummy function that calls the GenerateData() function to generate
   * output. It exists for future compatibility with ProcessObject
   * without streaming */
  void Update()
  { this->GenerateData(); }
protected:
  HistogramAlgorithmBase();
  virtual ~HistogramAlgorithmBase() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  virtual void GenerateData() = 0;

private:
  /** Target histogram data pointer */
  typename TInputHistogram::ConstPointer m_InputHistogram;
}; // end of class
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHistogramAlgorithmBase.txx"
#endif

#endif
