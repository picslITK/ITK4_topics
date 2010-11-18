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
#ifndef __itkObjectToObjectMetric_txx
#define __itkObjectToObjectMetric_txx

#include "itkObjectToObjectMetric.h"

namespace itk
{
/**
 * Constructor
 */
template< class TFixedImage, class TMovingImage >
ObjectToObjectMetric< TFixedImage, TMovingImage >
::ObjectToObjectMetric()
{
  m_Threader = MultiThreaderType::New();
  m_ThreaderParameter.metric = this;
  m_WithinThreadPreProcess = false;
  m_WithinThreadPostProcess = false;
  m_NumberOfThreads = m_Threader->GetNumberOfThreads();
}

template< class TFixedImage, class TMovingImage >
ObjectToObjectMetric< TFixedImage, TMovingImage >
::~ObjectToObjectMetric()
{
}

/**
 * Set the number of threads. This will be clamped by the
 * multithreader, so we must check to see if it is accepted.
 */
template< class TFixedImage, class TMovingImage >
void
ObjectToObjectMetric< TFixedImage, TMovingImage >
::SetNumberOfThreads(unsigned int numberOfThreads)
{
  m_Threader->SetNumberOfThreads(numberOfThreads);
  m_NumberOfThreads = m_Threader->GetNumberOfThreads();
}


/**
 * PrintSelf
 */
template< class TFixedImage, class TMovingImage >
void
ObjectToObjectMetric< TFixedImage, TMovingImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


} // end namespace itk

#endif
