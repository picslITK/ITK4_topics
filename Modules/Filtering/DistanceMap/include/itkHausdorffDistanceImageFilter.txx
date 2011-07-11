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
#ifndef __itkHausdorffDistanceImageFilter_txx
#define __itkHausdorffDistanceImageFilter_txx
#include "itkHausdorffDistanceImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkProgressAccumulator.h"
#include "itkDirectedHausdorffDistanceImageFilter.h"

namespace itk
{
template< class TInputImage1, class TInputImage2 >
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::HausdorffDistanceImageFilter()
{
  // this filter requires two input images
  this->SetNumberOfRequiredInputs(2);

  m_HausdorffDistance = NumericTraits< RealType >::Zero;
  m_AverageHausdorffDistance = NumericTraits< RealType >::Zero;
  m_UseImageSpacing = false;
}

template< class TInputImage1, class TInputImage2 >
void
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::SetInput1(const InputImage1Type *image)
{
  this->SetInput( image );
}

template< class TInputImage1, class TInputImage2 >
void
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::SetInput2(const TInputImage2 *image)
{
  this->SetNthInput( 1, const_cast< TInputImage2 * >( image ) );
}

template< class TInputImage1, class TInputImage2 >
const typename HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::InputImage1Type *
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::GetInput1()
{
  return this->GetInput();
}

template< class TInputImage1, class TInputImage2 >
const typename HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::InputImage2Type *
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::GetInput2()
{
  return static_cast< const TInputImage2 * >
         ( this->ProcessObject::GetInput(1) );
}

template< class TInputImage1, class TInputImage2 >
void
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // this filter requires:
  // - the largeset possible region of the first image
  // - the corresponding region of the second image
  if ( this->GetInput1() )
    {
    InputImage1Pointer image1 =
      const_cast< InputImage1Type * >( this->GetInput1() );
    image1->SetRequestedRegionToLargestPossibleRegion();

    if ( this->GetInput2() )
      {
      InputImage2Pointer image2 =
        const_cast< InputImage2Type * >( this->GetInput2() );
      image2->SetRequestedRegion(
        this->GetInput1()->GetRequestedRegion() );
      }
    }
}

template< class TInputImage1, class TInputImage2 >
void
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::EnlargeOutputRequestedRegion(DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}

template< class TInputImage1, class TInputImage2 >
void
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::GenerateData()
{
  ThreadIdType nbthreads = this->GetNumberOfThreads();

  // Pass the first input through as the output
  InputImage1Pointer image =
    const_cast< TInputImage1 * >( this->GetInput1() );

  this->GraftOutput(image);

  RealType distance12, distance21;

  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typedef DirectedHausdorffDistanceImageFilter< InputImage1Type, InputImage2Type >
  Filter12Type;

  typename Filter12Type::Pointer filter12 = Filter12Type::New();

  filter12->SetInput1( this->GetInput1() );
  filter12->SetInput2( this->GetInput2() );
  filter12->SetNumberOfThreads( nbthreads );

  typedef DirectedHausdorffDistanceImageFilter< InputImage2Type, InputImage1Type >
  Filter21Type;

  typename Filter21Type::Pointer filter21 = Filter21Type::New();

  filter21->SetInput1( this->GetInput2() );
  filter21->SetInput2( this->GetInput1() );
  filter21->SetNumberOfThreads( nbthreads );

  filter12->SetUseImageSpacing(m_UseImageSpacing);
  filter21->SetUseImageSpacing(m_UseImageSpacing);

  // Register the filter with the with progress accumulator using
  // equal weight proportion
  progress->RegisterInternalFilter(filter12, .5f);
  progress->RegisterInternalFilter(filter21, .5f);

  filter12->Update();
  distance12 = filter12->GetDirectedHausdorffDistance();
  filter21->Update();
  distance21 = filter21->GetDirectedHausdorffDistance();

  if ( distance12 > distance21 )
    {
    m_HausdorffDistance = distance12;
    }
  else
    {
    m_HausdorffDistance = distance21;
    }
  m_AverageHausdorffDistance =
    ( filter12->GetAverageHausdorffDistance() + filter21->GetAverageHausdorffDistance() ) / 2.0;
}

template< class TInputImage1, class TInputImage2 >
void
HausdorffDistanceImageFilter< TInputImage1, TInputImage2 >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "HausdorffDistance: "
     << m_HausdorffDistance << std::endl;
  os << indent << "AverageHausdorffDistance: "
     << m_AverageHausdorffDistance << std::endl;
  os << indent << "Use Image Spacing : "
     << m_UseImageSpacing << std::endl;
}
} // end namespace itk
#endif
