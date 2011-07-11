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
#ifndef __itkLabelMapContourOverlayImageFilter_txx
#define __itkLabelMapContourOverlayImageFilter_txx

#include "itkLabelMapContourOverlayImageFilter.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkObjectByObjectLabelMapFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkLabelUniqueLabelMapFilter.h"


namespace itk {

template<class TLabelMap, class TFeatureImage, class TOutputImage>
LabelMapContourOverlayImageFilter<TLabelMap, TFeatureImage, TOutputImage>
::LabelMapContourOverlayImageFilter()
{
  this->SetNumberOfRequiredInputs(2);
  m_Opacity = 0.5;
  m_Type = CONTOUR;
  m_Priority = HIGH_LABEL_ON_TOP;
  SizeType s;
  s.Fill( 1 );
  m_ContourThickness = SizeType( s );
  s.Fill( 0 );
  m_DilationRadius = SizeType( s );
  m_SliceDimension = ImageDimension - 1;
}

template<class TLabelMap, class TFeatureImage, class TOutputImage>
void
LabelMapContourOverlayImageFilter<TLabelMap, TFeatureImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // We need all the input.
  LabelMapPointer input = const_cast<LabelMapType *>(this->GetInput());
  if ( !input )
    { return; }
  input->SetRequestedRegion( input->GetLargestPossibleRegion() );
}

template <class TLabelMap, class TFeatureImage, class TOutputImage>
void
LabelMapContourOverlayImageFilter<TLabelMap, TFeatureImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()
    ->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template<class TLabelMap, class TFeatureImage, class TOutputImage>
void
LabelMapContourOverlayImageFilter<TLabelMap, TFeatureImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  typedef ObjectByObjectLabelMapFilter< LabelMapType, LabelMapType > OBOType;
  typename OBOType::Pointer obo = OBOType::New();
  obo->SetInput( this->GetInput() );
  SizeType rad = m_DilationRadius;
  for( int i=0; i<ImageDimension; i++ )
    {
    rad[i]++;
    }
  obo->SetPadSize( rad );
  // obo->SetInPlace( false );

  // dilate the image
  typedef typename OBOType::InternalInputImageType InternalImageType;
  typedef FlatStructuringElement< ImageDimension > KernelType;
  typedef BinaryDilateImageFilter< InternalImageType, InternalImageType, KernelType > DilateType;
  typename DilateType::Pointer dilate = DilateType::New();
  dilate->SetKernel( KernelType::Ball( m_DilationRadius ) );
  obo->SetInputFilter( dilate );

//   typedef typename CastImageFilter< InternalImageType, InternalImageType, KernelType > CastType;
//   typename CastType::Pointer cast = CastType::New();
//   cast->SetInPlace( false );

  typedef BinaryErodeImageFilter< InternalImageType, InternalImageType, KernelType > ErodeType;
  typename ErodeType::Pointer erode = ErodeType::New();
  erode->SetKernel( KernelType::Ball( m_ContourThickness ) );
  erode->SetInput( dilate->GetOutput() );

  typedef SubtractImageFilter< InternalImageType, InternalImageType > SubtractType;
  typename SubtractType::Pointer sub = SubtractType::New();
  sub->SetInput( 0, dilate->GetOutput() );
  sub->SetInput( 1, erode->GetOutput() );


  typedef SliceBySliceImageFilter< InternalImageType, InternalImageType > SliceType;
  typedef typename SliceType::InternalInputImageType SliceInternalImageType;
  typename SliceType::Pointer slice = SliceType::New();

  typedef CastImageFilter< SliceInternalImageType, SliceInternalImageType > SliceCastType;
  typename SliceCastType::Pointer scast = SliceCastType::New();
  scast->SetInPlace( false );
  slice->SetInputFilter( scast );

  typedef FlatStructuringElement< ImageDimension - 1 > SliceKernelType;
  typedef BinaryErodeImageFilter< SliceInternalImageType, SliceInternalImageType, SliceKernelType > SliceErodeType;
  typename SliceErodeType::Pointer serode = SliceErodeType::New();
  typedef typename SliceKernelType::RadiusType RadiusType;
  RadiusType srad;
  srad.Fill(NumericTraits<typename RadiusType::SizeValueType>::Zero);
  int j=0;
  for( int i=0; i<ImageDimension; i++ )
    {
    if( j != m_SliceDimension )
      {
      srad[j] = m_ContourThickness[i];
      j++;
      }
    }
  serode->SetKernel( SliceKernelType::Ball( srad ) );
  serode->SetInput( scast->GetOutput() );

  typedef SubtractImageFilter< SliceInternalImageType, SliceInternalImageType > SliceSubtractType;
  typename SliceSubtractType::Pointer ssub = SliceSubtractType::New();
  ssub->SetInput( 0, scast->GetOutput() );
  ssub->SetInput( 1, serode->GetOutput() );
  slice->SetOutputFilter( ssub );

  // search the contour, or not
  if( m_Type == PLAIN )
    {
    // nothing to do
    obo->SetOutputFilter( dilate );
    }
  else if( m_Type == CONTOUR )
    {
//     typedef BinaryContourImageFilter< InternalImageType, InternalImageType > ContourType;
//     typename ContourType::Pointer contour = ContourType::New();
//     contour->SetInput( dilate->GetOutput() );
//     obo->SetOutputFilter( contour );
     obo->SetOutputFilter( sub );
    }
  else if( m_Type == SLICE_CONTOUR )
    {
    slice->SetInput( dilate->GetOutput() );
    slice->SetDimension( m_SliceDimension );
    obo->SetOutputFilter( slice );

//     typedef typename SliceType::InternalInputImageType SliceInternalType;
//     typedef BinaryContourImageFilter< SliceInternalType, SliceInternalType > SliceContourType;
//     typename SliceContourType::Pointer slice_contour = SliceContourType::New();
//     slice->SetFilter( slice_contour );
    }
  else
    {
    itkExceptionMacro(<< "Unsupported Type: " << m_Type);
    }

  // choose which labels will be on top of the oters
  typedef LabelUniqueLabelMapFilter< LabelMapType > UniqueType;
  typename UniqueType::Pointer uniq = UniqueType::New();
  uniq->SetInput( obo->GetOutput() );
  uniq->SetReverseOrdering( m_Priority == LOW_LABEL_ON_TOP );

  m_TempImage = uniq->GetOutput();
  m_TempImage->Update();
  m_TempImage->DisconnectPipeline();


  ThreadIdType nbOfThreads = this->GetNumberOfThreads();
  if( itk::MultiThreader::GetGlobalMaximumNumberOfThreads() != 0 )
    {
    nbOfThreads = std::min( this->GetNumberOfThreads(), MultiThreader::GetGlobalMaximumNumberOfThreads() );
    }
  // number of threads can be constrained by the region size, so call the SplitRequestedRegion
  // to get the real number of threads which will be used
  typename TOutputImage::RegionType splitRegion;  // dummy region - just to call the following method
  nbOfThreads = this->SplitRequestedRegion(0, nbOfThreads, splitRegion);

  m_Barrier = Barrier::New();
  m_Barrier->Initialize( nbOfThreads );

  Superclass::BeforeThreadedGenerateData();

}


template<class TLabelMap, class TFeatureImage, class TOutputImage>
void
LabelMapContourOverlayImageFilter<TLabelMap, TFeatureImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId )
{
  OutputImageType * output = this->GetOutput();
  LabelMapType * input = const_cast<LabelMapType *>(this->GetInput());
  const FeatureImageType * input2 = this->GetFeatureImage();

  FunctorType function;
  function.SetBackgroundValue( input->GetBackgroundValue() );
  function.SetOpacity( m_Opacity );

  ImageRegionConstIterator< FeatureImageType > featureIt( input2, outputRegionForThread );
  ImageRegionIterator< OutputImageType > outputIt( output, outputRegionForThread );

  for ( featureIt.GoToBegin(), outputIt.GoToBegin();
        !featureIt.IsAtEnd();
        ++featureIt, ++outputIt )
    {
    outputIt.Set( function( featureIt.Get(), input->GetBackgroundValue() ) );
    }

  // wait for the other threads to complete that part
  m_Barrier->Wait();

  // and delegate to the superclass implementation to use the thread support for the label objects
  Superclass::ThreadedGenerateData( outputRegionForThread, threadId );
}


template<class TLabelMap, class TFeatureImage, class TOutputImage>
void
LabelMapContourOverlayImageFilter<TLabelMap, TFeatureImage, TOutputImage>
::ThreadedProcessLabelObject( LabelObjectType * labelObject )
{
  OutputImageType * output = this->GetOutput();
  LabelMapType * input = const_cast<LabelMapType *>(this->GetInput());
  const FeatureImageType * input2 = this->GetFeatureImage();

  FunctorType function;
  function.SetBackgroundValue( input->GetBackgroundValue() );
  function.SetOpacity( m_Opacity );

  const typename LabelObjectType::LabelType & label = labelObject->GetLabel();

  // the user want the mask to be the background of the label collection image
  typename LabelMapType::LabelObjectType::LineContainerType::const_iterator lit;
  typename LabelMapType::LabelObjectType::LineContainerType & lineContainer = labelObject->GetLineContainer();
  typedef typename LabelMapType::LabelObjectType::LengthType  LengthType;

  for( lit = lineContainer.begin(); lit != lineContainer.end(); lit++ )
    {
    IndexType idx = lit->GetIndex();
    LengthType length = lit->GetLength();
    for( LengthType i = 0; i < length; i++)
      {
      output->SetPixel( idx, function( input2->GetPixel(idx), label ) );
      idx[0]++;
      }
    }

}


template<class TLabelMap, class TFeatureImage, class TOutputImage>
void
LabelMapContourOverlayImageFilter<TLabelMap, TFeatureImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Opacity: " << m_Opacity << std::endl;
  os << indent << "Type: " << m_Type << std::endl;
  os << indent << "Priority: " << m_Priority << std::endl;
  os << indent << "ContourThickness: " << m_ContourThickness << std::endl;
  os << indent << "DilationRadius: " << m_DilationRadius << std::endl;
  os << indent << "SliceDimension: " << m_SliceDimension << std::endl;
}


}// end namespace itk
#endif
