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
#ifndef __itkCannyEdgeDetectionImageFilter_hxx
#define __itkCannyEdgeDetectionImageFilter_hxx
#include "itkCannyEdgeDetectionImageFilter.h"

#include "itkZeroCrossingImageFilter.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >::CannyEdgeDetectionImageFilter()
{
  unsigned int i;

  m_Variance.Fill(0.0);
  m_MaximumError.Fill(0.01);

  m_UpperThreshold = NumericTraits< OutputImagePixelType >::Zero;
  m_LowerThreshold = NumericTraits< OutputImagePixelType >::Zero;

  m_GaussianFilter      = GaussianImageFilterType::New();
  m_MultiplyImageFilter = MultiplyImageFilterType::New();
  m_UpdateBuffer1  = OutputImageType::New();

  // Set up neighborhood slices for all the dimensions.
  typename Neighborhood< OutputImagePixelType, ImageDimension >::RadiusType r;
  r.Fill(1);

  // Dummy neighborhood used to set up the slices.
  Neighborhood< OutputImagePixelType, ImageDimension > it;
  it.SetRadius(r);

  // Slice the neighborhood
  m_Center =  it.Size() / 2;

  for ( i = 0; i < ImageDimension; ++i )
    {
    m_Stride[i] = it.GetStride(i);
    }

  for ( i = 0; i < ImageDimension; ++i )
    {
    m_ComputeCannyEdgeSlice[i] =
      std::slice(m_Center - m_Stride[i], 3, m_Stride[i]);
    }

  // Allocate the derivative operator.
  m_ComputeCannyEdge1stDerivativeOper.SetDirection(0);
  m_ComputeCannyEdge1stDerivativeOper.SetOrder(1);
  m_ComputeCannyEdge1stDerivativeOper.CreateDirectional();

  m_ComputeCannyEdge2ndDerivativeOper.SetDirection(0);
  m_ComputeCannyEdge2ndDerivativeOper.SetOrder(2);
  m_ComputeCannyEdge2ndDerivativeOper.CreateDirectional();

  //Initialize the list
  m_NodeStore = ListNodeStorageType::New();
  m_NodeList = ListType::New();
}

template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::AllocateUpdateBuffer()
{
  // The update buffer looks just like the input.

  typename TInputImage::ConstPointer input = this->GetInput();

  m_UpdateBuffer1->CopyInformation(input);
  m_UpdateBuffer1->SetRequestedRegion( input->GetRequestedRegion() );
  m_UpdateBuffer1->SetBufferedRegion( input->GetBufferedRegion() );
  m_UpdateBuffer1->Allocate();
}

template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
throw( InvalidRequestedRegionError )
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  return;
  // get pointers to the input and output
  typename Superclass::InputImagePointer inputPtr =
    const_cast< TInputImage * >( this->GetInput() );
  typename Superclass::OutputImagePointer outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  //Set the kernel size.
  // fix me this needs to be based on the variance of the gaussian filter
  typename TInputImage::SizeValueType radius = 1;

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius(radius);

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() ) )
    {
    inputPtr->SetRequestedRegion(inputRequestedRegion);
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion(inputRequestedRegion);

    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    std::ostringstream          msg;
    msg << this->GetNameOfClass()
        << "::GenerateInputRequestedRegion()";
    e.SetLocation( msg.str().c_str() );
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}

template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::ThreadedCompute2ndDerivative(const OutputImageRegionType &
                               outputRegionForThread, ThreadIdType threadId)
{
  ZeroFluxNeumannBoundaryCondition< TInputImage > nbc;

  ImageRegionIterator< TOutputImage > it;

  void *globalData = 0;

  // Here input is the result from the gaussian filter
  //      output is the update buffer.
  typename  OutputImageType::Pointer input  = m_GaussianFilter->GetOutput();
  typename  OutputImageType::Pointer output = this->GetOutput();

  // set iterator radius
  Size< ImageDimension > radius; radius.Fill(1);

  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage >::
  FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage > bC;
  faceList = bC(input, outputRegionForThread, radius);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage >::
  FaceListType::iterator fit;

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels(), 100, 0.0f, 0.5f);

  // Process the non-boundady region and then each of the boundary faces.
  // These are N-d regions which border the edge of the buffer.
  for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    NeighborhoodType bit(radius, input, *fit);

    it = ImageRegionIterator< OutputImageType >(output, *fit);
    bit.OverrideBoundaryCondition(&nbc);
    bit.GoToBegin();

    while ( !bit.IsAtEnd() )
      {
      it.Value() = ComputeCannyEdge(bit, globalData);
      ++bit;
      ++it;
      progress.CompletedPixel();
      }
    }
}

template< class TInputImage, class TOutputImage >
typename CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::OutputImagePixelType
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::ComputeCannyEdge( const NeighborhoodType & it,
                    void *itkNotUsed(globalData) )
{

  NeighborhoodInnerProduct< OutputImageType > innerProduct;

  OutputImagePixelType dx[ImageDimension];
  OutputImagePixelType dxx[ImageDimension];
  OutputImagePixelType dxy[ImageDimension * ( ImageDimension - 1 ) / 2];

  //  double alpha = 0.01;

  //Calculate 1st & 2nd order derivative
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    dx[i] = innerProduct(m_ComputeCannyEdgeSlice[i], it,
                         m_ComputeCannyEdge1stDerivativeOper);
    dxx[i] = innerProduct(m_ComputeCannyEdgeSlice[i], it,
                          m_ComputeCannyEdge2ndDerivativeOper);
    }

  OutputImagePixelType deriv = NumericTraits< OutputImagePixelType >::Zero;

  int k = 0;
  //Calculate the 2nd derivative
  for ( unsigned int i = 0; i < ImageDimension - 1; i++ )
    {
    for ( unsigned int j = i + 1; j < ImageDimension; j++ )
      {
      dxy[k] = 0.25 * it.GetPixel(m_Center - m_Stride[i] - m_Stride[j])
               - 0.25 * it.GetPixel(m_Center - m_Stride[i] + m_Stride[j])
               - 0.25 * it.GetPixel(m_Center + m_Stride[i] - m_Stride[j])
               + 0.25 * it.GetPixel(m_Center + m_Stride[i] + m_Stride[j]);

      deriv += 2.0 * dx[i] * dx[j] * dxy[k];
      k++;
      }
    }

  OutputImagePixelType gradMag = static_cast<OutputImagePixelType>(0.0001); // alpha * alpha;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    deriv += dx[i] * dx[i] * dxx[i];
    gradMag += dx[i] * dx[i];
    }

  deriv = deriv / gradMag;

  return deriv;
}

// Calculate the second derivative
template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::Compute2ndDerivative()
{
  CannyThreadStruct str;

  str.Filter = this;

  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->Compute2ndDerivativeThreaderCallback, &str);

  this->GetMultiThreader()->SingleMethodExecute();
}

template< class TInputImage, class TOutputImage >
ITK_THREAD_RETURN_TYPE
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::Compute2ndDerivativeThreaderCallback(void *arg)
{
  CannyThreadStruct *str;

  ThreadIdType total, threadId, threadCount;

  threadId = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
  threadCount = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;

  str = (CannyThreadStruct *)( ( (MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  OutputImageRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  if ( threadId < total )
    {
    str->Filter->ThreadedCompute2ndDerivative(splitRegion, threadId);
    }

  return ITK_THREAD_RETURN_VALUE;
}

template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  // Allocate the output
  this->GetOutput()->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );
  this->GetOutput()->Allocate();

  typename  InputImageType::ConstPointer input  = this->GetInput();

  typename ZeroCrossingImageFilter< TOutputImage, TOutputImage >::Pointer
  zeroCrossFilter = ZeroCrossingImageFilter< TOutputImage, TOutputImage >::New();

  typename GradientMagnitudeImageFilter< TOutputImage, TOutputImage >::Pointer
  gradMag = GradientMagnitudeImageFilter< TOutputImage, TOutputImage >::New();

  typename MultiplyImageFilter< TOutputImage, TOutputImage, TOutputImage >::Pointer multFilter =
    MultiplyImageFilter< TOutputImage, TOutputImage, TOutputImage >::New();

  this->AllocateUpdateBuffer();

  // 1.Apply the Gaussian Filter to the input image.-------
  m_GaussianFilter->SetVariance(m_Variance);
  m_GaussianFilter->SetMaximumError(m_MaximumError);
  m_GaussianFilter->SetInput(input);
  // modify to force excution, due to grafting complications
  m_GaussianFilter->Modified();
  m_GaussianFilter->Update();

  //2. Calculate 2nd order directional derivative-------
  // Calculate the 2nd order directional derivative of the smoothed image.
  // The output of this filter will be used to store the directional
  // derivative.
  this->Compute2ndDerivative();

  this->Compute2ndDerivativePos();

  // 3. Non-maximum suppression----------

  // Calculate the zero crossings of the 2nd directional derivative and write
  // the result to output buffer.
  zeroCrossFilter->SetInput( this->GetOutput() );
  zeroCrossFilter->Update();

  // 4. Hysteresis Thresholding---------

  // First get all the edges corresponding to zerocrossings
  m_MultiplyImageFilter->SetInput1(m_UpdateBuffer1);
  m_MultiplyImageFilter->SetInput2( zeroCrossFilter->GetOutput() );

  // To save memory, we will graft the output of the m_GaussianFilter,
  // which is no longer needed, into the m_MultiplyImageFilter.
  m_MultiplyImageFilter->GraftOutput( m_GaussianFilter->GetOutput() );
  m_MultiplyImageFilter->Update();

  //Then do the double threshoulding upon the edge reponses
  this->HysteresisThresholding();
}

template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::HysteresisThresholding()
{
  // This is the Zero crossings of the Second derivative multiplied with the
  // gradients of the image. HysteresisThresholding of this image should give
  // the Canny output.
  typename OutputImageType::Pointer input = m_MultiplyImageFilter->GetOutput();
  float value;

  ListNodeType *node;

// fix me
  ImageRegionIterator< TOutputImage > oit( input, input->GetRequestedRegion() );

  oit.GoToBegin();

  ImageRegionIterator< TOutputImage > uit( this->GetOutput(),
                                           this->GetOutput()->GetRequestedRegion() );
  uit.GoToBegin();
  while ( !uit.IsAtEnd() )
    {
    uit.Value() = NumericTraits< OutputImagePixelType >::Zero;
    ++uit;
    }

  while ( !oit.IsAtEnd() )
    {
    value = oit.Value();

    if ( value > m_UpperThreshold )
      {
      node = m_NodeStore->Borrow();
      node->m_Value = oit.GetIndex();
      m_NodeList->PushFront(node);
      FollowEdge( oit.GetIndex() );
      }

    ++oit;
    }
}

template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::FollowEdge(IndexType index)
{
  // This is the Zero crossings of the Second derivative multiplied with the
  // gradients of the image. HysteresisThresholding of this image should give
  // the Canny output.
  typename OutputImageType::Pointer input = m_MultiplyImageFilter->GetOutput();
  InputImageRegionType inputRegion = input->GetRequestedRegion();

  IndexType     nIndex;
  IndexType     cIndex;
  ListNodeType *node;

  //assign iterator radius
  Size< ImageDimension > radius;
  radius.Fill(1);

  ConstNeighborhoodIterator< TOutputImage > oit( radius,
                                                 input, input->GetRequestedRegion() );
  ImageRegionIteratorWithIndex< TOutputImage > uit( this->GetOutput(),
                                                    this->GetOutput()->GetRequestedRegion() );

  uit.SetIndex(index);
  if ( uit.Get() == NumericTraits< OutputImagePixelType >::One )
    {
    // we must remove the node if we are not going to follow it!

    // Pop the front node from the list and read its index value.
    node = m_NodeList->Front(); // get a pointer to the first node
    m_NodeList->PopFront();     // unlink the front node
    m_NodeStore->Return(node);  // return the memory for reuse
    return;
    }

  int nSize = m_Center * 2 + 1;
  while ( !m_NodeList->Empty() )
    {
    // Pop the front node from the list and read its index value.
    node = m_NodeList->Front(); // get a pointer to the first node
    cIndex = node->m_Value;     // read the value of the first node
    m_NodeList->PopFront();     // unlink the front node
    m_NodeStore->Return(node);  // return the memory for reuse

    // Move iterators to the correct index position.
    oit.SetLocation(cIndex);
    uit.SetIndex(cIndex);
    uit.Value() = 1;

    // Search the neighbors for new indicies to add to the list.
    for ( int i = 0; i < nSize; i++ )
      {
      nIndex = oit.GetIndex(i);
      uit.SetIndex(nIndex);
      if ( inputRegion.IsInside(nIndex) )
        {
        if ( oit.GetPixel(i) > m_LowerThreshold && uit.Value() != NumericTraits< OutputImagePixelType >::One  )
          {
          node = m_NodeStore->Borrow();  // get a new node struct
          node->m_Value = nIndex;        // set its value
          m_NodeList->PushFront(node);   // add the new node to the list

          uit.SetIndex(nIndex);
          uit.Value() = NumericTraits< OutputImagePixelType >::One;
          }
        }
      }
    }
}

template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::ThreadedCompute2ndDerivativePos(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId)
{
  ZeroFluxNeumannBoundaryCondition< TInputImage > nbc;

  ConstNeighborhoodIterator< TInputImage > bit;
  ConstNeighborhoodIterator< TInputImage > bit1;

  ImageRegionIterator< TOutputImage > it;

  // Here input is the result from the gaussian filter
  //      input1 is the 2nd derivative result
  //      output is the gradient of 2nd derivative
  typename OutputImageType::Pointer input1 = this->GetOutput();
  typename OutputImageType::Pointer input = m_GaussianFilter->GetOutput();

  typename  InputImageType::Pointer output  = m_UpdateBuffer1;

  // set iterator radius
  Size< ImageDimension > radius; radius.Fill(1);

  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage >::
  FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage > bC;
  faceList = bC(input, outputRegionForThread, radius);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage >::
  FaceListType::iterator fit;

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels(), 100, 0.5f, 0.5f);

  InputImagePixelType zero = NumericTraits< InputImagePixelType >::Zero;

  OutputImagePixelType dx[ImageDimension];
  OutputImagePixelType dx1[ImageDimension];

  OutputImagePixelType directional[ImageDimension];
  OutputImagePixelType derivPos;

  OutputImagePixelType gradMag;

  // Process the non-boundary region and then each of the boundary faces.
  // These are N-d regions which border the edge of the buffer.

  NeighborhoodInnerProduct< OutputImageType > IP;

  for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    bit = ConstNeighborhoodIterator< InputImageType >(radius,
                                                      input, *fit);
    bit1 = ConstNeighborhoodIterator< InputImageType >(radius,
                                                       input1, *fit);
    it = ImageRegionIterator< OutputImageType >(output, *fit);
    bit.OverrideBoundaryCondition(&nbc);
    bit.GoToBegin();
    bit1.GoToBegin();
    it.GoToBegin();

    while ( !bit.IsAtEnd() )
      {
      gradMag = 0.0001;

      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        dx[i] = IP(m_ComputeCannyEdgeSlice[i], bit,
                   m_ComputeCannyEdge1stDerivativeOper);
        gradMag += dx[i] * dx[i];

        dx1[i] = IP(m_ComputeCannyEdgeSlice[i], bit1,
                    m_ComputeCannyEdge1stDerivativeOper);
        }

      gradMag = vcl_sqrt( (double)gradMag );
      derivPos = zero;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        //First calculate the directional derivative

        directional[i] = dx[i] / gradMag;

        //calculate gradient of 2nd derivative

        derivPos += dx1[i] * directional[i];
        }

      it.Value() = ( ( derivPos <= zero ) );
      it.Value() = it.Get() * gradMag;
      ++bit;
      ++bit1;
      ++it;
      progress.CompletedPixel();
      }
    }
}

//Calculate the second derivative
template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::Compute2ndDerivativePos()
{
  CannyThreadStruct str;

  str.Filter = this;

  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->Compute2ndDerivativePosThreaderCallback, &str);

  this->GetMultiThreader()->SingleMethodExecute();
}

template< class TInputImage, class TOutputImage >
ITK_THREAD_RETURN_TYPE
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::Compute2ndDerivativePosThreaderCallback(void *arg)
{
  CannyThreadStruct *str;

  ThreadIdType total, threadId, threadCount;

  threadId = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
  threadCount = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;

  str = (CannyThreadStruct *)( ( (MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.

  OutputImageRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  if ( threadId < total )
    {
    str->Filter->ThreadedCompute2ndDerivativePos(splitRegion, threadId);
    }

  return ITK_THREAD_RETURN_VALUE;
}

template< class TInputImage, class TOutputImage >
void
CannyEdgeDetectionImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "Variance: " << m_Variance << std::endl;
  os << "MaximumError: " << m_MaximumError << std::endl;
  os << indent << "UpperThreshold: "
     << static_cast< typename NumericTraits< OutputImagePixelType >::PrintType >( m_UpperThreshold )
     << std::endl;
  os << indent << "LowerThreshold: "
     << static_cast< typename NumericTraits< OutputImagePixelType >::PrintType >( m_LowerThreshold )
     << std::endl;
  os << "Center: "
     << m_Center << std::endl;
  os << "Stride: "
     << m_Stride << std::endl;
  os << "Gaussian Filter: " << std::endl;
  m_GaussianFilter->Print( os, indent.GetNextIndent() );
  os << "Multiply image Filter: " << std::endl;
  m_MultiplyImageFilter->Print( os, indent.GetNextIndent() );
  os << "UpdateBuffer1: " << std::endl;
  m_UpdateBuffer1->Print( os, indent.GetNextIndent() );
}
} //end of itk namespace
#endif
