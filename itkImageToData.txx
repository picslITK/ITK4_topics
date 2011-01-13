/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageSource.txx,v $
  Language:  C++
  Date:      $Date: 2009-11-03 12:24:24 $
  Version:   $Revision: 1.69 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToData_txx
#define __itkImageToData_txx
#include "itkImageToData.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 *
 */
template <unsigned int VDimension, class TDataHolder>
ImageToData<VDimension, TDataHolder>
::ImageToData()
{
  // Create the output. We use static_cast<> here because we know the default
  // output must be of type TOutputImage
//  typename TOutputImage::Pointer output
//    = static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer());
  this->ProcessObject::SetNumberOfRequiredOutputs(0);
//  this->ProcessObject::SetNthOutput(0, output.GetPointer());

  // Set the default behavior of an image source to NOT release its
  // output bulk data prior to GenerateData() in case that bulk data
  // can be reused (an thus avoid a costly deallocate/allocate cycle).
//  this->ReleaseDataBeforeUpdateFlagOff();
}

/**
 *
 */
//template<class TOutputImage>
//typename ImageToData<TOutputImage>::DataObjectPointer
//ImageToData<TOutputImage>
//::MakeOutput(unsigned int)
//{
//  return static_cast<DataObject*>(TOutputImage::New().GetPointer());
//}

/**
 *
 */
//template<class TOutputImage>
//typename ImageToData<TOutputImage>::OutputImageType *
//ImageToData<TOutputImage>
//::GetOutput()
//{
//  if (this->GetNumberOfOutputs() < 1)
//    {
//    return 0;
//    }
//
//  // we assume that the first output is of the templated type
//  return static_cast<TOutputImage*>
//    (this->ProcessObject::GetOutput(0));
//}
//
//
///**
// *
// */
//template<class TOutputImage>
//typename ImageToData<TOutputImage>::OutputImageType *
//ImageToData<TOutputImage>
//::GetOutput(unsigned int idx)
//{
//  TOutputImage* out = dynamic_cast<TOutputImage*>
//    (this->ProcessObject::GetOutput(idx));
//  if ( out == NULL ) {
//    itkWarningMacro ( << "dynamic_cast to output type failed" );
//  }
//  return out;
//}
//
///**
// *
// */
//template<class TOutputImage>
//void
//ImageToData<TOutputImage>
//::GraftOutput(DataObject *graft)
//{
//  this->GraftNthOutput(0, graft);
//}
//
//
///**
// *
// */
//template<class TOutputImage>
//void
//ImageToData<TOutputImage>
//::GraftNthOutput(unsigned int idx, DataObject *graft)
//{
//  if ( idx >= this->GetNumberOfOutputs() )
//    {
//    itkExceptionMacro(<<"Requested to graft output " << idx <<
//        " but this filter only has " << this->GetNumberOfOutputs() << " Outputs.");
//    }
//
//  if ( !graft )
//    {
//    itkExceptionMacro(<<"Requested to graft output that is a NULL pointer" );
//    }
//
//  // we use the process object method since all out output may not be
//  // of the same type
//  DataObject * output = this->ProcessObject::GetOutput(idx);
//
//  // Call GraftImage to copy meta-information, regions, and the pixel container
//  output->Graft( graft );
//}

//----------------------------------------------------------------------------
template <unsigned int VDimension, class TDataHolder>
int
ImageToData<VDimension, TDataHolder>
::SplitRequestedRegion(int i, int num, ImageRegionType &overallRegion, ImageRegionType& splitRegion)
{
  // Get the output pointer
//  ImageType * outputPtr = this->GetOutput();
//  const typename TOutputImage::SizeType& requestedRegionSize
//    = outputPtr->GetRequestedRegion().GetSize();

   const SizeType requestedRegionSize = overallRegion.GetSize();


  int splitAxis;
  IndexType splitIndex;
  SizeType splitSize;

  // Initialize the splitRegion to the output requested region
  splitRegion = overallRegion;
  splitIndex = splitRegion.GetIndex();
  splitSize = splitRegion.GetSize();

  // split on the outermost dimension available
  splitAxis = this->ImageDimension - 1;
  while (requestedRegionSize[splitAxis] == 1)
    {
    --splitAxis;
    if (splitAxis < 0)
      { // cannot split
      itkDebugMacro("  Cannot Split");
      return 1;
      }
    }

  // determine the actual number of pieces that will be generated
  typename SizeType::SizeValueType range = requestedRegionSize[splitAxis];
  int valuesPerThread = Math::Ceil<int>(range/(double)num);
  int maxThreadIdUsed = Math::Ceil<int>(range/(double)valuesPerThread) - 1;

  // Split the region
  if (i < maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    splitSize[splitAxis] = valuesPerThread;
    }
  if (i == maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    // last thread needs to process the "rest" dimension being split
    splitSize[splitAxis] = splitSize[splitAxis] - i*valuesPerThread;
    }

  // set the split region ivars
  splitRegion.SetIndex( splitIndex );
  splitRegion.SetSize( splitSize );

  itkDebugMacro("  Split Piece: " << splitRegion );

  return maxThreadIdUsed + 1;
}


//----------------------------------------------------------------------------
template <unsigned int VDimension, class TDataHolder>
void
ImageToData<VDimension, TDataHolder>
::GenerateData()
{
  // Set up the multithreaded processing
  ThreadStruct str;
  str.Filter = this;

  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);

  // multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

}



// Callback routine used by the threading library. This routine just calls
// the ThreadedGenerateData method after setting the correct region for this
// thread.
template <unsigned int VDimension, class TDataHolder>
ITK_THREAD_RETURN_TYPE
ImageToData<VDimension, TDataHolder>
::ThreaderCallback( void *arg )
{
  ThreadStruct *str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  ImageRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount, str->Filter->m_OverallRegion,
                                            splitRegion);


  if (threadId < total)
    {
    str->Filter->ThreadedGenerateData(splitRegion, threadId, str->Filter->m_Holder);
    }
  // else
  //   {
  //   otherwise don't use this thread. Sometimes the threads dont
  //   break up very well and it is just as efficient to leave a
  //   few threads idle.
  //   }






  return ITK_THREAD_RETURN_VALUE;
}


} // end namespace itk

#endif
