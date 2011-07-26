/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageSource.hxx,v $
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
#ifndef __itkObjectToDataBase_hxx
#define __itkObjectToDataBase_hxx
#include "itkObjectToDataBase.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 *
 */
template <class TInputObject, class TDataHolder>
ObjectToDataBase<TInputObject, TDataHolder>
::ObjectToDataBase()
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

  this->m_Holder = NULL;
  this->m_ThreadedGenerateData = NULL;
  this->m_OverallObjectHasBeenSet = false;
  this->m_NumberOfThreadsUsed = 0;
}

//----------------------------------------------------------------------------
template <class TInputObject, class TDataHolder>
void
ObjectToDataBase<TInputObject, TDataHolder>
::GenerateData()
{
  // Make sure hold and threader worker have been defined
  if( this->m_Holder == NULL )
    {
    itkExceptionMacro("m_Holder must be defined.");
    }
  if( this->m_ThreadedGenerateData == NULL )
    {
    itkExceptionMacro("m_ThreadedGenereateData must be defined.");
    }
  if( ! this->m_OverallObjectHasBeenSet )
    {
    itkExceptionMacro("m_OverallObject must be set.");
    }

  // Set up the multithreaded processing
  ThreadStruct str;
  str.Filter = this;

  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);

  // multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();
}

//----------------------------------------------------------------------------
// Callback routine used by the threading library. This routine just calls
// the ThreadedGenerateData method after setting the correct region for this
// thread.
template <class TInputObject, class TDataHolder>
ITK_THREAD_RETURN_TYPE
ObjectToDataBase<TInputObject, TDataHolder>
::ThreaderCallback( void *arg )
{
  ThreadStruct *str;
  ThreadIdType total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // first find out how many pieces extent can be split into.
  InputObjectType splitObject;
  total = str->Filter->SplitRequestedObject(threadId,
                                            threadCount,
                                            str->Filter->m_OverallObject,
                                            splitObject);

  // store the actual number of threads used
  if( threadId == 0 )
    {
    str->Filter->m_NumberOfThreadsUsed = total;
    }

  // execute the actual method with appropriate output region
  if (threadId < total)
    {
    str->Filter->m_ThreadedGenerateData(splitObject,
                                        threadId,
                                        str->Filter->m_Holder);
    }
  // else
  //   {
  //   otherwise don't use this thread. Sometimes the threads dont
  //   break up very well and it is just as efficient to leave a
  //   few threads idle.
  //   }

  return ITK_THREAD_RETURN_VALUE;
}

template <class TInputObject, class TDataHolder>
ThreadIdType
ObjectToDataBase<TInputObject, TDataHolder>
::DetermineNumberOfThreadsToUse(void)
{
  if( ! this->m_OverallObjectHasBeenSet )
    {
    itkExceptionMacro("m_OverallObject must be set.");
    }
  InputObjectType splitObject;
  return this->SplitRequestedObject( 0,
                                    this->GetNumberOfThreads(),
                                    this->m_OverallObject,
                                    splitObject);
}

} // end namespace itk

#endif
