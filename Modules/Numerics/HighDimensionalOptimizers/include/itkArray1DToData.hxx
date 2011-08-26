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
#ifndef __itkArray1DToData_hxx
#define __itkArray1DToData_hxx

#include "itkArray1DToData.h"
#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Default constructor
 */
template<class TDataHolder>
Array1DToData<TDataHolder>::Array1DToData()
{
  this->m_OverallObject.Fill(0);
}

/**
 * Set the overall range over which to thread.
 */
template<class TDataHolder>
void
Array1DToData<TDataHolder>
::SetOverallIndexRange(  IndexRangeType& range )
{
  if( range[0] > range[1] )
    {
    itkExceptionMacro("Error in range.  Begin is less than End: "
                      << range << ".");
    }
  this->SetOverallObject( range );
}

/**
 * Split the requested range into a subrange.
 */
template<class TDataHolder>
ThreadIdType
Array1DToData<TDataHolder>
::SplitRequestedObject(ThreadIdType i, ThreadIdType requestedTotal,
                       InputObjectType& overallIndexRange,
                       InputObjectType& splitIndexRange) const
{
  // overallIndexRange is expected to be inclusive

  // determine the actual number of pieces that will be generated
  IndexRangeType::IndexValueType count =
    overallIndexRange[1] - overallIndexRange[0] + 1;
  ThreadIdType valuesPerThread = Math::Ceil<ThreadIdType>(count/(double)requestedTotal);
  ThreadIdType maxThreadIdUsed = Math::Ceil<ThreadIdType>(count/(double)valuesPerThread) - 1;

  // Split the index range
  if (i < maxThreadIdUsed)
    {
    splitIndexRange[0] = overallIndexRange[0] + i * valuesPerThread;
    splitIndexRange[1] = splitIndexRange[0] + valuesPerThread - 1;
    }
  if (i == maxThreadIdUsed)
    {
    splitIndexRange[0] = overallIndexRange[0] + i * valuesPerThread;
    // last thread needs to process the "rest" of the range
    splitIndexRange[1] = overallIndexRange[1];
    }

  itkDebugMacro("Array1DToData:  Split : " << splitIndexRange );

  return maxThreadIdUsed + 1;
}

} // end namespace itk

#endif
