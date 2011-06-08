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
template <unsigned int VDimension, typename TInputObject>
ImageToData<VDimension,TInputObject>
::ImageToData()
{
}

//----------------------------------------------------------------------------
template <unsigned int VDimension, typename TInputObject>
int
ImageToData<VDimension,TInputObject>
::SplitRequestedObject(int i, int requestedTotal,
                       InputObjectType &overallRegion,
                       InputObjectType& splitRegion) const
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
  int valuesPerThread = Math::Ceil<int>(range/(double)requestedTotal);
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

} // end namespace itk

#endif
