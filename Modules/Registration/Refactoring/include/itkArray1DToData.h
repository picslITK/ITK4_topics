/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageSource.h,v $
  Language:  C++
  Date:      $Date: 2009-03-12 01:11:08 $
  Version:   $Revision: 1.59 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkArray1DToData_h
#define __itkArray1DToData_h

#include "itkObjectToDataBase.h"

namespace itk
{

/** \class Array1DToData
 *  \brief Class for custom threading setup and dispatch of 1D data.
 *
 * This class provides threading over 1D data.
 * The SplitRequestedObject method splits the provided IndexRange into
 * 1 or more inclusive sub-ranges.
 * IndexRange is expected to be an inclusive range, and need
 * not start at 0.
 *
 * Call SetOverallIndexRange to define the IndexRange over which to thread.
 * Call SetThreadedGenerateData to define the worker callback function,
 *  which during threading is called from each thread with a unique range
 *  to process.
 * Call SetHolder to provide a class instance...
 *
 * \sa ImageToData
 * \sa ObjectToDataBase
 * \ingroup DataSources
 */

namespace
{
  typedef Index<2>                           TIndexRange;
}

class ITK_EXPORT Array1DToData
  : public ObjectToDataBase< TIndexRange >
{
public:
  /** Standard class typedefs. */
  typedef Array1DToData                      Self;
  typedef ObjectToDataBase<IndexRangeType>   Superclass;
  typedef SmartPointer<Self>                 Pointer;
  typedef SmartPointer<const Self>           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Array1DToData,Superclass);

  /** Type of the object being threaded over */
  typedef typename Superclass::InputObjectType  InputObjectType;
  typedef TIndexRange                           IndexRangeType;

  /** Set the overall image region over which to operate.
   * This is equivalent to SetOverallObject, but named more intuitively
   * for this derived class. */
  virtual void SetOverallIndexRange(  IndexRangeType& range )
  {
    if( range[0] > range[1] )
      {
      itkExceptionMacro("Error in range. Begin is less than End: "
                        << range << ".");
      }
    this->SetOverallObject( range );
  }

protected:
  Array1DToData();
  virtual ~Array1DToData() {}

  /** Split the ImageRegion \c overallRegion into \c requestedTotal subregions,
   * returning subregion \c i as \c splitRegion.
   * This method is called \c requestedTotal times. The
   * pieces must not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedObject,
   * i.e. return value is less than or equal to \c requestedTotal. */
  virtual
  int SplitRequestedObject(int i,
                           int requestedTotal,
                           InputObjectType& overallIndexRange,
                           InputObjectType& splitIndexRange) const;

private:

  Array1DToData(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
# include "itkArray1DToData.txx"
#endif

#endif
