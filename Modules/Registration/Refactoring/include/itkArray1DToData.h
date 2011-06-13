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

/** \class Array1DToData
 *  \brief Class for custom threading setup and dispatch of 1D data.
 *
 * This class provides threading over 1D data.
 * The SplitRequestedObject method splits the provided IndexRange into
 * 1 or more inclusive sub-ranges.
 * IndexRange is treated as an inclusive range, and need
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
  typedef itk::Index<2>                      TIndexRange;
}

namespace itk
{
template<class TDataHolder>
class ITK_EXPORT Array1DToData
  : public ObjectToDataBase< TIndexRange, TDataHolder >
{
public:
  /** Standard class typedefs. */
  typedef Array1DToData                               Self;
  typedef ObjectToDataBase<TIndexRange, TDataHolder>  Superclass;
  typedef SmartPointer<Self>                          Pointer;
  typedef SmartPointer<const Self>                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Array1DToData,Superclass);

  /** Type of the object being threaded over */
  typedef TIndexRange                           IndexRangeType;
  /** Type for convenience of base class methods */
  typedef typename Superclass::InputObjectType  InputObjectType;

  /** Set the overall index range over which to operate.
   * This performs some error checking and is named more intuitively
   * for this derived class. */
  virtual void SetOverallIndexRange( IndexRangeType& range  );

protected:
  Array1DToData(); //use New() method instead of direct instantiation.
  virtual ~Array1DToData() {}

  /** Split the IndexRange \c overallIndexRange into
   * \c requestedTotal subranges, returning subrange \c i as \c splitIndex.
   * This method is called \c requestedTotal times. The
   * pieces will not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedObject,
   * i.e. return value is less than or equal to \c requestedTotal. */
  virtual
  ThreadIdType SplitRequestedObject(ThreadIdType i,
                           ThreadIdType requestedTotal,
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
