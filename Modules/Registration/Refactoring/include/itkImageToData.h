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
#ifndef __itkImageToData_h
#define __itkImageToData_h

#include "itkProcessObject.h"
#include "itkImage.h"

namespace itk
{

/** \class ImageToData
 *  \brief Base class for all process objects that output image data.
 *
 * ImageToData is the base class for all process objects that output
 * image data. Specifically, this class defines the GetOutput() method
 * that returns a pointer to the output image. The class also defines
 * some internal private data members that are used to manage streaming
 * of data.
 *
 * Memory management in an ImageToData is slightly different than a
 * standard ProcessObject.  ProcessObject's always release the bulk
 * data associated with their output prior to GenerateData() being
 * called. ImageSources default to not releasing the bulk data incase
 * that particular memory block is large enough to hold the new output
 * values.  This avoids unnecessary deallocation/allocation
 * sequences. ImageToData's can be forced to use a memory management
 * model similar to the default ProcessObject behaviour by calling
 * ProcessObject::ReleaseDataBeforeUpdateFlagOn().  A user may want to
 * set this flag to limit peak memory usage during a pipeline update.
 *
 * \ingroup DataSources
 */



template <unsigned int VDimension, class TDataHolder>
class ITK_EXPORT ImageToData : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ImageToData               Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);


  /** Smart Pointer type to a DataObject. */
  typedef DataObject::Pointer DataObjectPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToData,ProcessObject);

  /** Some convenient typedefs. */
  // typedef TImageRegion ImageRegionType;
  itkStaticConstMacro(ImageDimension, unsigned int, VDimension);

  typedef ImageRegion<VDimension> ImageRegionType;
  typedef Size<VDimension> SizeType;
  typedef Index<VDimension> IndexType;

// http://www.parashift.com/c++-faq-lite/pointers-to-members.html
//  [33.2] How do I pass a pointer-to-member-function to a signal handler, X event callback, system call that starts a thread/task, etc?
//  Don't.
//  so this call back function has to be define outside of the class OR use static function instead
//
//  Note: static member functions do not require an actual object to be invoked,
//  so pointers-to-static-member-functions are usually type-compatible with regular
//  pointers-to-functions. However, although it probably works on most compilers, it
//  actually would have to be an extern "C" non-member function to be correct,
//  since "C linkage" doesn't only cover things like name mangling, but also calling
//  conventions, which might be different between C and C++.
  typedef void (*ThreadedGenerateDataFuncType)(const ImageRegionType&,
                                    int ,
                                    TDataHolder *);

  ThreadedGenerateDataFuncType ThreadedGenerateData;
  virtual void GenerateData();


protected:
  ImageToData();
  virtual ~ImageToData() {}


  /** Split the output's RequestedRegion into "num" pieces, returning
   * region "i" as "splitRegion". This method is called "num" times. The
   * regions must not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedRegion,
   * i.e. return value is less than or equal to "num". */
  virtual
  int SplitRequestedRegion(int i, int num, ImageRegionType &overallRegion, ImageRegionType& splitRegion);



  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );

  /** Internal structure used for passing image data into the threading library */
  struct ThreadStruct
    {
    Pointer Filter;
    };

  // DataHolderType *m_Holder;

public:

  ImageRegionType m_OverallRegion;
  TDataHolder *m_Holder;

private:
  ImageToData(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
# include "itkImageToData.txx"
#endif

#endif
