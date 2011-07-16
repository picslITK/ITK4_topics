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
#ifndef __itkObjectToDataBase_h
#define __itkObjectToDataBase_h

#include "itkProcessObject.h"
#include "itkImage.h"

namespace itk
{

/** \class ObjectToDataBase
 *  \brief Virtual base class for threading setup and dispatch.
 *
 * The class is templated over the type of object over which threading
 * is performed, e.g. an image region. And it is templated over the
 * type of the data holder. The data holder is supplied to the threading
 * callback for the user, i.e. user data.
 *
 * SplitRequestedObject is a method to split the object into
 * non-overlapping pieces for threading. Must be overridden by derived
 * classes to provide the particular funcationality required for
 * TInputObject type.
 *
 * Call SetHolder to assign a user data object.
 *
 * Call SetOverallObject to assign the object to split into per-thread
 * regions.
 *
 * Call SetThreadedGenerateData to define the worker callback function,
 *  which is called from each thread with a unique region to process.
 * \warning This callback function must be \c static if it is a class method.
 *
 * Call GenerateData to begin threaded processing.
 *
 * \warning The actual number of threads used may be less than the
 * requested number of threads. Either because the requested number is
 * greater than the number available, or the SplitRequestedObject method
 * decides that fewer threads would be more efficient. After the threader
 * has run, m_NumberOfThreadsUsed holds the actual number used.
 * See \c DetermineNumberOfThreadsToUse to get the number of threads
 * before running.
 *
 * \note There is no test for this class yet. See Array1DToDataTest.
 * \note The name of this and derived classes should possibly be
 * changed for more clarity.
 *
 * \sa ImageToData
 * \ingroup ITK-RegistrationRefactoring
 */

template <class TInputObject, class TDataHolder>
class ITK_EXPORT ObjectToDataBase : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ObjectToDataBase          Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Smart Pointer type to a DataObject. */
  //typedef DataObject::Pointer DataObjectPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ObjectToDataBase,ProcessObject);

  /** Type of the input object that's split for threading */
  typedef TInputObject              InputObjectType;

  /** Type of the data holder, i.e. user data */
  typedef TDataHolder               DataHolderType;

  /** Type of callback function called by threader to do the work.
   */
  // http://www.parashift.com/c++-faq-lite/pointers-to-members.html
  //  [33.2] How do I pass a pointer-to-member-function to a signal handler,
  //  X event callback, system call that starts a thread/task, etc?
  //  Don't.
  //  so this call back function has to be define outside of the class OR use
  //  static function instead
  //
  //  Note: static member functions do not require an actual object to be
  //  invoked,
  //  so pointers-to-static-member-functions are usually type-compatible with
  //  regular pointers-to-functions. However, although it probably works on
  //  most compilers, it
  //  actually would have to be an extern "C" non-member function to be correct,
  //  since "C linkage" doesn't only cover things like name mangling, but also
  //  calling conventions, which might be different between C and C++.
  typedef void (*ThreadedGenerateDataFuncType)(const InputObjectType&,
                                                ThreadIdType threadId,
                                                DataHolderType * holder);

  /** Set the overall (i.e. complete) object over which to thread */
  void SetOverallObject( InputObjectType & object )
    {
    if( object != m_OverallObject )
      {
      m_OverallObject = object;
      m_OverallObjectHasBeenSet = true;
      this->Modified();
      }
    }

  /** Set the threaded worker callback. Used by the user class
   * to assign the worker callback.
   * \note This callback must be a static function if it is a class method. */
  void SetThreadedGenerateData( ThreadedGenerateDataFuncType func )
  { m_ThreadedGenerateData = func; }

  /** Set the object holder used during threading. */
  void SetHolder( DataHolderType* holder )
    {
    if( this->m_Holder != holder )
      {
      this->m_Holder = holder;
      this->Modified();
      }
    }

  /** Get the assigned holder, as a raw C-pointer */
  DataHolderType* GetHolder()
  { return this->m_Holder; }

  /** Accessor for number of threads actually used */
  itkGetMacro( NumberOfThreadsUsed, ThreadIdType );

  /** Start the threading process */
  virtual void GenerateData();

  /** Determine the number of threads that will be used by calling
   * SplitRequestedObject once and getting the return value. This uses
   * m_NumberOfThreads and requires that \c m_OverallObject has been set.
   * The number may be less than m_NumberOfThreads if SplitRequestedObject
   * determines that fewer threads would be more efficient. */
  ThreadIdType DetermineNumberOfThreadsToUse(void);

protected:
  ObjectToDataBase();
  virtual ~ObjectToDataBase() {}

  /** Split the output's RequestedObject into \c requestedTotal "pieces",
   * returning piece \c i as \c splitObject. "Pieces" may represent
   * an image region, or a index range for a parameter array, etc, depending
   * on the type of object over which this class is templated.
   * This method is called \c requestedTotal times, which is the number of
   * threads available for use. The
   * pieces must not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedObject,
   * i.e. return value is less than or equal to \c requestedTotal.
   * This must be overridden by derived classes to provide specialized
   * behavior. */
  virtual
  ThreadIdType SplitRequestedObject(ThreadIdType threadID,
                           ThreadIdType requestedTotal,
                           InputObjectType& overallObject,
                           InputObjectType& splitObject) const = 0;

  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to the function pointed to by m_ThreadedGenerateData. */
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );

  /** Internal structure used for passing image data into the
   * threading library */
  struct ThreadStruct
    {
    Pointer Filter;
    };

  /** The object over which to thread. */
  //This should probably be made into a SmartPointer
  InputObjectType               m_OverallObject;
  /** Flag is set when user calls SetOveralObject */
  bool                          m_OverallObjectHasBeenSet;
  /** Raw C-pointer to the data holder */
  DataHolderType *              m_Holder;

private:
  ThreadedGenerateDataFuncType  m_ThreadedGenerateData;
  /** Store the actual number of threads used, which may be less than
   * the number allocated by the threader if the object does not split
   * well into that number.
   * This value is valid once the threader has been run. */
  ThreadIdType                  m_NumberOfThreadsUsed;

  ObjectToDataBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkObjectToDataBase.txx"
#endif

#endif
