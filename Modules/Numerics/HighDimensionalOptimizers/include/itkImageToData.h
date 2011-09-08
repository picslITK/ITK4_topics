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
#ifndef __itkImageToData_h
#define __itkImageToData_h

#include "itkObjectToDataBase.h"
#include "itkImageRegion.h"

namespace itk
{

/** \class ImageToData
 *  \brief Class for custom threading setup and dispatch of Image data.
 *
 * This class provides threading over an ImageRegion. It provides a
 * \c SplitRequestedObject method that splits the provided ImageRegion into
 * subregions, along the z-axis. Additionally it makes setting up a threaded
 * operation easier.
 *
 * Call SetOverallRegion to define the ImageRegion over which to thread.
 * Call SetThreadedGenerateData to define the worker callback function,
 *  which is called from each thread with a unique region to process.
 *\warning This callback function must be \c static if it is a class method.
 *
 * Call \c SetHolder to provide a pointer to DataHolder object to
 * store arbitrary user data for use in the threader callback
 * (typically a class instance).
 *
 * Call GenerateData to begin threaded processing.
 *
 * This class is templated over image dimension and DataHolder type.
 * The third template parameter \c TInputObject should always be
 * left as default.
 *
 * \warning The actual number of threads used may be less than the
 * requested number of threads. Either because the requested number is
 * greater than the number available, or the SplitRequestedObject method
 * decides that fewer threads would be more efficient. After the threader
 * has run, m_NumberOfThreadsUsed holds the actual number used.
 * See \c DetermineNumberOfThreadsToUse to get the number of threads
 * before running.
 *
 * \sa ObjectToDataBase
 * \ingroup ITKHighDimensionalOptimizers
 */

template <unsigned int VDimension, class TDataHolder>
class ITK_EXPORT ImageToData
  : public ObjectToDataBase<ImageRegion<VDimension>, TDataHolder>
{
public:
  /** Standard class typedefs. */
  typedef ImageToData                                             Self;
  typedef ObjectToDataBase<ImageRegion<VDimension>, TDataHolder>  Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToData,ObjectToDataBase);

  /** Type of the object being threaded over */
  typedef typename Superclass::InputObjectType  InputObjectType;

  /** Some convenient typedefs. */
  // typedef TImageRegion ImageRegionType;
  itkStaticConstMacro(ImageDimension, unsigned int, VDimension);

  typedef ImageRegion<VDimension>   ImageRegionType;
  typedef Size<VDimension>          SizeType;
  typedef Index<VDimension>         IndexType;

  /** Set the overall image region over which to operate.
   * This is equivalent to SetOverallObject, but named more intuitively
   * for this derived class. */
  virtual void SetOverallRegion(  const ImageRegionType& region );

protected:
  ImageToData();
  virtual ~ImageToData();

  /** Split the ImageRegion \c overallRegion into \c requestedTotal subregions,
   * returning subregion \c i as \c splitRegion.
   * This method is called \c requestedTotal times. The
   * pieces must not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedObject,
   * i.e. return value is less than or equal to \c requestedTotal. */
  virtual
  ThreadIdType SplitRequestedObject(const ThreadIdType i,
                           const ThreadIdType requestedTotal,
                           const InputObjectType& overallRegion,
                           InputObjectType& splitRegion) const;

private:

  ImageToData(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkImageToData.hxx"
#endif

#endif
