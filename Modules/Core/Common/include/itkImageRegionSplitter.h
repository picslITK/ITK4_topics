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
#ifndef __itkImageRegionSplitter_h
#define __itkImageRegionSplitter_h

#include "itkImageRegion.h"
#include "itkObjectFactory.h"

namespace itk
{
/** \class ImageRegionSplitter
 * \brief Divide a region into several pieces.
 *
 * ImageRegionSplitter divides an ImageRegion into smaller regions.
 * ImageRegionSplitter is used by the StreamingImageFilter to divide a
 * requested output region into a series of smaller requests of the
 * pipeline.  This object has two basic methods: GetNumberOfSplits()
 * and GetSplit().
 *
 * GetNumberOfSplits() is used to determine how may subregions a given
 * region can be divided.  You call GetNumberOfSplits with an argument
 * that is the number of subregions you want.  If the image region can
 * support that number of subregions, that number is returned.
 * Otherwise, the maximum number of splits less then or equal to the
 * argumen  be returned.  For example, if a region splitter class only divides
 * a region into horizontal slabs, then the maximum number of splits
 * will be the number of rows in the region.
 *
 * GetSplit() returns the ith of N subregions (as an ImageRegion object).
 *
 * This ImageRegionSplitter class divides a region along the outermost
 * dimension. If the outermost dimension has size 1 (i.e. a volume
 * with a single slice), the ImageRegionSplitter will divide the
 * region along the next outermost dimension. If that dimension has size 1,
 * the process continues with the next outermost dimension.
 *
 * Other ImageRegionSplitter subclasses could divide an image into
 * more uniform shaped regions instead of slabs.
 *
 * \sa ImageRegionMultidimensionalSplitter
 *
 * \ingroup ITKSystemObjects
 * \ingroup DataProcessing
 * \ingroup ITK-Common
 */

template< unsigned int VImageDimension >
class ITK_EXPORT ImageRegionSplitter:public Object
{
public:
  /** Standard class typedefs. */
  typedef ImageRegionSplitter        Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageRegionSplitter, Object);

  /** Dimension of the image available at compile time. */
  itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

  /** Dimension of the image available at run time. */
  static unsigned int GetImageDimension()
  { return VImageDimension; }

  /** Index typedef support. An index is used to access pixel values. */
  typedef Index< VImageDimension >           IndexType;

  /** Size typedef support. A size is used to define region bounds. */
  typedef Size< VImageDimension >          SizeType;
  typedef typename SizeType::SizeValueType SizeValueType;

  /** Region typedef support.   */
  typedef ImageRegion< VImageDimension > RegionType;

  /** How many pieces can the specifed region be split? A given region
   * cannot always be divided into the requested number of pieces.  For
   * instance, if the numberOfPieces exceeds the number of pixels along
   * a certain dimensions, then some splits will not be possible. This
   * method returns a number less than or equal to the requested number
   * of pieces. */
  virtual unsigned int GetNumberOfSplits(const RegionType & region,
                                         unsigned int requestedNumber);

  /** Get a region definition that represents the ith piece a specified region.
   * The "numberOfPieces" must be equal to what
   * GetNumberOfSplits() returns. */
  virtual RegionType GetSplit(unsigned int i, unsigned int numberOfPieces,
                              const RegionType & region);

protected:
  ImageRegionSplitter() {}
  ~ImageRegionSplitter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  ImageRegionSplitter(const ImageRegionSplitter &); //purposely not implemented
  void operator=(const ImageRegionSplitter &);      //purposely not implemented
};
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_ImageRegionSplitter(_, EXPORT, TypeX, TypeY)                   \
  namespace itk                                                                     \
  {                                                                                 \
  _( 1 ( class EXPORT ImageRegionSplitter< ITK_TEMPLATE_1 TypeX > ) )               \
  namespace Templates                                                               \
  {                                                                                 \
  typedef ImageRegionSplitter< ITK_TEMPLATE_1 TypeX > ImageRegionSplitter##TypeY; \
  }                                                                                 \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkImageRegionSplitter+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkImageRegionSplitter.txx"
#endif

#endif
