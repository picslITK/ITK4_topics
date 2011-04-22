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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef __itkInPlaceImageFilter_h
#define __itkInPlaceImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
/** \class InPlaceImageFilter
 * \brief Base class for filters that take an image as input and overwrite that image as the output
 *
 * InPlaceImageFilter is the base class for all process objects whose
 * output image data is constructed by overwriting the input image
 * data. In other words, the output bulk data is the same block of
 * memory as the input bulk data.  This filter provides the mechanisms
 * for in place image processing while maintaining general pipeline
 * mechanics. InPlaceImageFilters use less memory than standard
 * ImageToImageFilters because the input buffer is reused as the
 * output buffer.  However, this benefit does not come without a cost.
 * Since the filter overwrites its input, the ownership of the bulk
 * data is transitioned from the input data object to the output data
 * object.  When a data object has multiple consumers with one
 * of the consumers being an in place filter, the in place filter
 * effectively destroys the bulk data for the data object. Upstream
 * filters will then have to re-execute to regenerate the data object's
 * bulk data for the remaining consumers.
 *
 * Since an InPlaceImageFilter reuses the input bulk data memory for the
 * output bulk data memory, the input image type must match the output
 * image type.  If the input and output image types are not identical,
 * the filter reverts to a traditional ImageToImageFilter behaviour
 * where an output image is allocated.  In place operation can also be
 * controlled (when the input and output image type match) via the
 * methods InPlaceOn() and InPlaceOff().
 *
 * Subclasses of InPlaceImageFilter must take extra care in how they
 * manage memory using (and perhaps overriding) the implementations of
 * ReleaseInputs() and AllocateOutputs() provided here.
 *
 * \ingroup ImageFilters
 * \ingroup ITK-Common
 */
template< class TInputImage, class TOutputImage = TInputImage >
class ITK_EXPORT InPlaceImageFilter:public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef InPlaceImageFilter                              Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(InPlaceImageFilter, ImageToImageFilter);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageType       OutputImageType;
  typedef typename Superclass::OutputImagePointer    OutputImagePointer;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;

  /** Some convenient typedefs. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** In place operation can be turned on and off. Asking for
   * in-place operation, i.e. calling SetInplace(true) or InplaceOn(),
   * will be effective only if CanRunInPlace also returns true.
   * By default CanRunInPlace checks whether the input and output
   * image type match. */
  itkSetMacro(InPlace, bool);
  itkGetConstMacro(InPlace, bool);
  itkBooleanMacro(InPlace);

  /** Can the filter run in place? To do so, the filter's first input
   * and output must have the same dimension and pixel type. This
   * method can be used in conjunction with the InPlace ivar to
   * determine whether a particular use of the filter is really
   * running in place. Some filters may be able to optimize their
   * operation if the InPlace is true and CanRunInPlace is true.
   * CanRunInPlace may also be overridded by InPlaceImageFilter
   * subclasses to fine tune its behavior. */
  virtual bool CanRunInPlace() const;

protected:
  InPlaceImageFilter();
  ~InPlaceImageFilter();

  virtual void PrintSelf(std::ostream & os, Indent indent) const;

  /** The GenerateData method normally allocates the buffers for all
   * of the outputs of a filter. Since InPlaceImageFilter's can use an
   * overwritten version of the input for its output, the output
   * buffer should not be allocated. When possible, we graft the input
   * to the filter to the output.  If an InPlaceFilter has multiple
   * outputs, then it would need to override this method to graft one
   * of its outputs and allocate the remaining. If a filter is
   * threaded (i.e. it provides an implementation of
   * ThreadedGenerateData()), this method is called automatically. If
   * an InPlaceFilter is not threaded (i.e. it provides an
   * implementation of GenerateData()), then this method (or
   * equivalent) must be called in GenerateData(). */
  virtual void AllocateOutputs();

  /** InPlaceImageFilter may transfer ownership of the input bulk data
   * to the output object.  Once the output object owns the bulk data
   * (done in AllocateOutputs()), the input object must release its
   * hold on the bulk data.  ProcessObject::ReleaseInputs() only
   * releases the input bulk data when the user has set the
   * ReleaseDataFlag.  InPlaceImageFilter::ReleaseInputs() also
   * releases the input that it has overwritten.
   *
   * \sa ProcessObject::ReleaseInputs() */
  virtual void ReleaseInputs();

private:
  InPlaceImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);     //purposely not implemented

  bool m_InPlace;
  /** borrowed from type_traits */
  struct true_type { };
  struct false_type{ };
  template<typename, typename>
    struct is_same
      : public false_type { };

  template<typename _Tp>
    struct is_same<_Tp, _Tp>
      : public true_type { };

  static bool CanRunInPlace(const true_type &) { return true; }
  static bool CanRunInPlace(const false_type &) { return false; }

};
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_InPlaceImageFilter(_, EXPORT, TypeX, TypeY)                  \
  namespace itk                                                                   \
  {                                                                               \
  _( 2 ( class EXPORT InPlaceImageFilter< ITK_TEMPLATE_2 TypeX > ) )              \
  namespace Templates                                                             \
  {                                                                               \
  typedef InPlaceImageFilter< ITK_TEMPLATE_2 TypeX > InPlaceImageFilter##TypeY; \
  }                                                                               \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkInPlaceImageFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkInPlaceImageFilter.txx"
#endif

#endif
