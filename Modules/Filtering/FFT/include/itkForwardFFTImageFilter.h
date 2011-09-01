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
#ifndef __itkForwardFFTImageFilter_h
#define __itkForwardFFTImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
/** \class ForwardFFTImageFilter
 *
 * \brief Base class for "Forward" Fast Fourier Transform.
 *
 * This is a base class for the "forward" or "direct" discrete Fourier
 * Transform.  This is an abstract base class: the actual implementation is
 * provided by the best child class available on the system when the object is
 * created via the object factory system.
 *
 * This class transforms a real input image into its complex Fourier Transform.
 * The transform of a real input image has complex conjugate symmetry.  That is,
 * values in the second half of the transform are the complex conjugates of
 * values in the first half.  Some implementations, e.g. FFTW, may take
 * advantage of this property and reduce the size of the output in one direction
 * to N/2+1, where N is the size of the input.  If this occurs, FullMatrix()
 * returns 'false'.
 *
 * \ingroup FourierTransform
 *
 * \sa InverseFFTImageFilter, FFTComplexToComplexImageFilter
 * \ingroup ITKFFT
 */
template< class TInputImage, class TOutputImage=Image< std::complex<typename TInputImage::PixelType>, TInputImage::ImageDimension> >
class ITK_EXPORT ForwardFFTImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename InputImageType::IndexType   InputIndexType;
  typedef typename InputImageType::SizeType    InputSizeType;
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::PixelType  OutputPixelType;
  typedef typename OutputImageType::IndexType  OutputIndexType;
  typedef typename OutputIndexType::SizeType   OutputSizeType;

  typedef ForwardFFTImageFilter                                 Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ForwardFFTImageFilter, ImageToImageFilter);

  /** Customized object creation methods that support configuration-based
    * selection of FFT implementation.
    *
    * Default implementation is VnlFFT. */
  static Pointer New(void);

protected:
  ForwardFFTImageFilter() {}
  virtual ~ForwardFFTImageFilter() {}

  /** The output may be a different size from the input if complex conjugate
   * symmetry is implicit. */
  virtual void GenerateOutputInformation();

  /** This class requires the entire input. */
  virtual void GenerateInputRequestedRegion();

  /** This class produces the entire output. */
  virtual void EnlargeOutputRequestedRegion(DataObject *output);

  /** Returns true if the outputs size is the same size as the input, i.e.
   * we do not take advantage of complex conjugate symmetry. */
  virtual bool FullMatrix() = 0; // must be implemented in child

private:
  ForwardFFTImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                       //purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#ifndef __itkVnlForwardFFTImageFilter_h
#ifndef __itkVnlForwardFFTImageFilter_hxx
#ifndef __itkFFTWForwardFFTImageFilter_h
#ifndef __itkFFTWForwardFFTImageFilter_hxx
#include "itkForwardFFTImageFilter.hxx"
#endif
#endif
#endif
#endif
#endif

#endif
