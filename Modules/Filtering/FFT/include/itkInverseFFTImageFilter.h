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
#ifndef __itkInverseFFTImageFilter_h
#define __itkInverseFFTImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
/** \class InverseFFTImageFilter
 *
 * \brief Base class for "Inverse" Fast Fourier Transform.
 *
 * This is a base class for the "inverse" or "reverse" Discrete Fourier
 * Transform.  This is an abstract base class: the actual implementation is
 * provided by the best child available on the system when the object is
 * created via the object factory system.
 *
 * This class transforms a complex conjugate symmetric image into its real
 * spatial domain representation.  If the input is not complex conjugate symmetric, the
 * imaginary component is discarded.  The transform of a real input image has
 * complex conjugate symmetry.  That is, values in the second half of the
 * transform are the complex conjugates of values in the first half.  Some
 * implementations, e.g. FFTW, may take advantage of this property and reduce
 * the size of the output in one direction during the forward transform  to
 * N/2+1, where N is the size of the input.  If this occurs, FullMatrix()
 * returns 'false'.  If this was the case, the size of the inverse output image
 * will be larger than the input.
 *
 * \ingroup FourierTransform
 *
 * \sa ForwardFFTImageFilter, InverseFFTImageFilter
 * \ingroup ITKFFT
 */
template< class TInputImage, class TOutputImage=Image< typename TInputImage::PixelType::value_type, TInputImage::ImageDimension> >
class ITK_EXPORT InverseFFTImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >

{
public:
  /** Standard class typedefs. */
  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  typedef InverseFFTImageFilter                                 Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      InputImageType::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(InverseFFTImageFilter, ImageToImageFilter);

  /** Customized object creation methods that support configuration-based
  * selection of FFT implementation.
  *
  * Default implementation is VnlFFT. */
  static Pointer New(void);

  /** The output may be a different size from the input if complex conjugate
   * symmetry is implicit. */
  virtual void GenerateOutputInformation();

  /** This class requires the entire input. */
  virtual void GenerateInputRequestedRegion();

  /** Returns true if the output's size is the same size as the input, i.e.
   * we do not take advantage of complex conjugate symmetry. */
  virtual bool FullMatrix() = 0; // must be implemented in child

  /** Was the original truncated dimension size odd? */
  void SetActualXDimensionIsOdd(bool isodd)
  {
    m_ActualXDimensionIsOdd = isodd;
  }

  void SetActualXDimensionIsOddOn()
  {
    this->SetActualXDimensionIsOdd(true);
  }

  void SetActualXDimensionIsOddOff()
  {
    this->SetActualXDimensionIsOdd(false);
  }

  bool ActualXDimensionIsOdd()
  {
    return m_ActualXDimensionIsOdd;
  }

protected:
  InverseFFTImageFilter():m_ActualXDimensionIsOdd(false) {}
  virtual ~InverseFFTImageFilter(){}

  /** Sets the output requested region to the largest possible output
   * region. */
  void EnlargeOutputRequestedRegion( DataObject *itkNotUsed(output) );

private:
  bool m_ActualXDimensionIsOdd;
  InverseFFTImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                       //purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#ifndef __itkVnlInverseFFTImageFilter_h
#ifndef __itkVnlInverseFFTImageFilter_hxx
#ifndef __itkFFTWInverseFFTImageFilter_h
#ifndef __itkFFTWInverseFFTImageFilter_hxx
#include "itkInverseFFTImageFilter.hxx"
#endif
#endif
#endif
#endif
#endif

#endif
