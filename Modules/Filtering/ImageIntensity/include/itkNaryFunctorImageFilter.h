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
#ifndef __itkNaryFunctorImageFilter_h
#define __itkNaryFunctorImageFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkImageIterator.h"
#include "itkArray.h"

namespace itk
{
/** \class NaryFunctorImageFilter
 * \brief Implements pixel-wise generic operation of Nth similar images.
 *
 * This class is parameterized over the types of the input images
 * and the type of the output image.  It is also parameterized by the
 * operation to be applied.  A Functor style is used to represent the
 * function.
 *
 * All the input images are of the same type.
 *
 * \ingroup IntensityImageFilters   Multithreaded
 * \ingroup ITK-ImageIntensity
 */

template< class TInputImage, class TOutputImage, class TFunction >
class ITK_EXPORT NaryFunctorImageFilter:
  public InPlaceImageFilter< TInputImage, TOutputImage >

{
public:
  /** Standard class typedefs. */
  typedef NaryFunctorImageFilter                          Self;
  typedef InPlaceImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NaryFunctorImageFilter, InPlaceImageFilter);

  /** Some typedefs. */
  typedef TFunction                            FunctorType;
  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename InputImageType::PixelType   InputImagePixelType;
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;
  typedef std::vector< InputImagePixelType >   NaryArrayType;

  /** Get the functor object.  The functor is returned by reference.
   * (Functors do not have to derive from itk::LightObject, so they do
   * not necessarily have a reference count. So we cannot return a
   * SmartPointer). */
  FunctorType & GetFunctor() { return m_Functor; }

  /** Set the functor object.  This replaces the current Functor with a
   * copy of the specified Functor. This allows the user to specify a
   * functor that has ivars set differently than the default functor.
   * This method requires an operator!=() be defined on the functor
   * (or the compiler's default implementation of operator!=() being
   * appropriate). */
  void SetFunctor(FunctorType & functor)
  {
    if ( m_Functor != functor )
      {
      m_Functor = functor;
      this->Modified();
      }
  }

  /** ImageDimension constants */
  itkStaticConstMacro(
    InputImageDimension, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(
    OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( SameDimensionCheck,
                   ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
  itkConceptMacro( OutputHasZeroCheck,
                   ( Concept::HasZero< OutputImagePixelType > ) );
  /** End concept checking */
#endif
protected:
  NaryFunctorImageFilter();
  virtual ~NaryFunctorImageFilter() {}

  /** NaryFunctorImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData() routine
   * which is called for each processing thread. The output image data is
   * allocated automatically by the superclass prior to calling
   * ThreadedGenerateData().  ThreadedGenerateData can only write to the
   * portion of the output image specified by the parameter
   * "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            int threadId);

private:
  NaryFunctorImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);         //purposely not implemented

  FunctorType m_Functor;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNaryFunctorImageFilter.txx"
#endif

#endif
