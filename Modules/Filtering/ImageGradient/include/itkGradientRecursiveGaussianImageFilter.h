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
#ifndef __itkGradientRecursiveGaussianImageFilter_h
#define __itkGradientRecursiveGaussianImageFilter_h

#include "itkRecursiveGaussianImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkImage.h"
#include "itkCovariantVector.h"
#include "itkPixelTraits.h"
#include "itkProgressAccumulator.h"
#include <vector>

namespace itk
{
/** \class GradientRecursiveGaussianImageFilter
 * \brief Computes the gradient of an image by convolution
 *        with the first derivative of a Gaussian.
 *
 * This filter is implemented using the recursive gaussian
 * filters
 *
 *
 * \ingroup GradientFilters
 * \ingroup Singlethreaded
 * \ingroup ITK-ImageGradient
 * \wikiexample{EdgesAndGradients/GradientRecursiveGaussianImageFilter,Compute the gradient of an imabe by convolution with the first derivative of a Gaussian}
 */
// NOTE that the ITK_TYPENAME macro has to be used here in lieu
// of "typename" because VC++ doesn't like the typename keyword
// on the defaults of template parameters
template< typename TInputImage,
          typename TOutputImage = Image< CovariantVector<
                                           ITK_TYPENAME NumericTraits< ITK_TYPENAME TInputImage::PixelType >::RealType,
                                           ::itk::GetImageDimension< TInputImage >::ImageDimension >,
                                         ::itk::GetImageDimension< TInputImage >::ImageDimension > >
class ITK_EXPORT GradientRecursiveGaussianImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef GradientRecursiveGaussianImageFilter            Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Pixel Type of the input image */
  typedef TInputImage                                   InputImageType;
  typedef typename TInputImage::PixelType               PixelType;
  typedef typename NumericTraits< PixelType >::RealType RealType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Define the image type for internal computations
      RealType is usually 'double' in NumericTraits.
      Here we prefer float in order to save memory.  */

  typedef float InternalRealType;
  typedef Image< InternalRealType,
                 itkGetStaticConstMacro(ImageDimension) >   RealImageType;

  /**  Output Image Nth Element Adaptor
   *  This adaptor allows to use conventional scalar
   *  smoothing filters to compute each one of the
   *  components of the gradient image pixels. */
  typedef NthElementImageAdaptor< TOutputImage,
                                  InternalRealType >  OutputImageAdaptorType;

  typedef typename OutputImageAdaptorType::Pointer OutputImageAdaptorPointer;

  /**  Smoothing filter type */
  typedef RecursiveGaussianImageFilter<
    RealImageType,
    RealImageType
    >    GaussianFilterType;

  /**  Derivative filter type, it will be the first in the pipeline  */
  typedef RecursiveGaussianImageFilter<
    InputImageType,
    RealImageType
    >    DerivativeFilterType;

  /**  Pointer to a gaussian filter.  */
  typedef typename GaussianFilterType::Pointer GaussianFilterPointer;

  /**  Pointer to a derivative filter.  */
  typedef typename DerivativeFilterType::Pointer DerivativeFilterPointer;

  /**  Pointer to the Output Image */
  typedef typename TOutputImage::Pointer OutputImagePointer;

  /** Type of the output Image */
  typedef TOutputImage                                       OutputImageType;
  typedef typename          OutputImageType::PixelType       OutputPixelType;
  typedef typename PixelTraits< OutputPixelType >::ValueType OutputComponentType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(GradientRecursiveGaussianImageFilter,
               ImageToImageFilter);

  /** Set Sigma value. Sigma is measured in the units of image spacing.  */
  void SetSigma(RealType sigma);

  /** Define which normalization factor will be used for the Gaussian
   *  \sa  RecursiveGaussianImageFilter::SetNormalizeAcrossScale
   */
  void SetNormalizeAcrossScale(bool normalizeInScaleSpace);
  itkGetConstMacro(NormalizeAcrossScale, bool);

  /** GradientRecursiveGaussianImageFilter needs all of the input to produce an
   * output. Therefore, GradientRecursiveGaussianImageFilter needs to provide
   * an implementation for GenerateInputRequestedRegion in order to inform
   * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion()
  throw( InvalidRequestedRegionError );

  /** The UseImageDirection flag determines whether the gradients are
   * computed with respect to the image grid or with respect to the physical
   * space. When this flag is ON the gradients are computed with respect to
   * the coodinate system of physical space. The difference is whether we take
   * into account the image Direction or not. The flag ON will take into
   * account the image direction and will result in an extra matrix
   * multiplication compared to the amount of computation performed when the
   * flag is OFF.
   * The default value of this flag is On.
   */
  itkSetMacro(UseImageDirection, bool);
  itkGetConstMacro(UseImageDirection, bool);
  itkBooleanMacro(UseImageDirection);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< PixelType > ) );
  itkConceptMacro( OutputHasPixelTraitsCheck,
                   ( Concept::HasPixelTraits< OutputPixelType > ) );
  /** End concept checking */
#endif
protected:
  GradientRecursiveGaussianImageFilter();
  virtual ~GradientRecursiveGaussianImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Generate Data */
  void GenerateData(void);

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(DataObject *output);

private:
  GradientRecursiveGaussianImageFilter(const Self &); //purposely not
                                                      // implemented
  void operator=(const Self &);                       //purposely not

  // implemented

  std::vector< GaussianFilterPointer > m_SmoothingFilters;
  DerivativeFilterPointer              m_DerivativeFilter;
  OutputImageAdaptorPointer            m_ImageAdaptor;

  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale;

  /** Take into account image orientation when computing the Gradient */
  bool m_UseImageDirection;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGradientRecursiveGaussianImageFilter.txx"
#endif

#endif
