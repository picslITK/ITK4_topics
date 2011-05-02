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
#ifndef __itkLabelMapOverlayImageFilter_h
#define __itkLabelMapOverlayImageFilter_h

#include "itkLabelMapFilter.h"
#include "itkBarrier.h"
#include "itkLabelOverlayFunctor.h"
#include "itkRGBPixel.h"

namespace itk {

/** \class LabelMapOverlayImageFilter
* \brief Apply a colormap to a label map and put it on top of the
 *  input image
 *
 * Apply a colormap to a label map and put it on top of the input image.
 * The set of colors is a good selection of distinct colors. The opacity of
 * the label map can be defined by the user.
 * A background label produce a gray pixel with the same intensity
 * than the input one.
 *
 * \author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.
 *
 * This implementation was taken from the Insight Journal paper:
 * http://hdl.handle.net/1926/584  or
 * http://www.insight-journal.org/browse/publication/176
 *
 * \sa LabelMapToBinaryImageFilter, LabelMapToLabelImageFilter
 * \ingroup ImageEnhancement  MathematicalMorphologyImageFilters
 * \ingroup ITK-Review
*/
template<class TInputImage, class TFeatureImage, class TOutputImage=Image< RGBPixel< typename TFeatureImage::PixelType >, TFeatureImage::ImageDimension > >
class ITK_EXPORT LabelMapOverlayImageFilter :
    public LabelMapFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef LabelMapOverlayImageFilter                Self;
  typedef LabelMapFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::RegionType      InputImageRegionType;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename InputImageType::LabelObjectType LabelObjectType;
  typedef typename LabelObjectType::LabelType      LabelType;
  typedef typename LabelObjectType::LengthType     LengthType;

  typedef TFeatureImage                             FeatureImageType;
  typedef typename FeatureImageType::Pointer        FeatureImagePointer;
  typedef typename FeatureImageType::ConstPointer   FeatureImageConstPointer;
  typedef typename FeatureImageType::RegionType     FeatureImageRegionType;
  typedef typename FeatureImageType::PixelType      FeatureImagePixelType;

  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::RegionType     OutputImageRegionType;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename OutputImageType::IndexType      IndexType;
  typedef typename OutputImageType::SizeType       SizeType;
  typedef typename OutputImageType::RegionType     RegionType;

  typedef typename Functor::LabelOverlayFunctor<FeatureImagePixelType, InputImagePixelType, OutputImagePixelType> FunctorType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(LabelMapOverlayImageFilter,
               ImageToImageFilter);

   /** Set the feature image */
  void SetFeatureImage(const TFeatureImage *input)
    {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast<TFeatureImage *>(input) );
    }

  /** Get the feature image */
  const FeatureImageType * GetFeatureImage()
    {
    return static_cast<FeatureImageType*>((this->ProcessObject::GetInput(1)));
    }

   /** Set the input image */
  void SetInput1(const TInputImage *input)
    {
    this->SetInput( input );
    }

  /** Set the feature image */
  void SetInput2(const TFeatureImage *input)
    {
    this->SetFeatureImage( input );
    }

  /** Set/Get the opacity of the colored label image. The value must be
   * between 0 and 1
   */
  itkSetMacro( Opacity, double );
  itkGetConstReferenceMacro( Opacity, double );

protected:
  LabelMapOverlayImageFilter();
  ~LabelMapOverlayImageFilter() {};

  /** LabelMapOverlayImageFilter needs the entire input be
   * available. Thus, it needs to provide an implementation of
   * GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion();

  /** LabelMapOverlayImageFilter will produce the entire output. */
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));

  virtual void BeforeThreadedGenerateData();

  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId );

  virtual void ThreadedProcessLabelObject( LabelObjectType * labelObject );

  virtual void GenerateOutputInformation();

  void PrintSelf(std::ostream& os, Indent indent) const;


private:
  LabelMapOverlayImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  double                    m_Opacity;
  typename Barrier::Pointer m_Barrier;

}; // end of class

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelMapOverlayImageFilter.txx"
#endif

#endif
