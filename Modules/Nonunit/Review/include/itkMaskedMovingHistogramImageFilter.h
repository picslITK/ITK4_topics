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
#ifndef __itkMaskedMovingHistogramImageFilter_h
#define __itkMaskedMovingHistogramImageFilter_h

#include "itkMovingHistogramImageFilterBase.h"
#include <list>
#include <map>
#include <set>

namespace itk
{
/**
 * \class MaskedMovingHistogramImageFilter
 *
 * \author Richard Beare
 * \author Gaetan Lehmann
 * \ingroup ITK-Review
 */

template< class TInputImage, class TMaskImage, class TOutputImage, class TKernel, class THistogram >
class ITK_EXPORT MaskedMovingHistogramImageFilter:
  public MovingHistogramImageFilterBase< TInputImage, TOutputImage, TKernel >
{
public:
  /** Standard class typedefs. */
  typedef MaskedMovingHistogramImageFilter                                     Self;
  typedef MovingHistogramImageFilterBase< TInputImage, TOutputImage, TKernel > Superclass;
  typedef SmartPointer< Self >                                                 Pointer;
  typedef SmartPointer< const Self >                                           ConstPointer;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(MaskedMovingHistogramImageFilter, MovingHistogramImageFilter);

  /** Image related typedefs. */
  typedef TInputImage                                InputImageType;
  typedef TOutputImage                               OutputImageType;
  typedef TMaskImage                                 MaskImageType;
  typedef typename TInputImage::RegionType           RegionType;
  typedef typename TInputImage::SizeType             SizeType;
  typedef typename TInputImage::IndexType            IndexType;
  typedef typename TInputImage::PixelType            PixelType;
  typedef typename TInputImage::OffsetType           OffsetType;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename TOutputImage::PixelType           OutputPixelType;
  typedef typename TInputImage::PixelType            InputPixelType;
  typedef typename MaskImageType::PixelType          MaskPixelType;
  typedef THistogram                                 HistogramType;

  /** Set the marker image */
  void SetMaskImage(MaskImageType *input)
  {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast< TMaskImage * >( input ) );
  }

  /** Get the marker image */
  MaskImageType * GetMaskImage()
  {
    return static_cast< MaskImageType * >( const_cast< DataObject * >( this->ProcessObject::GetInput(1) ) );
  }

  /** Set the input image */
  void SetInput1(InputImageType *input)
  {
    this->SetInput(input);
  }

  /** Set the marker image */
  void SetInput2(MaskImageType *input)
  {
    this->SetMaskImage(input);
  }

  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Kernel typedef. */
  typedef TKernel KernelType;

  /** Kernel (structuring element) iterator. */
  typedef typename KernelType::ConstIterator KernelIteratorType;

  /** n-dimensional Kernel radius. */
  typedef typename KernelType::SizeType RadiusType;

  typedef typename std::list< OffsetType > OffsetListType;

  typedef typename std::map< OffsetType, OffsetListType,
                             typename Functor::OffsetLexicographicCompare< itkGetStaticConstMacro(ImageDimension) > >
  OffsetMapType;

  /** Get the modified mask image */
  MaskImageType * GetOutputMask();

  void AllocateOutputs();

  DataObject::Pointer MakeOutput(unsigned int idx);

  itkSetMacro(FillValue, OutputPixelType);
  itkGetConstMacro(FillValue, OutputPixelType);

  itkSetMacro(MaskValue, MaskPixelType);
  itkGetConstMacro(MaskValue, MaskPixelType);

  itkSetMacro(BackgroundMaskValue, MaskPixelType);
  itkGetConstMacro(BackgroundMaskValue, MaskPixelType);

  void SetGenerateOutputMask(bool);

  itkGetConstMacro(GenerateOutputMask, bool);
//   itkBooleanMacro(GenerateOutputMask);

  /** ConfigurewHistogram can be used to configure the histogram. The default version just do nothing. */
  virtual void ConfigureHistogram(THistogram &) {}

protected:
  MaskedMovingHistogramImageFilter();
  ~MaskedMovingHistogramImageFilter() {}

  /** Multi-thread version GenerateData. */
  void  ThreadedGenerateData(const OutputImageRegionType &
                             outputRegionForThread,
                             int threadId);

  void PrintSelf(std::ostream & os, Indent indent) const;

  void pushHistogram(HistogramType & histogram,
                     const OffsetListType *addedList,
                     const OffsetListType *removedList,
                     const RegionType & inputRegion,
                     const RegionType & kernRegion,
                     const InputImageType *inputImage,
                     const MaskImageType *maskImage,
                     const IndexType currentIdx);

private:
  MaskedMovingHistogramImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                   //purposely not implemented

  bool m_GenerateOutputMask;

  OutputPixelType m_FillValue;

  MaskPixelType m_MaskValue;

  MaskPixelType m_BackgroundMaskValue;
}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMaskedMovingHistogramImageFilter.txx"
#endif

#endif
