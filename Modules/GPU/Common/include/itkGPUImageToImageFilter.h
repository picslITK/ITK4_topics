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
#ifndef __itkGPUImageToImageFilter_h
#define __itkGPUImageToImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkGPUKernelManager.h"

namespace itk
{

/** \class GPUImageToImageFilter
 *
 * \brief class to abstract the behaviour of the GPU filters.
 *
 * FIXME   Won-Ki to write more documentation here...
 *
 * \ingroup ITKGPUCommon
 */
template< class TInputImage, class TOutputImage, class TParentImageFilter = ImageToImageFilter< TInputImage, TOutputImage > >
class ITK_EXPORT GPUImageToImageFilter: public TParentImageFilter
{
public:
  /** Standard class typedefs. */
  typedef GPUImageToImageFilter       Self;
  typedef TParentImageFilter          Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GPUImageToImageFilter, TParentImageFilter);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;

  /** Some convenient typedefs. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

  // macro to set if GPU is used
  itkSetMacro(GPUEnabled, bool);
  itkGetConstMacro(GPUEnabled, bool);
  itkBooleanMacro(GPUEnabled);

  void GenerateData();

protected:
  GPUImageToImageFilter();
  ~GPUImageToImageFilter();

  virtual void PrintSelf(std::ostream & os, Indent indent) const;

  virtual void GPUGenerateData() {};

  // GPU kernel manager
  typename GPUKernelManager::Pointer m_KernelManager;

private:
  GPUImageToImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);        //purposely not implemented

  bool m_GPUEnabled;
};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_GPUImageToImageFilter(_, EXPORT, TypeX, TypeY)                  \
  namespace itk                                                                   \
{                                                                               \
  _( 2 ( class EXPORT GPUImageToImageFilter< ITK_TEMPLATE_2 TypeX > ) )              \
  namespace Templates                                                             \
{                                                                               \
  typedef GPUImageToImageFilter< ITK_TEMPLATE_2 TypeX > GPUImageToImageFilter##TypeY; \
}                                                                               \
}

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkGPUImageToImageFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkGPUImageToImageFilter.hxx"
#endif

#endif
