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

#ifndef __itkConvolutionImageFilter_hxx
#define __itkConvolutionImageFilter_hxx

#include "itkConvolutionImageFilter.h"

#include "itkConstantPadImageFilter.h"
#include "itkCropImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageBase.h"
#include "itkImageKernelOperator.h"
#include "itkNeighborhoodOperatorImageFilter.h"
#include "itkNormalizeToConstantImageFilter.h"

namespace itk
{
template< class TInputImage, class TKernelImage, class TOutputImage >
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::ConvolutionImageFilter()
{
  m_Normalize = false;
  m_BoundaryCondition = &m_DefaultBoundaryCondition;
  m_OutputRegionMode = Self::SAME;
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::GenerateData()
{
  // Allocate the output
  this->AllocateOutputs();

  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter( this );

  // Build a mini-pipeline that involves a
  // NeighborhoodOperatorImageFilter to compute the convolution, a
  // normalization filter for the kernel, and a pad filter for making
  // the kernel an odd size.
  if ( m_Normalize )
    {
    typedef typename NumericTraits< InputPixelType >::RealType RealPixelType;
    typedef Image< RealPixelType, ImageDimension >             RealImageType;

    typedef NormalizeToConstantImageFilter< KernelImageType, RealImageType > NormalizeFilterType;
    typename NormalizeFilterType::Pointer normalizeFilter = NormalizeFilterType::New();
    normalizeFilter->SetConstant( NumericTraits< RealPixelType >::One );
    normalizeFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    normalizeFilter->SetInput( this->GetKernelImage() );
    normalizeFilter->ReleaseDataFlagOn();
    progress->RegisterInternalFilter( normalizeFilter, 0.1f );
    normalizeFilter->UpdateLargestPossibleRegion();

    this->ComputeConvolution( normalizeFilter->GetOutput(), progress );
    }
  else
    {
    this->ComputeConvolution( this->GetKernelImage(), progress );
    }
}

template< class TInputImage, class TKernelImage, class TOutputImage >
template< class TImage >
void
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::ComputeConvolution( const TImage * kernelImage,
                      ProgressAccumulator * progress )
{
  typedef typename TImage::PixelType KernelImagePixelType;
  typedef ImageKernelOperator< KernelImagePixelType, ImageDimension > KernelOperatorType;
  KernelOperatorType kernelOperator;

  bool kernelNeedsPadding = this->GetKernelNeedsPadding();

  float optionalFilterWeights = 0.0f;
  if ( m_Normalize )
    {
    optionalFilterWeights += 0.1f;
    }
  if ( this->GetKernelNeedsPadding() )
    {
    optionalFilterWeights += 0.1f;
    }
  if ( m_OutputRegionMode == Self::VALID )
    {
    optionalFilterWeights += 0.1f;
    }

  // Flip the kernel
  typedef FlipImageFilter< TImage > FlipperType;
  typename FlipperType::Pointer flipper = FlipperType::New();
  typename FlipperType::FlipAxesArrayType axesArray;
  axesArray.Fill( true );
  flipper->SetFlipAxes( axesArray );
  flipper->SetInput( kernelImage );

  if ( kernelNeedsPadding )
    {
    // Pad the kernel if necessary to an odd size in each dimension.
    typedef ConstantPadImageFilter< TImage, TImage > PadImageFilterType;
    typename PadImageFilterType::Pointer kernelPadImageFilter = PadImageFilterType::New();
    kernelPadImageFilter->SetConstant( NumericTraits< KernelImagePixelType >::ZeroValue() );
    kernelPadImageFilter->SetPadLowerBound( this->GetKernelPadSize() );
    kernelPadImageFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    kernelPadImageFilter->ReleaseDataFlagOn();
    kernelPadImageFilter->SetInput( flipper->GetOutput() );
    progress->RegisterInternalFilter( kernelPadImageFilter, 0.1f );
    kernelPadImageFilter->UpdateLargestPossibleRegion();

    kernelOperator.SetImageKernel( kernelPadImageFilter->GetOutput() );
    }
  else
    {
    flipper->UpdateLargestPossibleRegion();
    kernelOperator.SetImageKernel( flipper->GetOutput() );
    }

  KernelSizeType radius = this->GetKernelRadius( kernelImage );
  kernelOperator.CreateToRadius( radius );

  typename InputImageType::Pointer localInput = InputImageType::New();
  localInput->Graft( this->GetInput() );

  // The NeighborhoodOperatorImageFilter does all the work.
  typedef NeighborhoodOperatorImageFilter< InputImageType, OutputImageType, KernelImagePixelType >
    ConvolutionFilterType;
  typename ConvolutionFilterType::Pointer convolutionFilter = ConvolutionFilterType::New();
  convolutionFilter->SetOperator( kernelOperator );
  convolutionFilter->OverrideBoundaryCondition( m_BoundaryCondition );
  convolutionFilter->SetInput( localInput );
  convolutionFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
  convolutionFilter->ReleaseDataFlagOn();
  progress->RegisterInternalFilter( convolutionFilter, 1.0f - optionalFilterWeights );

  if ( m_OutputRegionMode == Self::SAME )
    {
    // Graft the output of the convolution filter onto this filter's
    // output.
    convolutionFilter->GraftOutput( this->GetOutput() );
    convolutionFilter->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );
    convolutionFilter->Update();
    this->GraftOutput( convolutionFilter->GetOutput() );
    }
  else // m_OutputRegionMode == Self::VALID
    {
    typedef CropImageFilter< OutputImageType, OutputImageType > CropFilterType;
    typedef typename CropFilterType::Pointer                    CropFilterPointer;
    typedef typename CropFilterType::SizeType                   CropSizeType;

    // Set up the crop sizes.
    CropSizeType upperCropSize( radius );
    CropSizeType lowerCropSize( radius );

    // For the lower crop, the crop size can be reduced by 1 in a
    // dimension when the kernel size is odd in that dimension.
    lowerCropSize -= this->GetKernelPadSize();

    // Set up the crop filter.
    CropFilterPointer cropFilter = CropFilterType::New();
    cropFilter->SetLowerBoundaryCropSize( lowerCropSize );
    cropFilter->SetUpperBoundaryCropSize( upperCropSize );
    cropFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    cropFilter->ReleaseDataFlagOn();
    progress->RegisterInternalFilter( cropFilter, 0.1f );
    cropFilter->SetInput( convolutionFilter->GetOutput() );

    // Graft the minipipeline output to this filter.
    cropFilter->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );
    cropFilter->GraftOutput( this->GetOutput() );
    cropFilter->Update();

    // Graft the output of the crop filter back onto this
    // filter's output.
    this->GraftOutput( cropFilter->GetOutput() );
    }
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::GenerateOutputInformation()
{
  // Call the superclass' implementation of this method. Default
  // behavior corresponds to SAME output region mode.
  Superclass::GenerateOutputInformation();

  if ( m_OutputRegionMode == Self::VALID )
    {
    OutputRegionType validRegion = this->GetValidRegion();

    typename OutputImageType::Pointer outputPtr  = this->GetOutput();
    outputPtr->SetLargestPossibleRegion( validRegion );
    }
}

template< class TInputImage, class TKernelImage, class TOutputImage >
bool
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::GetKernelNeedsPadding() const
{
  const KernelImageType *kernel = this->GetKernelImage();
  InputRegionType kernelRegion = kernel->GetLargestPossibleRegion();
  InputSizeType kernelSize = kernelRegion.GetSize();

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( kernelSize[i] % 2 == 0 ) // Check if dimension is even
      {
      return true;
      }
    }

  return false;
}

template< class TInputImage, class TKernelImage, class TOutputImage >
typename ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >::KernelSizeType
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::GetKernelPadSize() const
{
  const KernelImageType *kernel = this->GetKernelImage();
  KernelRegionType kernelRegion = kernel->GetLargestPossibleRegion();
  KernelSizeType kernelSize = kernelRegion.GetSize();
  KernelSizeType padSize;

  for ( unsigned int i = 0; i < ImageDimension; i++)
    {
    // Pad by 1 if the size fo the image in this dimension is even.
    padSize[i] = 1 - (kernelSize[i] % 2);
    }

  return padSize;
}

template< class TInputImage, class TKernelImage, class TOutputImage >
template< class TImage >
typename ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >::KernelSizeType
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::GetKernelRadius(const TImage *kernelImage) const
{
  // Compute the kernel radius.
  KernelSizeType radius;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    radius[i] = kernelImage->GetLargestPossibleRegion().GetSize()[i] / 2;
    }

  return radius;
}

template< class TInputImage, class TKernelImage, class TOutputImage >
typename ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >::OutputRegionType
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::GetValidRegion() const
{
  typename InputImageType::ConstPointer inputPtr  = this->GetInput();

  InputRegionType inputLargestPossibleRegion = inputPtr->GetLargestPossibleRegion();

  OutputIndexType validIndex = inputLargestPossibleRegion.GetIndex();
  OutputSizeType validSize = inputLargestPossibleRegion.GetSize();

  // Shrink the output largest possible region by the kernel radius.
  KernelSizeType kernelSize = this->GetKernelImage()->GetLargestPossibleRegion().GetSize();
  KernelSizeType radius = this->GetKernelRadius( this->GetKernelImage() );
  for ( unsigned int i = 0; i < ImageDimension; ++i )
    {
    if ( validSize[i] < 2*radius[i] )
      {
      validIndex[i] = 0;
      validSize[i] = 0;
      }
    else
      {
      validIndex[i] = validIndex[i] + static_cast< typename OutputIndexType::IndexValueType >
        ( radius[i] );
      validSize[i] = validSize[i] - 2*radius[i];

      // If the kernel is has an even size in one dimension, then we
      // need to expand the image size by one and subtract one from
      // the index in that dimension to account for the zero-padding
      // on the low-index side of the image.
      if ( kernelSize[i] % 2 == 0 )
        {
        validIndex[i] -= 1;
        validSize[i] += 1;
        }
      }
    }

  OutputRegionType validRegion( validIndex, validSize );

  return validRegion;
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::SetOutputRegionModeToSame()
{
  this->SetOutputRegionMode( Self::SAME );
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::SetOutputRegionModeToValid()
{
  this->SetOutputRegionMode( Self::VALID );
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  typedef ImageBase< ImageDimension > ImageBaseType;

  // Pad the input image with the radius of the kernel.
  if ( this->GetInput() )
    {
    InputRegionType inputRegion = this->GetOutput()->GetRequestedRegion();

    // Pad the output request region by the kernel radius.
    KernelSizeType radius = this->GetKernelRadius( this->GetKernelImage() );
    inputRegion.PadByRadius( radius );

    // Crop the output request region to fit within the largest
    // possible region.
    typename InputImageType::Pointer inputPtr =
      const_cast< InputImageType * >( this->GetInput() );
    bool cropped = inputRegion.Crop( inputPtr->GetLargestPossibleRegion() );
    if ( !cropped )
      {
      InvalidRequestedRegionError e( __FILE__, __LINE__ );
      e.SetLocation( ITK_LOCATION );
      e.SetDescription( "Requested region is (at least partially) outside the largest possible region." );
      e.SetDataObject( inputPtr );
      throw e;
      }

    // Input is an image, cast away the constness so we can set
    // the requested region.
    inputPtr->SetRequestedRegion( inputRegion );
    }

  // Request the largest possible region for the kernel image.
  if ( this->GetKernelImage() )
    {
    // Input kernel is an image, cast away the constness so we can set
    // the requested region.
    typename KernelImageType::Pointer kernelPtr =
      const_cast< KernelImageType * >( this->GetKernelImage() );
    kernelPtr->SetRequestedRegionToLargestPossibleRegion();
    }
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
ConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Normalize: "  << m_Normalize << std::endl;
  os << indent << "OutputRegionMode: ";
  switch ( m_OutputRegionMode )
    {
    case SAME:
      os << "SAME";
      break;

    case VALID:
      os << "VALID";
      break;

    default:
      os << "unknown";
      break;
    }
  os << std::endl;
}
}
#endif
