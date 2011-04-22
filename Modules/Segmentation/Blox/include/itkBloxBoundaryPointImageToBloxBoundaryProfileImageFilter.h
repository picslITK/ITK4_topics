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
#ifndef __itkBloxBoundaryPointImageToBloxBoundaryProfileImageFilter_h
#define __itkBloxBoundaryPointImageToBloxBoundaryProfileImageFilter_h

#include "itkBloxBoundaryProfileImage.h"
#include "itkImageToImageFilter.h"
#include "itkSize.h"

namespace itk
{
/** \class BloxBoundaryPointImageToBloxBoundaryProfileImageFilter
 * \brief Converts a BloxImage of BloxBoundaryPoints to a BloxImage of
 * BloxBoundaryProfiles
 *
 * Samples the BloxBoundaryPointImage to form a BloxBoundaryProfileImage by
 * sampling voxels in an ellipsoidal region, where the center of the
 * ellipoid is the location of each boundary point. Voxels within these
 * regions are splatted onto the major axis of the ellipsoid in bins
 * to form a profile of average intensity traversing blurred boundaries.
 * Using curve fitting techniques, a cumulative Gaussian is fit to this
 * intensity profile to yield estimates of actual boundary location,
 * intensities on both sides of the boundary, and blurred boundary width.
 *
 * References:
 * Robert J. Tamburo, George D. Stetten: Gradient-Oriented Profiles
 * for Boundary Parameterization and Their Application to Core Atoms
 * Towards Shape Analysis. International Journal of Image and
 * Graphics 1(4): 659-680 (2001)
 *
 * Robert J.Tamburo and George D.Stetten,M.D.,Ph.D.: Gradient-Oriented
 * Profiles for Unsupervised Boundary Classification. Proceedings of
 * the 29th Applied Imagery Pattern Recognition Workshop. October 2000:
 * Washington, D.C.
 *
 * \ingroup ImageEnhancement
 * \ingroup ITK-Blox
 */
template< typename TSourceImage >
class ITK_EXPORT BloxBoundaryPointImageToBloxBoundaryProfileImageFilter:
  public ImageToImageFilter< TSourceImage,
                             BloxBoundaryProfileImage< ::itk::GetImageDimension< TSourceImage >::ImageDimension > >
{
public:
  /** Number of dimensions */
  itkStaticConstMacro(NDimensions, unsigned int, TSourceImage::ImageDimension);

  /** Standard class typedefs */
  typedef BloxBoundaryPointImageToBloxBoundaryProfileImageFilter Self;
  typedef ImageToImageFilter< TSourceImage,
                              BloxBoundaryProfileImage< itkGetStaticConstMacro(NDimensions) > > Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(BloxBoundaryPointImageToBloxBoundaryProfileImageFilter, ImageToImageFilter);

  /** Typedef for boundary point image */
  typedef BloxBoundaryPointImage< itkGetStaticConstMacro(NDimensions) >
  BoundaryPointImageType;
  typedef typename BoundaryPointImageType::Pointer      BoundaryPointImagePointer;
  typedef typename BoundaryPointImageType::RegionType   BoundaryPointImageRegionType;
  typedef typename BoundaryPointImageType::PixelType    BoundaryPointImagePixelType;
  typedef typename BoundaryPointImageType::ConstPointer BoundaryPointImageConstPointer;

  /** Typedef for blurred source image */
  typedef TSourceImage                           SourceImageType;
  typedef typename SourceImageType::Pointer      SourceImagePointer;
  typedef typename SourceImageType::RegionType   SourceImageRegionType;
  typedef typename SourceImageType::PixelType    SourceImagePixelType;
  typedef typename SourceImageType::ConstPointer SourceImageConstPointer;

  /** Typedef for profile image */
  typedef BloxBoundaryProfileImage< itkGetStaticConstMacro(NDimensions) >
  OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;

  /** Image index typedef */
  typedef typename BloxBoundaryProfileImage< itkGetStaticConstMacro(NDimensions) >::IndexType IndexType;

  /** Image pixel value typedef */
  typedef typename BloxBoundaryProfileImage< itkGetStaticConstMacro(NDimensions) >::PixelType PixelType;

  /** The type of vector used to convert between physical and blox space */
  typedef Point< double, itkGetStaticConstMacro(NDimensions) > PositionType;

  /** Vector typedef */
  typedef typename PositionType::VectorType VectorType;

  /** Set the blurred original image */
  void SetInput1(const SourceImageType *image1);

  /** Set the boundary point image */
  void SetInput2(const BoundaryPointImageType *image2);

  /** Find maximum in accumulator */
  double FindAccumulatorMaximum();

  /** Find minimum in accumulator */
  double FindAccumulatorMinimum();

  /** Find boundary profiles from input images and store them */
  void FindBoundaryProfiles();

  /** Add weighted pixel value to appropriate bin number in splat accumulator
    and normalizer */
  bool AddSplatToAccumulatorAndNormalizer(int binNumber, double weight, double sourcePixelValue);

  /** Normalize the splat accumulator by the normalizer */
  void NormalizeSplatAccumulator();

  /** Fit the boundary profile to a cumulative Gaussian */
  int FitProfile();

  /** Parameters required to find boundary profiles */
  void Initialize(double setUniqueAxis, double setSymmetricAxes, unsigned int numberOfBins,
                  unsigned int splatMethod, unsigned int spaceDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( SourceConvertibleToDoubleCheck,
                   ( Concept::Convertible< typename TSourceImage::PixelType, double > ) );
  /** End concept checking */
#endif
protected:
  BloxBoundaryPointImageToBloxBoundaryProfileImageFilter();
  ~BloxBoundaryPointImageToBloxBoundaryProfileImageFilter();

  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Method for forming the BloxBoundaryProfileImage */
  void GenerateData();

private:
  //purposely not implemented
  BloxBoundaryPointImageToBloxBoundaryProfileImageFilter(const Self &);
  void operator=(const Self &);

  /** Length of major axis of ellipsoid */
  double m_UniqueAxis;

  /** Lengths of minor axes of ellipsoid */
  double m_SymmetricAxes;

  /** Number of bins in splat accumulator and normalizer */
  unsigned int m_NumberOfBins;

  /** Type of method to splat. 0: Gaussian 1: Triangle  */
  unsigned int m_SplatMethod;

  /** Count of number of boundary profiles found */
  unsigned long int m_NumBoundaryProfiles;

  /** Weight pixel values */
  double *m_Accumulator;

  /** Count of pixels added to accumulator */
  double *m_Normalizer;

  /** Normalized accumulator */
  double *m_NormalizedAccumulator;

  /** Final parameters delivered by FitProfile() */
  double *m_FinalParameters;

  /** Number of parameters in cost function */
  unsigned int m_SpaceDimension;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBloxBoundaryPointImageToBloxBoundaryProfileImageFilter.txx"
#endif

#endif
