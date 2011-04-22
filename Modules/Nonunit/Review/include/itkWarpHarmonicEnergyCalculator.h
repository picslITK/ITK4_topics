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
#ifndef __itkWarpHarmonicEnergyCalculator_h
#define __itkWarpHarmonicEnergyCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkVector.h"

namespace itk
{
/** \class WarpHarmonicEnergyCalculator
 *
 * \brief Compute the harmonic energy of a deformation field
 *
 * This class computes the harmonic energy of a deformation
 * field which is a measure inversely related to the smoothness
 * of the deformation field
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 *
 * This implementation was taken from the Insight Journal paper:
 * http://hdl.handle.net/1926/510
 *
 * \ingroup Operators
 * \ingroup ITK-Review
 */
template< class TInputImage >
class ITK_EXPORT WarpHarmonicEnergyCalculator:public Object
{
public:
  /** Standard class typedefs. */
  typedef WarpHarmonicEnergyCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer< Self >         Pointer;
  typedef SmartPointer< const Self >   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(WarpHarmonicEnergyCalculator, Object);

  /** Type definition for the input image. */
  typedef TInputImage ImageType;

  /** Pointer type for the image. */
  typedef typename TInputImage::Pointer ImagePointer;

  /** Const Pointer type for the image. */
  typedef typename TInputImage::ConstPointer ImageConstPointer;

  /** Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType PixelType;

  /** Type definition for the input image index type. */
  typedef typename TInputImage::IndexType IndexType;

  /** Type definition for the input image region type. */
  typedef typename TInputImage::RegionType RegionType;

  /** The dimensionality of the input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Length of the vector pixel type of the input image. */
  itkStaticConstMacro(VectorDimension, unsigned int,
                      PixelType::Dimension);

  /** Type of the iterator that will be used to move through the image.  Also
      the type which will be passed to the evaluate function */
  typedef ConstNeighborhoodIterator< ImageType >             ConstNeighborhoodIteratorType;
  typedef typename ConstNeighborhoodIteratorType::RadiusType RadiusType;

  /** Set the derivative weights according to the spacing of the input image
   *  (1/spacing). Use this option if you want to calculate the Jacobian
   *  determinant in the space in which the data was acquired. */
  void SetUseImageSpacingOn()
  {
    this->SetUseImageSpacing(true);
  }

  /** Reset the derivative weights to ignore image spacing.  Use this option if
   *  you want to calculate the Jacobian determinant in the image space.
   *  Default is ImageSpacingOn. */
  void SetUseImageSpacingOff()
  {
    this->SetUseImageSpacing(false);
  }

  /** Set/Get whether or not the filter will use the spacing of the input
   *  image in its calculations */
  void SetUseImageSpacing(bool);

  itkGetConstMacro(UseImageSpacing, bool);

  typedef FixedArray< double, ImageDimension > WeightsType;

  /** Directly Set/Get the array of weights used in the gradient calculations.
   *  Note that calling UseImageSpacingOn will clobber these values. */
  itkSetMacro(DerivativeWeights, WeightsType);
  itkGetConstReferenceMacro(DerivativeWeights, WeightsType);

  /** Set the input image. */
  itkSetConstObjectMacro(Image, ImageType);

  /** Compute the minimum and maximum values of intensity of the input image. */
  void Compute(void);

  /** Return the smoothness value. */
  itkGetConstMacro(HarmonicEnergy, double);

  /** Set the region over which the values will be computed */
  void SetRegion(const RegionType & region);

protected:
  WarpHarmonicEnergyCalculator();
  virtual ~WarpHarmonicEnergyCalculator() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Get/Set the neighborhood radius used for gradient computation */
  itkGetConstReferenceMacro(NeighborhoodRadius, RadiusType);
  itkSetMacro(NeighborhoodRadius, RadiusType);

  double EvaluateAtNeighborhood(ConstNeighborhoodIteratorType & it) const;

private:
  WarpHarmonicEnergyCalculator(const Self &); //purposely not implemented
  void operator=(const Self &);               //purposely not implemented

  double            m_HarmonicEnergy;
  ImageConstPointer m_Image;

  RegionType m_Region;
  bool       m_RegionSetByUser;

  bool m_UseImageSpacing;

  /** The weights used to scale partial derivatives during processing */
  WeightsType m_DerivativeWeights;

  RadiusType m_NeighborhoodRadius;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWarpHarmonicEnergyCalculator.txx"
#endif

#endif /* __itkWarpHarmonicEnergyCalculator_h */
