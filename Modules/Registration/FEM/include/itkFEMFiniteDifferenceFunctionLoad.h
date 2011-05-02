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
#ifndef __itkFEMFiniteDifferenceFunctionLoad_h
#define __itkFEMFiniteDifferenceFunctionLoad_h

#include "itkFEMLoadElementBase.h"

#include "itkImage.h"
#include "itkTranslationTransform.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkDerivativeOperator.h"
#include "itkForwardDifferenceOperator.h"
#include "vnl/vnl_math.h"

#include "itkDemonsRegistrationFunction.h"
#include "itkMeanSquareRegistrationFunction.h"
#include "itkNCCRegistrationFunction.h"
#include "itkMIRegistrationFunction.h"

namespace itk
{
namespace fem
{
/**
 * \class FiniteDifferenceFunctionLoad
 * \brief General image pair load that uses the itkFiniteDifferenceFunctions.
 *
 * This load computes FEM gravity loads by using derivatives provided
 * by itkFiniteDifferenceFunctions (e.g. mean squares intensity difference.)
 * The function responsible for this is called Fg, as required by the FEMLoad
 * standards.  It takes a vnl_vector as input.
 * We assume the vector input is of size 2*ImageDimension.
 * The 0 to ImageDimension-1 elements contain the position, p,
 * in the reference (moving) image.  The next ImageDimension to 2*ImageDimension-1
 * elements contain the value of the vector field at that point, v(p).
 * The metrics return both a scalar similarity value and vector-valued derivative.
 * The derivative is what gives us the force to drive the FEM registration.
 * These values are computed with respect to some region in the Fixed image.
 * This region size may be set by the user by calling SetMetricRadius.
 * As the metric derivative computation evolves, performance should improve
 * and more functionality will be available (such as scale selection).
 * \ingroup ITK-FEMRegistration
 */
template< class TMoving, class TFixed >
class ITK_EXPORT FiniteDifferenceFunctionLoad:public LoadElement
{
  FEM_CLASS(FiniteDifferenceFunctionLoad, LoadElement)
public:

  // Necessary typedefs for dealing with images BEGIN
  typedef typename LoadElement::Float Float;

  typedef TMoving                                MovingImageType;
  typedef typename MovingImageType::ConstPointer MovingConstPointer;
  typedef MovingImageType *                      MovingPointer;
  typedef TFixed                                 FixedImageType;
  typedef FixedImageType *                       FixedPointer;
  typedef typename FixedImageType::ConstPointer  FixedConstPointer;

  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      MovingImageType::ImageDimension);

  typedef ImageRegionIteratorWithIndex< MovingImageType > MovingRegionIteratorType;
  typedef ImageRegionIteratorWithIndex< FixedImageType >  FixedRegionIteratorType;

  typedef NeighborhoodIterator< MovingImageType >
  MovingNeighborhoodIteratorType;
  typedef typename MovingNeighborhoodIteratorType::IndexType
  MovingNeighborhoodIndexType;
  typedef typename MovingNeighborhoodIteratorType::RadiusType
  MovingRadiusType;
  typedef typename MovingNeighborhoodIteratorType::RadiusType
  RadiusType;
  typedef NeighborhoodIterator< FixedImageType >
  FixedNeighborhoodIteratorType;
  typedef typename FixedNeighborhoodIteratorType::IndexType
  FixedNeighborhoodIndexType;
  typedef typename FixedNeighborhoodIteratorType::RadiusType
  FixedRadiusType;

// IMAGE DATA
  typedef   typename  MovingImageType::PixelType MovingPixelType;
  typedef   typename  FixedImageType::PixelType  FixedPixelType;
  typedef   Float                                PixelType;
  typedef   Float                                ComputationType;
  typedef   Image< PixelType, itkGetStaticConstMacro(ImageDimension) >
  ImageType;
  typedef   itk::Vector< float, itkGetStaticConstMacro(ImageDimension) >
  VectorType;
  typedef   vnl_vector< Float > FEMVectorType;
  typedef   Image< VectorType, itkGetStaticConstMacro(ImageDimension) >
  DeformationFieldType;
  typedef   typename DeformationFieldType::Pointer DeformationFieldTypePointer;

  typedef NeighborhoodIterator< DeformationFieldType >
  FieldIteratorType;

// Necessary typedefs for dealing with images END

  /** PDEDeformableRegistrationFilterFunction type. */
  typedef PDEDeformableRegistrationFunction< FixedImageType, MovingImageType,
                                             DeformationFieldType >
  FiniteDifferenceFunctionType;
  typedef typename FiniteDifferenceFunctionType::Pointer FiniteDifferenceFunctionTypePointer;

  typedef typename FiniteDifferenceFunctionType::TimeStepType TimeStepType;

  typedef MeanSquareRegistrationFunction< FixedImageType, MovingImageType,
                                          DeformationFieldType >  MeanSquareRegistrationFunctionType;

  typedef DemonsRegistrationFunction< FixedImageType, MovingImageType,
                                      DeformationFieldType >  DemonsRegistrationFunctionType;

  typedef NCCRegistrationFunction< FixedImageType, MovingImageType,
                                   DeformationFieldType >  NCCRegistrationFunctionType;

  typedef MIRegistrationFunction< FixedImageType, MovingImageType,
                                  DeformationFieldType >  MIRegistrationFunctionType;

// FUNCTIONS

  /* This method sets the pointer to a FiniteDifferenceFunction object that
   * will be used by the filter to calculate updates at image pixels.
   * \returns A FiniteDifferenceObject pointer. */
  void SetDifferenceFunction(FiniteDifferenceFunctionTypePointer drfp)
  {
    drfp->SetFixedImage(m_FixedImage);
    drfp->SetMovingImage(m_MovingImage);
    drfp->SetRadius(m_MetricRadius);
    drfp->SetDeformationField(m_DeformationField);
    drfp->InitializeIteration();
    this->m_DifferenceFunction = drfp;
  }

  void SetMetric(FiniteDifferenceFunctionTypePointer drfp)
  {
    this->SetDifferenceFunction( static_cast< FiniteDifferenceFunctionType * >(
                                   drfp.GetPointer() ) );

    m_FixedSize = m_DeformationField->GetLargestPossibleRegion().GetSize();
  }

  /** Define the reference (moving) image. */
  void SetMovingImage(MovingImageType *R)
  {
    m_MovingImage = R;
    m_MovingSize = m_MovingImage->GetLargestPossibleRegion().GetSize();
    if ( this->m_DifferenceFunction ) { this->m_DifferenceFunction->SetMovingImage(m_MovingImage); }
  }

  /** Define the target (fixed) image. */
  void SetFixedImage(FixedImageType *T)
  {
    m_FixedImage = T;
    m_FixedSize = T->GetLargestPossibleRegion().GetSize();
    if ( this->m_DifferenceFunction )
      {
      this->m_DifferenceFunction->SetFixedImage(m_MovingImage);
      }
  }

  MovingPointer GetMovingImage()
  {
    return m_MovingImage;
  }

  FixedPointer GetFixedImage()
  {
    return m_FixedImage;
  }

  /** Define the metric region size. */
  void SetMetricRadius(MovingRadiusType T)
  {
    m_MetricRadius  = T;
  }

  /** Get the metric region size. */
  MovingRadiusType GetMetricRadius()
  {
    return m_MetricRadius;
  }

  /** Set/Get methods for the number of integration points to use
   * in each 1-dimensional line integral when evaluating the load.
   * This value is passed to the load implementation.
   */
  void SetNumberOfIntegrationPoints(unsigned int i)
  {
    m_NumberOfIntegrationPoints = i;
  }

  unsigned int GetNumberOfIntegrationPoints()
  {
    return m_NumberOfIntegrationPoints;
  }

  /** Set the direction of the gradient (uphill or downhill).
    * E.g. the mean squares metric should be minimized while NCC and PR should be maximized.
    */
  void SetSign(Float s)
  {
    m_Sign = s;
  }

  /** Set the sigma in a gaussian measure. */
  void SetTemp(Float s)
  {
    m_Temp = s;
  }

  /** Scaling of the similarity energy term */
  void SetGamma(Float s)
  {
    m_Gamma = s;
  }

  void SetSolution(Solution::ConstPointer ptr)
  {
    m_Solution = ptr;
  }

  Solution::ConstPointer GetSolution()
  {
    return m_Solution;
  }

  /** \todo FIXME - WE ASSUME THE 2ND VECTOR (INDEX 1) HAS THE INFORMATION WE WANT */
  Float GetSolution(unsigned int i, unsigned int which = 0)
  {
    return m_Solution->GetSolutionValue(i, which);
  }

  FiniteDifferenceFunctionLoad(); // cannot be private until we always use smart
                                  // pointers
  Float EvaluateMetricGivenSolution(Element::ArrayType *el, Float step = 1.0);

  /**
   * Compute the image based load - implemented with ITK metric derivatives.
   */
  VectorType    Fe1(VectorType);
  FEMVectorType Fe(FEMVectorType, FEMVectorType);

  static Baseclass * NewFiniteDifferenceFunctionLoad(void)
  { return new FiniteDifferenceFunctionLoad; }

  /** Set the  */
  void SetDeformationField(DeformationFieldTypePointer df)
  { m_DeformationField = df; }

  /** Get the  */
  DeformationFieldTypePointer GetDeformationField() { return m_DeformationField; }
  void InitializeIteration();

  void InitializeMetric();

  void PrintCurrentEnergy();

  double GetCurrentEnergy();

  void  SetCurrentEnergy(double e = 0.0);

protected:
private:
  MovingPointer    m_MovingImage;
  FixedPointer     m_FixedImage;
  MovingRadiusType m_MetricRadius;                                   /** used by
                                                                       the
                                                                       metric to
                                                                       set
                                                                       region
                                                                       size for
                                                                       fixed
                                                                       image*/
  typename MovingImageType::SizeType m_MovingSize;

  typename FixedImageType::SizeType m_FixedSize;

  unsigned int m_NumberOfIntegrationPoints;
  unsigned int m_SolutionIndex;
  unsigned int m_SolutionIndex2;
  Float        m_Temp;
  Float        m_Gamma;

  typename Solution::ConstPointer m_Solution;

  float                               m_GradSigma;
  float                               m_Sign;
  float                               m_WhichMetric;
  FiniteDifferenceFunctionTypePointer m_DifferenceFunction;

  typename DeformationFieldType::Pointer m_DeformationField;
  /** Dummy static int that enables automatic registration
      with FEMObjectFactory. */
  static const int m_DummyCLID;
};
}
} // end namespace fem/itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFEMFiniteDifferenceFunctionLoad.txx"
#endif

#endif
