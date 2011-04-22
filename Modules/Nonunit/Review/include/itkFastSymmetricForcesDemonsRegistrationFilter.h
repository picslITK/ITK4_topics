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
#ifndef __itkFastSymmetricForcesDemonsRegistrationFilter_h
#define __itkFastSymmetricForcesDemonsRegistrationFilter_h

#include "itkPDEDeformableRegistrationFilter.h"
#include "itkESMDemonsRegistrationFunction.h"

#include "itkMultiplyByConstantImageFilter.h"
#include "itkExponentialDeformationFieldImageFilter.h"

namespace itk
{
/** \class FastSymmetricForcesDemonsRegistrationFilter
 * \brief Deformably register two images using a symmetric forces demons algorithm.
 *
 * This class was contributed by Tom Vercauteren, INRIA & Mauna Kea Technologies
 * based on a variation of the DemonsRegistrationFilter.
 *
 * FastSymmetricForcesDemonsRegistrationFilter implements the demons deformable algorithm that
 * register two images by computing the deformation field which will map a
 * moving image onto a fixed image.
 *
 * A deformation field is represented as a image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the fixed image. The vector type must support element access via operator
 * []. It is assumed that the vector elements behave like floating point
 * scalars.
 *
 * This class is templated over the fixed image type, moving image type
 * and the deformation field type.
 *
 * The input fixed and moving images are set via methods SetFixedImage
 * and SetMovingImage respectively. An initial deformation field maybe set via
 * SetInitialDeformationField or SetInput. If no initial field is set,
 * a zero field is used as the initial condition.
 *
 * The output deformation field can be obtained via methods GetOutput
 * or GetDeformationField.
 *
 * This class make use of the finite difference solver hierarchy. Update
 * for each iteration is computed in DemonsRegistrationFunction.
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 *
 * This implementation was taken from the Insight Journal paper:
 * http://hdl.handle.net/1926/510
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa DemonsRegistrationFilter
 * \sa DemonsRegistrationFunction
 * \ingroup DeformableImageRegistration MultiThreaded
 * \ingroup ITK-Review
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
class ITK_EXPORT FastSymmetricForcesDemonsRegistrationFilter:
  public PDEDeformableRegistrationFilter< TFixedImage, TMovingImage,
                                          TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef FastSymmetricForcesDemonsRegistrationFilter                                     Self;
  typedef PDEDeformableRegistrationFilter< TFixedImage, TMovingImage, TDeformationField > Superclass;
  typedef SmartPointer< Self >                                                            Pointer;
  typedef SmartPointer< const Self >                                                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FastSymmetricForcesDemonsRegistrationFilter,
               PDEDeformableRegistrationFilter);

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType    FixedImageType;
  typedef typename Superclass::FixedImagePointer FixedImagePointer;

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename Superclass::MovingImagePointer MovingImagePointer;

  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType    DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer DeformationFieldPointer;

  /** Get the metric value. The metric value is the mean square difference
   * in intensity between the fixed image and transforming moving image
   * computed over the the overlapping region between the two images.
   * This value is calculated for the current iteration */
  virtual double GetMetric() const;

  virtual const double & GetRMSChange() const;

  /** DemonsRegistrationFilterFunction type.
   *
   *  FIXME: Why is this the only permissible function ?
   *
   */
  typedef ESMDemonsRegistrationFunction<
    FixedImageType,
    MovingImageType, DeformationFieldType >                DemonsRegistrationFunctionType;

  typedef typename DemonsRegistrationFunctionType::GradientType GradientType;
  virtual void SetUseGradientType(GradientType gtype);

  virtual GradientType GetUseGradientType() const;

  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void SetIntensityDifferenceThreshold(double);

  virtual double GetIntensityDifferenceThreshold() const;

  virtual void SetMaximumUpdateStepLength(double);

  virtual double GetMaximumUpdateStepLength() const;

protected:
  FastSymmetricForcesDemonsRegistrationFilter();
  ~FastSymmetricForcesDemonsRegistrationFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Initialize the state of filter and equation before each iteration. */
  virtual void InitializeIteration();

  /** This method allocates storage in m_UpdateBuffer.  It is called from
   * FiniteDifferenceFilter::GenerateData(). */
  virtual void AllocateUpdateBuffer();

  /** FiniteDifferenceFunction type. */
  typedef typename
  Superclass::FiniteDifferenceFunctionType FiniteDifferenceFunctionType;

  /** Take timestep type from the FiniteDifferenceFunction. */
  typedef typename
  FiniteDifferenceFunctionType::TimeStepType TimeStepType;

  /** Apply update. */
  virtual void ApplyUpdate(TimeStepType dt);

  /** other typedefs */
  typedef MultiplyByConstantImageFilter<
    DeformationFieldType,
    TimeStepType, DeformationFieldType >                  MultiplyByConstantType;

  typedef AddImageFilter<
    DeformationFieldType,
    DeformationFieldType, DeformationFieldType >          AdderType;

  typedef typename MultiplyByConstantType::Pointer MultiplyByConstantPointer;
  typedef typename AdderType::Pointer              AdderPointer;
private:
  FastSymmetricForcesDemonsRegistrationFilter(const Self &); //purposely not
                                                             // implemented
  void operator=(const Self &);                              //purposely not

  // implemented

  /** Downcast the DifferenceFunction using a dynamic_cast to ensure that it is of the correct type.
   * this method will throw an exception if the function is not of the expected type. */
  DemonsRegistrationFunctionType *  DownCastDifferenceFunctionType();

  const DemonsRegistrationFunctionType *  DownCastDifferenceFunctionType() const;

  MultiplyByConstantPointer m_Multiplier;
  AdderPointer              m_Adder;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFastSymmetricForcesDemonsRegistrationFilter.txx"
#endif

#endif
