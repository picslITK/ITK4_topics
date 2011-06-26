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
#ifndef __itkTransformVirtualDomainCalculator_h
#define __itkTransformVirtualDomainCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"

namespace itk
{
/** \class TransformVirtualDomainCalculator
 *
 * \ingroup Transforms
 */
template<class TInputImage1, class TInputImage2>
class ITK_EXPORT TransformVirtualDomainCalculator:public Object
{
public:
  /** Standard class typedefs. */
  typedef TransformVirtualDomainCalculator      Self;
  typedef Object                                Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TransformVirtualDomainCalculator, Object );

  /** Dimension of parameters. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage1::ImageDimension );

  /** Image Types to use in the initialization of the transform */
  typedef TInputImage1                           Image1Type;
  typedef typename Image1Type::ConstPointer      Image1Pointer;
  typedef TInputImage2                           Image2Type;
  typedef typename Image2Type::ConstPointer      Image2Pointer;

  typedef typename Image1Type::PointType         OriginType;
  typedef typename Image1Type::DirectionType     DirectionType;
  typedef typename Image1Type::IndexType         IndexType;

  /** Set image1 */
  itkSetConstObjectMacro( InputImage1, Image1Type );

  /** Set image2 */
  itkSetConstObjectMacro( InputImage2, Image2Type );

  /** Get the virtual domain origin */
  itkGetConstMacro( VirtualDomainOrigin, OriginType )

  /** Get the virtual domain direction */
  itkGetConstMacro( VirtualDomainDirection, DirectionType );

  /** Get the physically consistent origin for image 1 */
  itkGetConstMacro( PhysicallyConsistentOrigin1, OriginType )

  /** Get the physically consistent direction for image 1 */
  itkGetConstMacro( PhysicallyConsistentDirection1, DirectionType );

  /** Get the physically consistent origin for image 2 */
  itkGetConstMacro( PhysicallyConsistentOrigin2, OriginType )

  /** Get the physically consistent direction for image 2 */
  itkGetConstMacro( PhysicallyConsistentDirection2, DirectionType );

  /**
   * Calculate the virtual domain transform parameters once both images
   * have been set.
   */
  void CalculateVirtualDomainParameters();

  /**
   * The user might want to use physically consistent origin and direction
   * parameters in estimating a virtual domain.  For example, a given image,
   * when passed through FlipAxisImageFilter or PermuteAxesImageFilter,
   * will reside in physical space the same way before and after the filtering
   * operation.  However, the origin and direction will be different.
   * Estimating a virtual domain that resides between the two images would
   * produce undesirable results.  Default = true.
   */
  itkSetMacro( UsePhysicalConsistency, bool );
  itkGetConstMacro( UsePhysicalConsistency, bool );
  itkBooleanMacro( UsePhysicalConsistency );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( SameDimensionCheck1,
    ( Concept::SameDimension<ImageDimension, TInputImage1::ImageDimension> ) );
  itkConceptMacro( SameDimensionCheck2,
    ( Concept::SameDimension<ImageDimension, TInputImage2::ImageDimension> ) );
  /** End concept checking */
#endif

protected:
  TransformVirtualDomainCalculator();
  ~TransformVirtualDomainCalculator(){}

  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  void operator=( const Self & );                   //purposely not
                                                    // implemented

  /**
   * Calculate the unbiased direction cosine matrix that resides between
   * two direction cosines matrices.  This uses the method of
   * Denman, E., Beavers, A., 1976. The matrix sign function and computations
   * in systems. Appl. Math. Comput. 2 (1), 63–94.
   */
  DirectionType CalculateUnbiasedDirectionCosineMatrix(
    DirectionType, DirectionType );

  /**
   * Calculate a physically consistent origin and directions based on
   * how the image resides in physical space.
   */
  void CalculatePhysicallyConsistentDomainParameters( Image1Pointer,
    OriginType &, DirectionType & );

  Image1Pointer                                 m_InputImage1;
  Image2Pointer                                 m_InputImage2;

  OriginType                                    m_VirtualDomainOrigin;
  DirectionType                                 m_VirtualDomainDirection;

  bool                                          m_UsePhysicalConsistency;
  OriginType                                    m_PhysicallyConsistentOrigin1;
  OriginType                                    m_PhysicallyConsistentOrigin2;
  DirectionType                                 m_PhysicallyConsistentDirection1;
  DirectionType                                 m_PhysicallyConsistentDirection2;

}; //class TransformVirtualDomainCalculator
}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTransformVirtualDomainCalculator.txx"
#endif

#endif /* __itkTransformVirtualDomainCalculator_h */
