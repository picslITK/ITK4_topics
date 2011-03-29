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
#ifndef __itkBSplineDeformationFieldTransform_h
#define __itkBSplineDeformationFieldTransform_h

#include "itkDeformationFieldTransform.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"

namespace itk
{

/** \class BSplineDeformationFieldTransform
 *
 * \ingroup ITK-Transform
 */
template
  <class TScalar, unsigned int NDimensions>
class ITK_EXPORT BSplineDeformationFieldTransform :
  public DeformationFieldTransform<TScalar, NDimensions>
{
public:
  /** Standard class typedefs. */
  typedef BSplineDeformationFieldTransform                  Self;
  typedef DeformationFieldTransform<TScalar, NDimensions>   Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( BSplineDeformationFieldTransform, Transform );

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro( Self );

  /** Dimension of the domain spaces. */
  itkStaticConstMacro( Dimension, unsigned int, NDimensions );

  /** Use the deformation field type */
  typedef typename Superclass::DeformationFieldType DeformationFieldType;

  /**
   * typedefs for projecting the input deformation field onto a B-spline field.
   */
  typedef typename DeformationFieldType::PixelType DisplacementVectorType;
  typedef DeformationFieldType                     ControlPointLatticeType;
  typedef PointSet<DisplacementVectorType,
    itkGetStaticConstMacro( Dimension )>           PointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <PointSetType, DeformationFieldType>           BSplineFilterType;
  typedef typename BSplineFilterType::PointDataImageType
                                                   DeformationFieldControlPointLatticeType;
  typedef typename BSplineFilterType::ArrayType    ArrayType;

  /** Get/Set the deformation field. */
  virtual void SetDeformationField( DeformationFieldType * );

  /**
   * Set the spline order defining the bias field estimate.  Default = 3.
   */
  itkSetMacro( SplineOrder, unsigned int );

  /**
   * Get the spline order defining the bias field estimate.  Default = 3.
   */
  itkGetConstMacro( SplineOrder, unsigned int );

  /**
   * Set the control point grid size definining the B-spline estimate of the
   * scalar bias field.  In each dimension, the B-spline mesh size is equal
   * to the number of control points in that dimension minus the spline order.
   * Default = 4 control points in each dimension for a mesh size of 1 in each
   * dimension.
   */
  itkSetMacro( NumberOfControlPoints, ArrayType );

  /**
   * Get the control point grid size definining the B-spline estimate of the
   * scalar bias field.  In each dimension, the B-spline mesh size is equal
   * to the number of control points in that dimension minus the spline order.
   * Default = 4 control points in each dimension for a mesh size of 1 in each
   * dimension.
   */
  itkGetConstMacro( NumberOfControlPoints, ArrayType );

  /**
   * Set the mesh size which is another way of specifying the control point
   * grid size.  The mesh size in each dimension is calculated as the difference
   * between the control point grid size and the spline order, i.e.
   * meshSize = controlPointGridSize - SplineOrder.
   */
  void SetMeshSize( const ArrayType meshSize )
    {
    ArrayType numberOfControlPoints;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      numberOfControlPoints[d] = meshSize[d] + this->m_SplineOrder;
      }
    this->SetNumberOfControlPoints( numberOfControlPoints );
    }

  /**
   * Set the number of fitting levels.  One of the contributions of N4 is the
   * introduction of a multi-scale approach to fitting. This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.  Default = 3 levels.
   */
  itkSetMacro( NumberOfFittingLevels, ArrayType );

  /**
   * Set the number of fitting levels.  One of the contributions of N4 is the
   * introduction of a multi-scale approach to fitting. This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.  Default = 3 levels.
   */
  void SetNumberOfFittingLevels( const unsigned int n )
    {
    ArrayType nlevels;

    nlevels.Fill( n );
    this->SetNumberOfFittingLevels( nlevels );
    }

  /**
   * Get the number of fitting levels.  One of the contributions is the
   * introduction of a multi-scale approach to fitting. This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.  Default = 3 levels.
   */
  itkGetConstMacro( NumberOfFittingLevels, ArrayType );


  /**
   * Boolean value to calculate an approximate inverse of the deformation field.
   */
  itkSetMacro( CalculateApproximateInverseDeformationField, bool );

  /**
   * Boolean value to calculate an approximate inverse of the deformation field.
   */
  itkGetConstMacro( CalculateApproximateInverseDeformationField, bool );

  /**
   * Boolean value to calculate an approximate inverse of the deformation field.
   */
  itkBooleanMacro( CalculateApproximateInverseDeformationField );

protected:
  BSplineDeformationFieldTransform();
  virtual ~BSplineDeformationFieldTransform();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  bool                   m_CalculateApproximateInverseDeformationField;

  typename DeformationFieldControlPointLatticeType::Pointer
                                              m_DeformationFieldControlPointLattice;
  typename DeformationFieldControlPointLatticeType::Pointer
                                              m_InverseDeformationFieldControlPointLattice;

  unsigned int                                m_SplineOrder;
  ArrayType                                   m_NumberOfControlPoints;
  ArrayType                                   m_NumberOfFittingLevels;

private:
  BSplineDeformationFieldTransform( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkBSplineDeformationFieldTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkBSplineDeformationFieldTransform.hxx"
#endif

#endif // __itkBSplineDeformationFieldTransform_h
