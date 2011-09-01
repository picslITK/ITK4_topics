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
#ifndef __itkBSplineSmoothingOnUpdateDisplacementFieldTransform_h
#define __itkBSplineSmoothingOnUpdateDisplacementFieldTransform_h

#include "itkDisplacementFieldTransform.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"

namespace itk
{

/** \class BSplineSmoothingOnUpdateDisplacementFieldTransform
 * \brief Representation of a smooth deformation field  with B-splines.
 *
 * Although there already exists a B-spline transform in ITK which can be used
 * for processes such as image registration, if these processes involve a dense
 * sampling of an image a significant computational speed-up can be achieved
 * by densely sampling the B-spline transform prior to invoking transformations.
 *
 * This class takes as input a displacement field, smooths it on demand using
 * the specified B-spline parameters.  This represents an alternative approach
 * to B-spline (FFD) registration and is explained more in detail in the
 * reference given below.
 *
 * \author Nicholas J. Tustison
 *
 * \par REFERENCE
 * NJ Tustison, BB Avants, JC Gee, "Directly Manipulated Free-Form Deformation
 * Image Registration", IEEE Transactions on Image Processing, 18(3):624-635,
 * 2009.
 *
 * \ingroup ITKDisplacementField
 */
template<class TScalar, unsigned int NDimensions>
class ITK_EXPORT BSplineSmoothingOnUpdateDisplacementFieldTransform :
  public DisplacementFieldTransform<TScalar, NDimensions>
{
public:
  /** Standard class typedefs. */
  typedef BSplineSmoothingOnUpdateDisplacementFieldTransform    Self;
  typedef DisplacementFieldTransform<TScalar, NDimensions>      Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( BSplineSmoothingOnUpdateDisplacementFieldTransform, DisplacementFieldTransform );

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro( Self );

  /** Dimension of the domain spaces. */
  itkStaticConstMacro( Dimension, unsigned int, NDimensions );

  /** Types from superclass */
  typedef typename Superclass::ScalarType               ScalarType;
  typedef typename Superclass::DerivativeType           DerivativeType;
  typedef typename Superclass::DisplacementFieldType    DisplacementFieldType;

  /**
   * typedefs for projecting the input displacement field onto a
   * B-spline field.
   */
  typedef typename DisplacementFieldType::PixelType                                         DisplacementVectorType;
  typedef DisplacementFieldType                                                             ControlPointLatticeType;
  typedef PointSet<DisplacementVectorType, Dimension>                                       PointSetType;
  typedef unsigned int                                                                      SplineOrderType;
  typedef BSplineScatteredDataPointSetToImageFilter<PointSetType, DisplacementFieldType>    BSplineFilterType;
  typedef typename BSplineFilterType::WeightsContainerType                                  WeightsContainerType;
  typedef typename BSplineFilterType::PointDataImageType                                    DisplacementFieldControlPointLatticeType;
  typedef typename BSplineFilterType::ArrayType                                             ArrayType;
  typedef typename BSplineFilterType::OutputImageType                                       OutputImageType;
  typedef typename ArrayType::ValueType                                                     ArrayValueType;

  /**
   * Update the transform's parameters by the values in \c update.  We
   * assume \c update is of the same length as Parameters. Throw exception
   * otherwise. The update process performs an smoothing on the displacement
   * field by using BSplines.
   * \c factor is a scalar multiplier for each value in update.
   * \c BSplineSmoothDisplacementField is called after the update is
   * added to the field.
   * See base class for more details.
   */
  virtual void UpdateTransformParameters( DerivativeType & update, ScalarType factor = 1.0 );

  /**
   * Set the spline order defining the bias field estimate.  Default = 3.
   */
  itkSetMacro( SplineOrder, SplineOrderType );

  /**
   * Get the spline order defining the bias field estimate.  Default = 3.
   */
  itkGetConstMacro( SplineOrder, SplineOrderType );

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
   * grid size.  The mesh size in each dimension is calculated as the
   * difference between the control point grid size and the spline order, i.e.
   * meshSize = controlPointGridSize - SplineOrder.
   */
  void SetMeshSize( const ArrayType & );

  /**
   * Set the number of fitting levels.  This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.  Default = 1 levels.
   */
  itkSetMacro( NumberOfFittingLevelsPerDimension, ArrayType );

  /**
   * Set the number of fitting levels.  This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.  Default = 1 level.
   */
  void SetNumberOfFittingLevels( const ArrayValueType );

  /**
   * Get the number of fitting levels.  One of the contributions is the
   * introduction of a multi-scale approach to fitting. This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.
   * Default = 1 levels.
   */
  itkGetConstMacro( NumberOfFittingLevelsPerDimension, ArrayType );

  /**
   * Enforce zero motion on the transform domain boundaries. Default = true.
   */
  itkSetMacro( EnforceStationaryBoundary, bool );

  /**
   * Enforce zero motion on the transform domain boundaries. Default = true.
   */
  itkGetConstMacro( EnforceStationaryBoundary, bool );

  /**
   * Enforce zero motion on the transform domain boundaries. Default = true.
   */
  itkBooleanMacro( EnforceStationaryBoundary );


protected:
  BSplineSmoothingOnUpdateDisplacementFieldTransform();
  virtual ~BSplineSmoothingOnUpdateDisplacementFieldTransform();

  void PrintSelf( std::ostream& os, Indent indent ) const;

  /**
   * Smooth the displacement field in-place using B-splines.
   * \warning Not thread safe. Does its own threading.
   */
  virtual void BSplineSmoothDisplacementField();

private:
  BSplineSmoothingOnUpdateDisplacementFieldTransform( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  SplineOrderType             m_SplineOrder;
  ArrayType                   m_NumberOfControlPoints;
  ArrayType                   m_NumberOfFittingLevelsPerDimension;
  bool                        m_EnforceStationaryBoundary;
};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkBSplineSmoothingOnUpdateDisplacementFieldTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkBSplineSmoothingOnUpdateDisplacementFieldTransform.hxx"
#endif

#endif // __itkBSplineSmoothingOnUpdateDisplacementFieldTransform_h
