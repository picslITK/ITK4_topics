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
#ifndef __itkPointSetToPointSetMetric_h
#define __itkPointSetToPointSetMetric_h

#include "itObjectToObjectMetric.h"

#include "itkFixedArray.h"
#include "itkTransform.h"

namespace itk
{
/** \class PointSetToPointSetMetric
 * \brief Computes similarity between two point sets.
 *
 * This Class is templated over the type of the two point-sets.  It
 * expects a Transform to be plugged in.  This particular
 * class is the base class for a hierarchy of point-set to point-set metrics.
 *
 * This class computes a value that measures the similarity between the fixed
 * point-set and the transformed moving point-set.
 *
 * \ingroup RegistrationMetrics
 *
 * \ingroup ITK-RegistrationCommon
 */

template<class TFixedPointSet,  class TMovingPointSet>
class ITK_EXPORT PointSetToPointSetMetric
: public ObjectToObjectMetric
{
public:

  /** Standard class typedefs. */
  typedef PointSetToPointSetMetric      Self;
  typedef ObjectToObjectMetric          Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Type used for representing point components  */
  typedef Superclass::ParametersValueType CoordinateRepresentationType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( PointSetToPointSetMetric, ObjectToObjectMetric );

  /**  Type of the moving Pointset. */
  typedef TMovingPointSet                           MovingPointSetType;
  typedef TMovingPointSet::PointType                MovingPointType;
  typedef typename TMovingPointSet::PixelType       MovingPointSetPixelType;
  typedef typename MovingPointSetType::ConstPointer MovingPointSetConstPointer;

  /**  Type of the fixed Pointset. */
  typedef TFixedPointSet                           FixedPointSetType;
  typedef TFixedPointSet::PointType                FixedPointType;
  typedef typename FixedPointSetType::ConstPointer FixedPointSetConstPointer;

  /** Constants for the pointset dimensions */
  itkStaticConstMacro( MovingPointSetDimension, unsigned int,
    TMovingPointSet::PointDimension );
  itkStaticConstMacro( FixedPointSetDimension, unsigned int,
    TFixedPointSet::PointDimension );

  typedef typename FixedPointSetType::PointsContainer::ConstIterator
    FixedPointSetIteratorType;
  typedef typename FixedPointSetType::PointDataContainer::ConstIterator
    FixedPointSetDataIteratorType;

  typedef typename MovingPointSetType::PointsContainer::ConstIterator
    MovingPointSetIteratorType;
  typedef typename MovingPointSetType::PointDataContainer::ConstIterator
    MovingPointSetDataIteratorType;

  /**  Type of the Transform Base class */
  typedef Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( MovingPointSetDimension ),
    itkGetStaticConstMacro( FixedPointSetDimension )>  MovingTransformType;

  typedef Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( FixedPointSetDimension ),
    itkGetStaticConstMacro( MovingPointSetDimension )> FixedTransformType;

  typedef typename TransformType::Pointer         TransformPointer;
  typedef typename TransformType::InputPointType  InputPointType;
  typedef typename TransformType::OutputPointType OutputPointType;
  typedef typename TransformType::ParametersType  TransformParametersType;
  typedef typename TransformType::JacobianType    TransformJacobianType;

  /**  Type of the measure. */
  typedef Superclass::MeasureType MeasureType;

  /**  Type of the derivative. */
  typedef Superclass::DerivativeType DerivativeType;

  typedef FixedArray<typename DerivativeType::ValueType,
    itkGetStaticConstMacro( FixedPointSetDimension )>  LocalFixedDerivativeType;

  typedef MovingArray<typename DerivativeType::ValueType,
    itkGetStaticConstMacro( MovingPointSetDimension )> LocalMovingDerivativeType;

  /**  Type of the parameters. */
  typedef Superclass::ParametersType ParametersType;

  /** Type of coordinate system used to calculate values, derivatives */
  enum CoordinateSystemType { Fixed, Moving, Both };

  /** Connect the fixed pointset.  */
  itkSetConstObjectMacro( FixedPointSet, FixedPointSetType );

  /** Get the fixed point set. */
  itkGetConstObjectMacro( FixedPointSet, FixedPointSetType );

  /** Connect the fixed transform. */
  itkSetConstOjbectMacro( FixedTransform, FixedTransformType );

  /** Get a pointer to the fixed transform.  */
  itkGetObjectMacro( FixedTransform, FixedTransformType );

  /** Connect the moving point set.  */
  itkSetConstObjectMacro( MovingPointSet, MovingPointSetType );

  /** Get the moving point set. */
  itkGetConstObjectMacro( MovingPointSet, MovingPointSetType );

  /** Connect the moving transform. */
  itkSetConstObjectMacro( MovingTransform, MovingTransformType );

  /** Get a pointer to the moving transform.  */
  itkGetObjectMacro( MovingTransform, MovingTransformType );

  /**
   * This method returns the value of the metric based on the current
   * transformation(s).  This function can be redefined in derived classes
   * but many point set metrics follow the same structure---one iterates
   * through the points and, for each point a metric value is calculated.
   * The summation of these individual point metric values gives the total
   * value of the metric.  Note that this might not be applicable to all
   * point set metrics.  For those cases, the developer will have to redefine
   * the GetValue() function.
   */
  virtual MeasureType GetValue();

  /**
   * This method returns the derivative based on the current
   * transformation(s).  This function can be redefined in derived classes
   * but many point set metrics follow the same structure---one iterates
   * through the points and, for each point a derivative is calculated.
   * The set of all these local derivatives constitutes the total derivative.
   * Note that this might not be applicable to all point set metrics.  For
   * those cases, the developer will have to redefine the GetDerivative()
   * function.
   */
  virtual void GetDerivative( DerivativeType & derivative ) const;

  /**
   * This method returns the derivative and value based on the current
   * transformation(s).  This function can be redefined in derived classes
   * but many point set metrics follow the same structure---one iterates
   * through the points and, for each point a derivative and value is calculated.
   * The set of all these local derivatives/values constitutes the total
   * derivative and value.  Note that this might not be applicable to all
   * point set metrics.  For those cases, the developer will have to redefine
   * the GetValue() and GetDerivative() functions.
   */
  virtual void GetValueAndDerivative( MeasureType & value,
    DerivativeType & derivative ) const;

  /**
   * Function to be defined in the appropriate derived classes.  Calculates
   * the local metric value for a single point within the fixed point set.
   */
  virtual MeasureType GetLocalFixedNeighborhoodValue( const FixedPointType )
    const = 0;

  /**
   * Function to be defined in the appropriate derived classes.  Calculates
   * the local metric value for a single point within the moving point set.
   */
  virtual MeasureType GetLocalMovingNeighborhoodValue( const MovingPointType )
    const = 0;

  /**
   * Function to be defined in the appropriate derived classes.  Calculates
   * the local derivative for a single point within the fixed point set.
   */
  virtual LocalFixedDerivativeType GetLocalFixedNeighborhoodDerivative(
    const FixedPointType ) const = 0;

  /**
   * Function to be defined in the appropriate derived classes.  Calculates
   * the local derivative for a single point within the moving point set.
   */
  virtual LocalMovingDerivativeType GetLocalMovingNeighborhoodDerivative(
    const MovingPointType ) const = 0;

  /**
   * Calculates the local value/derivative for a single point within the fixed
   * point set.
   */
  virtual void GetLocalFixedNeighborhoodDerivative( const FixedPointType,
    MeasuresType &, LocalFixedDerivativeType & ) const = 0;

  /**
   * Calculates the local value/derivative for a single point within the moving
   * point set.
   */
  virtual void GetLocalMovingNeighborhoodDerivative( const MovingPointType,
    MeasuresType &, LocalMovingDerivativeType & ) const = 0;

  /**
   * Set coordinate system type.  This variable allows the user to switch
   * between calculating the value and derivative with respect to the fixed
   * point set or moving point set.
   */
  itkSetMacro( CoordinateSystem, CoordinateSystemType );

  /**
   * Get coordinate system type.
   */
  itGetConstMacro( CoordinateSystem, CoordinateSystemType );

  /** Initialize the metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize( void )
    throw ( ExceptionObject );

protected:
  PointSetToPointSetMetric();
  virtual ~PointSetToPointSetMetric() {}
  void PrintSelf( std::ostream & os, Indent indent ) const;

  FixedPointSetConstPointer  m_FixedPointSet;
  MovingPointSetPointer      m_FixedTransformedPointSet;
  mutable TransformPointer   m_FixedTransform;

  MovingPointSetConstPointer m_MovingPointSet;
  FixedPointSetConstPointer  m_MovingTransformedPointSet;
  mutable TransformPointer   m_MovingTransform;

  CoordinateSystemType       m_CoordinateSystem;

  /**
   * Warp the fixed point set based on the fixed transform.  Note that the
   * warped moving point set is of type FixedPointSetType since the transform
   * takes the points from the fixed to the moving domain.
   */
  void TransformFixedPointSet();

  /**
   * Warp the moving point set based on the moving transform.  Note that the
   * warped moving point set is of type FixedPointSetType since the transform
   * takes the points from the moving to the fixed domain.
   */
  void TransformMovingPointSet();

private:
  PointSetToPointSetMetric( const Self & ); //purposely not implemented
  void operator=( const Self & );           //purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetToPointSetMetric.txx"
#endif

#endif
