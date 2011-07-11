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
#ifndef __itkExpectationBasedPointSetMetric_txx
#define __itkExpectationBasedPointSetMetric_txx

#include "itkExpectationBasedPointSetMetric.h"

#include "itkArray.h"

namespace itk {

template<class TFixedPointSet, class TMovingPointSet>
ExpectationBasedPointSetMetric<TFixedPointSet, TMovingPointSet>
::ExpectationBasedPointSetMetric() :
  m_Sigma( 1.0 ),
  m_EvaluationKNeighborhood( 50 )
{
}

template<class TFixedPointSet, class TMovingPointSet>
typename ExpectationBasedPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
ExpectationBasedPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue( const PointType point ) const
{
  MeasureType localValue = 0.0;
  MeasureType preFactor = 1.0 / ( vcl_sqrt( 2 * vnl_math::pi ) * this->m_Sigma );

  if( this->GetDerivativeSource() == Superclass::Fixed )
    {
    NeighborsIdentifierType neighborhood;
    this->m_MovingTransformedPointsLocator->FindClosestNPoints( point,
      this->m_EvaluationKNeighborhood, neighborhood );

    typename NeighborsIdentifierType::const_iterator it;
    for( it = neighborhood.begin(); it != neighborhood.end(); ++it )
      {
      PointType neighbor = this->m_MovingTransformedPointSet->GetPoint( *it );
      MeasureType distance = point.SquaredEuclideanDistanceTo( neighbor );
      localValue += preFactor *
        vcl_exp( -distance / ( 2.0 * vnl_math_sqr( this->m_Sigma ) ) );
      }
    }
  else
    {
    NeighborsIdentifierType neighborhood;
    this->m_MovingTransformedPointsLocator->FindClosestNPoints( point,
      this->m_EvaluationKNeighborhood, neighborhood );

    typename NeighborsIdentifierType::const_iterator it;
    for( it = neighborhood.begin(); it != neighborhood.end(); ++it )
      {
      PointType neighbor = this->m_FixedTransformedPointSet->GetPoint( *it );
      MeasureType distance = point.SquaredEuclideanDistanceTo( neighbor );
      localValue += preFactor *
        vcl_exp( -distance / ( 2.0 * vnl_math_sqr( this->m_Sigma ) ) );
      }
    }
  return localValue;
}

template<class TFixedPointSet, class TMovingPointSet>
void
ExpectationBasedPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative( const PointType point,
  MeasureType &measure, LocalDerivativeType &localDerivative ) const
{
  Array<MeasureType> measureValues;
  measureValues.SetSize( this->m_EvaluationKNeighborhood );
  measureValues.Fill( 0.0 );

  measure = 0.0;
  MeasureType preFactor = 1.0 / ( vcl_sqrt( 2 * vnl_math::pi ) * this->m_Sigma );

  localDerivative.Fill( 0.0 );

  PointType weightedPoint;
  weightedPoint.Fill( 0.0 );

  NeighborsIdentifierType neighborhood;

  if( this->GetDerivativeSource() == Superclass::Fixed )
    {
    this->m_MovingTransformedPointsLocator->FindClosestNPoints( point,
      this->m_EvaluationKNeighborhood, neighborhood );

    typename NeighborsIdentifierType::const_iterator it;
    for( it = neighborhood.begin(); it != neighborhood.end(); ++it )
      {
      PointType neighbor = this->m_MovingTransformedPointSet->GetPoint( *it );
      MeasureType distance = point.SquaredEuclideanDistanceTo( neighbor );
      measureValues[it - neighborhood.begin()] = preFactor *
        vcl_exp( -distance / ( 2.0 * vnl_math_sqr( this->m_Sigma ) ) );
      measure += measureValues[it - neighborhood.begin()];
      }

    for( it = neighborhood.begin(); it != neighborhood.end(); ++it )
      {
      PointType neighbor = this->m_MovingTransformedPointSet->GetPoint( *it );
      typename PointType::VectorType neighborVector =
        neighbor.GetVectorFromOrigin();
      weightedPoint += ( neighborVector *
        measureValues[it - neighborhood.begin()] / measure );
      }
    }
  else
    {
    this->m_FixedTransformedPointsLocator->FindClosestNPoints( point,
      this->m_EvaluationKNeighborhood, neighborhood );

    typename NeighborsIdentifierType::const_iterator it;
    for( it = neighborhood.begin(); it != neighborhood.end(); ++it )
      {
      PointType neighbor = this->m_FixedTransformedPointSet->GetPoint( *it );
      MeasureType distance = point.SquaredEuclideanDistanceTo( neighbor );
      measureValues[it - neighborhood.begin()] = preFactor *
        vcl_exp( -distance / ( 2.0 * vnl_math_sqr( this->m_Sigma ) ) );
      measure += measureValues[it - neighborhood.begin()];
      }

    for( it = neighborhood.begin(); it != neighborhood.end(); ++it )
      {
      PointType neighbor = this->m_FixedTransformedPointSet->GetPoint( *it );
      typename PointType::VectorType neighborVector =
        neighbor.GetVectorFromOrigin();
      weightedPoint += ( neighborVector *
        measureValues[it - neighborhood.begin()] / measure );
      }
    }

  MeasureType distance = point.SquaredEuclideanDistanceTo( weightedPoint );
  MeasureType weight = preFactor *
    vcl_exp( -distance / ( 2.0 * vnl_math_sqr( this->m_Sigma ) ) ) /
    measure;
  typename PointType::VectorType force =
    ( weightedPoint - point ) * weight;

  for( unsigned int d = 0; d < localDerivative.Size(); d++ )
    {
    localDerivative[d] = force[d];
    }
}

template<class TFixedPointSet, class TMovingPointSet>
void
ExpectationBasedPointSetMetric<TFixedPointSet, TMovingPointSet>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Sigma: " << this->m_Sigma << std::endl;
}
} // end namespace itk


#endif
