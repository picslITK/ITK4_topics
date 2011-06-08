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
#ifndef __itkPointSetToPointSetMetric_txx
#define __itkPointSetToPointSetMetric_txx

#include "itkPointSetToPointSetMetric.h"

#include "itkIdentityTransform.h"
#include "itkVector.h"

namespace itk
{
/** Constructor */
template<class TFixedPointSet, class TMovingPointSet>
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::PointSetToPointSetMetric()
{
  this->m_FixedPointSet = NULL;    // has to be provided by the user.
  this->m_MovingPointSet = NULL;    // has to be provided by the user.

  this->m_FixedTransformedPointSet = NULL;
  this->m_MovingTransformedPointSet = NULL;

  this->m_CoordinateSystem = Fixed;

  // Set transforms to identity as default

  typedef IdentityTransform<CoordinateRepresentationType, FixedDimension>
    FixedIdentityTransformType;
  this->m_FixedTransform = FixedIdentityTransformType::New();
  this->m_FixedTransform->SetIdentity();

  this->m_FixedTransformModificationTime = this->m_FixedTransform->GetMTime();

  typedef IdentityTransform<CoordinateRepresentationType, MovingDimension>
    MovingIdentityTransformType;
  this->m_MovingIdentityTransform = MovingIdentityTransformType::New();
  this->m_MovingIdentityTransform->SetIdentity();

  this->m_MovingTransformModificationTime = this->m_MovingTransform->GetMTime();
}

/** Initialize the metric */
template<class TFixedPointSet, class TMovingPointSet>
void
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::Initialize( void )
throw ( ExceptionObject )
{
  if ( !this->m_FixedTransform )
    {
    itkExceptionMacro( "Fixed transform is not present" );
    }

  if ( !this->m_MovingTransform )
    {
    itkExceptionMacro( "Moving transform is not present" );
    }

  if ( !this->m_MovingPointSet )
    {
    itkExceptionMacro( "Moving point set is not present" );
    }

  if ( !this->m_FixedPointSet )
    {
    itkExceptionMacro( "Fixed point set is not present" );
    }

  // If the PointSet is provided by a source, update the source.
  if ( this->m_MovingPointSet->GetSource() )
    {
    this->m_MovingPointSet->GetSource()->Update();
    this->TransformMovingPointSet();
    }

  // If the point set is provided by a source, update the source.
  if ( this->m_FixedPointSet->GetSource() )
    {
    this->m_FixedPointSet->GetSource()->Update();
    this->TransformFixedPointSet();
    }
}

template<class TFixedPointSet, class TMovingPointSet>
unsigned int
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetNumberOfComponents()
{
  numberOfComponents = 0;
  if( this->m_CoordinateSystem == Fixed || this->m_CoordinateSystem == Both )
    {
    numberOfComponents += this->m_FixedPointSet->GetNumberOfPoints();
    }
  if( this->m_CoordinateSystem == Moving || this->m_CoordinateSystem == Both )
    {
    numberOfComponents += this->m_MovingPointSet->GetNumberOfPoints();
    }
  return numberOfComponents;
}

template<class TFixedPointSet, class TMovingPointSet>
typename PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetValue()
{
  this->TransformFixedPointSet();

  this->TransformMovingPointSet();

  MeasuresType measure = 0.0;

  if( this->m_CoordinateSystem == Fixed || this->m_CoordinateSystem == Both )
    {
    MovingPointSetIteratorType It = this->m_FixedTransformedPointSet->Begin();
    while( It != this->m_FixedTransformedPointSet->End() )
      {
      measure += this->GetLocalMovingNeighborhoodValue( It.Value() );
      ++It;
      }
    }
  if( this->m_CoordinateSystem == Moving || this->m_CoordinateSystem == Both )
    {
    FixedPointSetIteratorType It = this->m_MovingTransformedPointSet->Begin();
    while( It != this->m_MovingTransformedPointSet->End() )
      {
      measure += this->GetLocalFixedNeighborhoodValue( It.Value() );
      ++It;
      }
    }
  measure /= static_cast<MeasureType>( this->GetNumberOfComponents() );

  return measure;
}

template<class TFixedPointSet, class TMovingPointSet>
void
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetDerivative( DerivativeType &derivative )
{
  this->TransformFixedPointSet();

  this->TransformMovingPointSet();

  derivative.SetSize( this->GetNumberOfComponents(), vnl_math_max(
    FixedPointSetDimension, MovingPointSetDimension ) );
  derivative.Fill( 0 );

  unsigned long index = 0;

  if( this->m_CoordinateSystem == Fixed || this->m_CoordinateSystem == Both )
    {
    MovingPointSetIteratorType It = this->m_FixedTransformedPointSet->Begin();
    while( It != this->m_FixedTransformedPointSet->End() )
      {
      LocalMovingDerivativeType localDerivative =
        this->GetLocalMovingNeighborhoodDerivative( It.Value() );
      for( unsigned int d = 0; d < MovingPointSetDimension; d++ )
        {
        derivative(index, d) = localDerivative[d];
        }
      ++index;

      ++It;
      }
    }

  if( this->m_CoordinateSystem == Moving || this->m_CoordinateSystem == Both )
    {
    FixedPointSetIteratorType It = this->m_MovingTransformedPointSet->Begin();
    while( It != this->m_MovingTransformedPointSet->End() )
      {
      LocaFixedDerivativeType localDerivative =
        this->GetLocalFixedNeighborhoodDerivative( It.Value() );
      for( unsigned int d = 0; d < FixedPointSetDimension; d++ )
        {
        derivative(index, d) = localDerivative[d];
        }
      ++index

      ++It;
      }
    }
}

template<class TFixedPointSet, class TMovingPointSet>
void
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetValueAndDerivative( MeasureType &value, DerivativeType &derivative )
{
  this->TransformFixedPointSet();

  this->TransformMovingPointSet();

  derivative.SetSize( this->GetNumberOfComponents(), vnl_math_max(
    FixedPointSetDimension, MovingPointSetDimension ) );
  derivative.Fill( 0 );

  value = 0.0;

  unsigned long index = 0;

  if( this->m_CoordinateSystem == Fixed || this->m_CoordinateSystem == Both )
    {
    MovingPointSetIteratorType It = this->m_FixedTransformedPointSet->Begin();
    while( It != this->m_FixedTransformedPointSet->End() )
      {
      MeasureType localValue = 0.0;
      LocalMovingDerivativeType localDerivative;
      this->GetLocalMovingNeighborhoodValueAndDerivative( It.Value(), localValue,
        localDerivative );
      for( unsigned int d = 0; d < MovingPointSetDimension; d++ )
        {
        derivative(index, d) = localDerivative[d];
        }
      ++index;

      value += localValue;

      ++It;
      }
    }

  if( this->m_CoordinateSystem == Moving || this->m_CoordinateSystem == Both )
    {
    FixedPointSetIteratorType It = this->m_MovingTransformedPointSet->Begin();
    while( It != this->m_MovingTransformedPointSet->End() )
      {
      MeasureType localValue = 0.0;
      LocalFixedDerivativeType localDerivative;
      this->GetLocalFixedNeighborhoodValueAndDerivative( It.Value(), localValue,
        localDerivative );
      for( unsigned int d = 0; d < FixedPointSetDimension; d++ )
        {
        derivative(index, d) = localDerivative[d];
        }
      ++index;

      value += localValue;

      ++It;
      }
    }

  value /= static_cast<RealType>( index );
}

template<class TFixedPointSet, class TMovingPointSet>
void
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::TransformMovingPointSet()
{
  if( this->m_MovingTransform == NULL )
    {
    itkExceptionMacro( "No moving transform set to transform points." );
    }

  if( this->m_MovingTransform->GetMTime() >
    this->GetMTime() ||
    this->m_MovingTransformedPointSet == NULL )
    {
    this->m_MovingTransformedPointSet = MovingPointSetType::New();
    this->m_MovingTransformedPointSet->Initialize();

    MovingPointSetIteratorType It = this->m_MovingPointSet->Begin();
    while( It != this->m_MovingPointSet->End() )
      {
      FixedPointType point =
        this->m_MovingTransform->TransformPoint( It.Value() );

      this->m_MovingTransformedPointSet->InsertElement( point, It.Value() );

      ++It;
      }
    this->m_MovingTransformedPointSet->SetPointData(
      this->m_MovingPointSet->GetPointData() );
    }
}

template<class TFixedPointSet, class TMovingPointSet>
void
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::TransformFixedPointSet()
{
  if( this->m_FixedTransform == NULL )
    {
    itkExceptionMacro( "No fixed transform set to transform points." );
    }

  if( this->m_FixedTransform->GetMTime() >
    this->GetMTime() ||
    this->m_FixedTransformedPointSet == NULL )
    {
    this->m_FixedTransformedPointSet = FixedPointSetType::New();
    this->m_FixedTransformedPointSet->Initialize();

    FixedPointSetIteratorType It = this->m_FixedPointSet->Begin();
    while( It != this->m_FixedPointSet->End() )
      {
      MovingPointType point =
        this->m_FixedTransform->TransformPoint( It.Value() );

      this->m_FixedTransformedPointSet->InsertElement( point, It.Value() );

      ++It;
      }
    this->m_FixedTransformedPointSet->SetPointData(
      this->m_FixedPointSet->GetPointData() );
    }
}

/** PrintSelf */
template<class TFixedPointSet, class TMovingPointSet>
void
PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Fixed PointSet: " << this->m_FixedPointSet.GetPointer()
    << std::endl;
  os << indent << "Fixed Transform: " << this->m_FixedTransform.GetPointer()
    << std::endl;
  os << indent << "Moving PointSet: " << this->m_MovingPointSet.GetPointer()
    << std::endl;
  os << indent << "Moving Transform: " << this->m_MovingTransform.GetPointer()
    << std::endl;
}
} // end namespace itk

#endif
