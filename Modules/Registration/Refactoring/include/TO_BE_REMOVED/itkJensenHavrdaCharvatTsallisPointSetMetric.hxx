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
#ifndef __itkJensenHavrdaCharvatTsallisPointSetMetric_hxx
#define __itkJensenHavrdaCharvatTsallisPointSetMetric_hxx

#include "itkJensenHavrdaCharvatTsallisPointSetMetric.h"

namespace itk {

template<class TPointSet>
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>
::JensenHavrdaCharvatTsallisPointSetMetric() :
  m_UseRegularizationTerm( false ),
  m_UseAnisotropicCovariances( false ),
  m_PointSetSigma( 1.0 ),
  m_KernelSigma( 10.0 ),
  m_CovarianceKNeighborhood( 5.0 ),
  m_EvaluationKNeighborhood( 50 ),
  m_Alpha( 1 )
{
}

/** Initialize the metric */
template<class TPointSet>
void
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>
::Initialize( void ) throw ( ExceptionObject )
{
  if ( !this->m_FixedTransform )
    {
    itkExceptionMacro( "Fixed transform is not present" );
    }

  if ( !this->m_MovingTransform )
    {
    itkExceptionMacro( "Moving transform is not present" );
    }

  if ( !this->m_FixedPointSet )
    {
    itkExceptionMacro( "Fixed point set is not present" );
    }

  if ( !this->m_MovingPointSet )
    {
    itkExceptionMacro( "Moving point set is not present" );
    }

  // If the PointSet is provided by a source, update the source.
  if( this->m_MovingPointSet->GetSource() )
    {
    this->m_MovingPointSet->GetSource()->Update();
    }
  this->TransformMovingPointSet();

  // If the point set is provided by a source, update the source.
  if( this->m_FixedPointSet->GetSource() )
    {
    this->m_FixedPointSet->GetSource()->Update();
    }
  this->TransformFixedPointSet();

  // Initialize the fixed density function
  this->m_FixedDensityFunction = DensityFunctionType::New();
  this->m_FixedDensityFunction->SetKernelSigma( this->m_KernelSigma );
  this->m_FixedDensityFunction->SetRegularizationSigma(
    this->m_PointSetSigma );
  this->m_FixedDensityFunction->SetNormalize( true );
  this->m_FixedDensityFunction->SetUseAnisotropicCovariances(
    this->m_UseAnisotropicCovariances );
  this->m_FixedDensityFunction->SetCovarianceKNeighborhood(
    this->m_CovarianceKNeighborhood );
  this->m_FixedDensityFunction->SetEvaluationKNeighborhood(
    this->m_EvaluationKNeighborhood );
  this->m_FixedDensityFunction->SetInputPointSet(
    this->m_FixedTransformedPointSet );

  // Initialize the moving density function
  this->m_MovingDensityFunction = DensityFunctionType::New();
  this->m_MovingDensityFunction->SetKernelSigma( this->m_KernelSigma );
  this->m_MovingDensityFunction->SetRegularizationSigma(
    this->m_PointSetSigma );
  this->m_MovingDensityFunction->SetNormalize( true );
  this->m_MovingDensityFunction->SetUseAnisotropicCovariances(
    this->m_UseAnisotropicCovariances );
  this->m_MovingDensityFunction->SetCovarianceKNeighborhood(
    this->m_CovarianceKNeighborhood );
  this->m_MovingDensityFunction->SetEvaluationKNeighborhood(
    this->m_EvaluationKNeighborhood );
  this->m_MovingDensityFunction->SetInputPointSet(
    this->m_MovingTransformedPointSet );
}

/** Get the match Measure */
template<class TPointSet>
typename JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>
::MeasureType
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>
::GetValue() const
{
  DensityFunctionPointer densityFunctions[2];
  densityFunctions[0] = this->m_MovingDensityFunction;
  densityFunctions[1] = this->m_FixedDensityFunction;

  RealType totalNumberOfPoints = 0;
  for( unsigned int d = 0; d < 2; d++ )
    {
    totalNumberOfPoints +=
      static_cast<RealType>( densityFunctions[d]->GetInputPointSet()->
      GetNumberOfPoints() );
    }
  RealType totalNumberOfSamples = totalNumberOfPoints;

  MeasureType measure = 0.0;

  MeasureType energyTerm1 = 0.0;
  MeasureType energyTerm2 = 0.0;

  /**
    * first term
    */
  RealType prefactor = -1.0 / totalNumberOfSamples;
  if( this->m_Alpha != 1.0 )
    {
    prefactor /= ( this->m_Alpha - 1.0 );
    }
  for( unsigned int d = 0; d < 2; d++ )
    {
    PointsContainerConstIterator It =
      densityFunctions[d]->GetInputPointSet()->GetPoints()->Begin();
    while( It != densityFunctions[d]->GetInputPointSet()->GetPoints()->End() )
      {
      PointType samplePoint = It.Value();

      RealType probabilityStar =
        densityFunctions[1-d]->Evaluate( samplePoint ) *
        static_cast<RealType>( densityFunctions[1-d]->GetInputPointSet()->
        GetNumberOfPoints() );
      probabilityStar /= totalNumberOfPoints;

      if( probabilityStar == 0 )
        {
        ++It;
        continue;
        }

      if( this->m_Alpha == 1.0 )
        {
        energyTerm1 += vcl_log( probabilityStar );
        }
      else
        {
        energyTerm1 += vcl_pow( probabilityStar,
          static_cast<RealType>( this->m_Alpha - 1.0 ) );
        }
      ++It;
      }
    if( this->m_Alpha != 1.0 )
      {
      energyTerm1 -= 1.0;
      }
    }
  energyTerm1 *= prefactor;

  /**
    * second term, i.e. regularization term
    */
  if( this->m_UseRegularizationTerm )
    {
    for( unsigned int d = 0; d < 2; d++ )
      {
      RealType prefactor2 = -static_cast<RealType>(
        densityFunctions[d]->GetInputPointSet()->GetNumberOfPoints() ) /
        ( totalNumberOfPoints * static_cast<RealType>(
        densityFunctions[d]->GetInputPointSet()->GetNumberOfPoints() ) );
      if( this->m_Alpha != 1.0 )
        {
        prefactor2 /= ( this->m_Alpha - 1.0 );
        }
      PointsContainerConstIterator It =
        densityFunctions[d]->GetInputPointSet()->GetPoints()->Begin();
      while( It != densityFunctions[d]->GetInputPointSet()->GetPoints()->End() )
        {
        PointType samplePoint = It.Value();

        RealType probability = densityFunctions[d]->Evaluate( samplePoint );

        if( probability == 0 )
          {
          ++It;
          continue;
          }

        if( this->m_Alpha == 1.0 )
          {
          energyTerm2 += ( prefactor2 * vcl_log( probability ) );
          }
        else
          {
          energyTerm2 += ( prefactor2 * vcl_pow( probability,
            static_cast<RealType>( this->m_Alpha - 1.0 ) ) );
          }
        ++It;
        }
      if( this->m_Alpha != 1.0 )
        {
        energyTerm2 -= 1.0;
        }
      energyTerm2 *= prefactor2;
      }
    }

  measure = energyTerm1 - energyTerm2;

  return measure;
}

/** Get the Derivative Measure */
template<class TPointSet>
void
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>
::GetDerivative( DerivativeType &derivative ) const
{
  MeasureType value;
  this->GetValueAndDerivative( value, derivative );
}

/** Get both the match Measure and the Derivative Measure  */
template<class TPointSet>
void
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>
::GetValueAndDerivative( MeasureType &value, DerivativeType &derivative ) const
{
  DensityFunctionPointer densityFunctions[2];
  densityFunctions[0] = this->m_FixedDensityFunction;
  densityFunctions[1] = this->m_MovingDensityFunction;

  unsigned int start = 0;
  unsigned int end = 0;
  if( this->GetDerivativeSource() == Superclass::Moving )
    {
    start = 1;
    end = 1;
    }
  else if( this->GetDerivativeSource() == Superclass::Fixed )
    {
    start = 0;
    end = 0;
    }
  else
    {
    start = 0;
    end = 1;
    }

  RealType totalNumberOfPoints = 0;
  for( unsigned int d = 0; d < 2; d++ )
    {
    totalNumberOfPoints +=
      static_cast<RealType>( densityFunctions[d]->GetInputPointSet()->
      GetNumberOfPoints() );
    }
  RealType totalNumberOfSamples = totalNumberOfPoints;

  derivative.SetSize( this->GetNumberOfComponents() * PointDimension );
  derivative.Fill( 0 );

  value = 0;

  /**
   * first term
   */
  RealType energyTerm1 = 0.0;
  RealType energyTerm2 = 0.0;

  RealType prefactor[2];
  prefactor[0] = -1.0 / totalNumberOfSamples;
  if( this->m_Alpha != 1.0 )
    {
    prefactor[0] /= ( this->m_Alpha - 1.0 );
    }
  prefactor[1] = 1.0 / ( totalNumberOfSamples * totalNumberOfPoints );

  for( unsigned int d = start; d <= end; d++ )
    {
    PointsContainerConstIterator It =
      densityFunctions[d]->GetInputPointSet()->GetPoints()->Begin();
    while( It != densityFunctions[d]->GetInputPointSet()->GetPoints()->End() )
      {
      PointType samplePoint = It.Value();

      RealType probabilityStar =
          densityFunctions[1-d]->Evaluate( samplePoint ) *
          static_cast<RealType>( densityFunctions[1-d]->GetInputPointSet()->
          GetNumberOfPoints() );
      probabilityStar /= totalNumberOfPoints;

      if( probabilityStar == 0 )
        {
        ++It;
        continue;
        }

      if( this->m_Alpha == 1.0 )
        {
        energyTerm1 += ( prefactor[0] * vcl_log( probabilityStar ) );
        }
      else
        {
        energyTerm1 += ( prefactor[0] * vcl_pow( probabilityStar,
          static_cast<RealType>( this->m_Alpha - 1.0 ) ) );
        }

      RealType probabilityStarFactor = vcl_pow( probabilityStar,
        static_cast<RealType>( 2.0 - this->m_Alpha ) );

      typename DensityFunctionType::NeighborsIdentifierType neighbors;
      densityFunctions[1-d]->GetPointsLocator()->FindClosestNPoints(
        samplePoint, this->m_EvaluationKNeighborhood, neighbors );

      for( unsigned int n = 0; n < neighbors.size(); n++ )
        {
        RealType gaussian = densityFunctions[1-d]->
          GetGaussian( neighbors[n] )->Evaluate( samplePoint );

        if( gaussian == 0 )
          {
          continue;
          }

        typename GaussianType::MeanType mean =
          densityFunctions[1-d]->GetGaussian( neighbors[n] )->GetMean();


        Array<CoordRepType> diffMean( PointDimension );
        for( unsigned int i = 0; i < PointDimension; i++ )
          {
          diffMean[i] = mean[i] - samplePoint[i];
          }

        if( this->m_UseAnisotropicCovariances )
          {
          typename GaussianType::CovarianceType Ci =
            densityFunctions[1-d]->GetGaussian( neighbors[n] )->
            GetInverseCovariance();
          diffMean = Ci * diffMean;
          }
        else
          {
          diffMean /= densityFunctions[1-d]->GetGaussian( neighbors[n] )->
            GetCovariance()(0, 0);
          }

        for( unsigned int i = 0; i < PointDimension; i++ )
          {
          derivative( ( end - start ) * densityFunctions[0]->GetInputPointSet()->
            GetNumberOfPoints() * PointDimension +
            It.Index() * PointDimension + i ) -=
            diffMean[i] * ( prefactor[1] * gaussian / probabilityStarFactor );
          }
        }
      ++It;
      }
    }
  if( this->m_Alpha != 1.0 )
    {
    energyTerm1 -= 1.0;
    }
  energyTerm1 *= prefactor[0];

  /**
   * second term, i.e. regularization term
   */
  if( this->m_UseRegularizationTerm )
    {
    for( unsigned int d = start; d <= end; d++ )
      {
      RealType prefactor2 = -1.0 / ( static_cast<RealType>(
        densityFunctions[1-d]->GetInputPointSet()->GetNumberOfPoints() ) *
        totalNumberOfPoints );

      typename PointSetType::PointsContainerConstIterator It
        = densityFunctions[1-d]->GetInputPointSet()->GetPoints()->Begin();
      while( It != densityFunctions[1-d]->GetInputPointSet()->GetPoints()->End() )
        {
        PointType samplePoint = It.Value();

        RealType probability = densityFunctions[1-d]->Evaluate( samplePoint );

        if( probability == 0 )
          {
          ++It;
          continue;
          }

        if( this->m_Alpha == 1.0 )
          {
          energyTerm2 += ( prefactor2 * vcl_log( probability ) );
          }
        else
          {
          energyTerm2 += ( prefactor2 * vcl_pow( probability,
            static_cast<RealType>( this->m_Alpha - 1.0 ) ) );
          }

        RealType probabilityFactor = vcl_pow( probability,
          static_cast<RealType>( 2.0 - this->m_Alpha ) );
        probabilityFactor *= ( densityFunctions[1-d]->GetInputPointSet()->
          GetNumberOfPoints() / totalNumberOfSamples );

        typename DensityFunctionType::NeighborsIdentifierType neighbors;
        densityFunctions[1-d]->GetPointsLocator()->FindClosestNPoints(
          samplePoint, this->m_EvaluationKNeighborhood, neighbors );

        for( unsigned int n = 0; n < neighbors.size(); n++ )
          {
          RealType gaussian = densityFunctions[1-d]->
            GetGaussian( neighbors[n] )->Evaluate( samplePoint );
          if( gaussian == 0 )
            {
            continue;
            }

          typename GaussianType::MeanType mean
            = densityFunctions[1-d]->GetGaussian( neighbors[n] )->GetMean();

          Array<CoordRepType> diffMean( PointDimension );
          for( unsigned int i = 0; i < PointDimension; i++ )
            {
            diffMean[i] = mean[i] - samplePoint[i];
            }

          if( this->m_UseAnisotropicCovariances )
            {
            typename GaussianType::CovarianceType Ci =
              densityFunctions[1-d]->GetGaussian( neighbors[n] )->
              GetInverseCovariance();
            diffMean = Ci * diffMean;
            }
          else
            {
            diffMean /= densityFunctions[1-d]->GetGaussian( neighbors[n] )->
              GetCovariance()(0, 0);
            }

          for( unsigned int i = 0; i < PointDimension; i++ )
            {
            derivative( ( end - start ) * densityFunctions[0]->
              GetInputPointSet()->GetNumberOfPoints() * PointDimension +
              It.Index() * PointDimension + i ) +=
              diffMean[i] * ( prefactor2 * gaussian / probabilityFactor );
            }
          }
        ++It;
        }
      }
    if( this->m_Alpha != 1.0 )
      {
      energyTerm2 -= 1.0;
      }
    energyTerm2 *= prefactor[1];
    }

  derivative *= -1.0;
  value = energyTerm1 - energyTerm2;
}


template<class TPointSet>
void
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Use regularization term: "
     << this->m_UseRegularizationTerm << std::endl;
  os << indent << "Alpha: "
     << this->m_Alpha << std::endl;

  os << indent << "Point set sigma: "
     << this->m_PointSetSigma << std::endl;

  if( this->m_UseAnisotropicCovariances )
    {
    os << indent << "Kernel sigma: "
       << this->m_KernelSigma << std::endl;
    os << indent << "Covariance k-neighborhood: "
       << this->m_CovarianceKNeighborhood << std::endl;
    }
  else
    {
    os << indent << "Isotropic covariances are used." << std::endl;
    }

}
} // end namespace itk


#endif
