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
#ifndef __itkBSplineDeformationFieldTransform_hxx
#define __itkBSplineDeformationFieldTransform_hxx

#include "itkBSplineDeformationFieldTransform.h"

#include "itkContinuousIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

/**
 * Constructor
 */
template<class TScalar, unsigned int NDimensions>
BSplineDeformationFieldTransform<TScalar, NDimensions>::
BSplineDeformationFieldTransform() : Superclass()
{
  this->m_SplineOrder = 3;
  this->m_NumberOfFittingLevels.Fill( 3 );
  this->m_NumberOfControlPoints.Fill( 4 );

  this->m_CalculateApproximateInverseDeformationField = false;
}

/**
 * Destructor
 */
template<class TScalar, unsigned int NDimensions>
BSplineDeformationFieldTransform<TScalar, NDimensions>::
~BSplineDeformationFieldTransform()
{
}

/**
 * set deformation field and project it onto the space of b-spline transforms
 */
template<class TScalar, unsigned int NDimensions>
void
BSplineDeformationFieldTransform<TScalar, NDimensions>
::SetDeformationField( DeformationFieldType *deformationField )
{
  typename PointSetType::Pointer fieldPoints = PointSetType::New();
  fieldPoints->Initialize();

  ImageRegionConstIteratorWithIndex<DeformationFieldType>
    It( deformationField, deformationField->GetRequestedRegion() );

  itkDebugMacro( "Extracting points from input deformation field. " )

  // Temporarily set the direction cosine to identity since the B-spline
  // approximation algorithm works in parametric space and not physical
  // space.

  typename DeformationFieldType::DirectionType identity;
  identity.SetIdentity();

  typename DeformationFieldType::DirectionType originalDirection =
    deformationField->GetDirection();

  deformationField->SetDirection( identity );

  unsigned int N = 0;
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    DisplacementVectorType data = It.Get();

    typename PointSetType::PointType point;
    deformationField->TransformIndexToPhysicalPoint( It.GetIndex(), point );

    fieldPoints->SetPointData( N, data );
    fieldPoints->SetPoint( N, point );
    N++;
    }
  deformationField->SetDirection( originalDirection );

  itkDebugMacro( "Calculating the B-spline deformation field. " );

  ArrayType close;
  close.Fill( false );

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetOrigin( deformationField->GetOrigin() );
  bspliner->SetSpacing( deformationField->GetSpacing() );
  bspliner->SetSize( deformationField->GetLargestPossibleRegion().GetSize() );
  bspliner->SetDirection( deformationField->GetDirection() );
  bspliner->SetNumberOfLevels( this->m_NumberOfFittingLevels );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetNumberOfControlPoints( this->m_NumberOfControlPoints );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( fieldPoints );
  bspliner->SetGenerateOutputImage( true );
  bspliner->Update();

  this->m_DeformationFieldControlPointLattice = bspliner->GetPhiLattice();

  typename DeformationFieldType::Pointer bsplineDeformationField =
    bspliner->GetOutput();
  bsplineDeformationField->DisconnectPipeline();

  Superclass::SetDeformationField( bsplineDeformationField );

  if( this->m_CalculateApproximateInverseDeformationField )
    {
    typename PointSetType::Pointer inverseFieldPoints = PointSetType::New();
    inverseFieldPoints->Initialize();

    typedef typename DeformationFieldType::PointType DeformationFieldPointType;
    typedef typename DeformationFieldPointType::CoordRepType CoordRepType;

    deformationField->SetDirection( identity );

    N = 0;
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      DisplacementVectorType data = It.Get();

      typename PointSetType::PointType point;
      deformationField->TransformIndexToPhysicalPoint( It.GetIndex(), point );
      DeformationFieldPointType inversePoint = point + data;

      ContinuousIndex<CoordRepType> cidx;
      deformationField->TransformPhysicalPointToContinuousIndex(
        inversePoint, cidx );
      if( !deformationField->GetRequestedRegion().IsInside( cidx ) )
        {
        continue;
        }

      inverseFieldPoints->SetPointData( N, -data );
      inverseFieldPoints->SetPoint( N, inversePoint );
      N++;
      }
    deformationField->SetDirection( originalDirection );

    itkDebugMacro( "Calculating the inverse approximation of the "
      << "B-spline deformation field. " );

    typename BSplineFilterType::Pointer inverseBSpliner = BSplineFilterType::New();
    inverseBSpliner->SetOrigin( deformationField->GetOrigin() );
    inverseBSpliner->SetSpacing( deformationField->GetSpacing() );
    inverseBSpliner->SetSize( deformationField->GetLargestPossibleRegion().GetSize() );
    inverseBSpliner->SetDirection( deformationField->GetDirection() );
    inverseBSpliner->SetNumberOfLevels( this->m_NumberOfFittingLevels );
    inverseBSpliner->SetSplineOrder( this->m_SplineOrder );
    inverseBSpliner->SetNumberOfControlPoints( this->m_NumberOfControlPoints );
    inverseBSpliner->SetCloseDimension( close );
    inverseBSpliner->SetInput( inverseFieldPoints );
    inverseBSpliner->SetGenerateOutputImage( true );
    inverseBSpliner->Update();

    this->m_InverseDeformationFieldControlPointLattice =
      inverseBSpliner->GetPhiLattice();

    typename DeformationFieldType::Pointer inverseBSplineDeformationField =
      inverseBSpliner->GetOutput();
    inverseBSplineDeformationField->DisconnectPipeline();

    Superclass::SetInverseDeformationField( inverseBSplineDeformationField );
    }
}

template <class TScalar, unsigned int NDimensions>
void
BSplineDeformationFieldTransform<TScalar, NDimensions>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os,indent );

  std::cout << indent << "B-spline parameters: " << std::endl;
  std::cout << indent << "  spline order = " << this->m_SplineOrder << std::endl;
  std::cout << indent << "  number of control points = "
    << this->m_NumberOfControlPoints << std::endl;
  std::cout << indent << "  number of fitting levels = "
    << this->m_NumberOfFittingLevels << std::endl;
}
} // namespace itk

#endif
