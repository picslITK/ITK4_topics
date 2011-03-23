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
#ifndef __itkTransformVirtualDomainCalculator_txx
#define __itkTransformVirtualDomainCalculator_txx

#include "itkTransformVirtualDomainCalculator.h"

#include "itkCastImageFilter.h"
#include "itkPointSet.h"

namespace itk
{
template<class TInputImage1, class TInputImage2>
TransformVirtualDomainCalculator<TInputImage1, TInputImage2>
::TransformVirtualDomainCalculator()
{
  this->m_InputImage1 = NULL;
  this->m_InputImage2 = NULL;

  this->m_VirtualDomainOrigin.Fill( 0.0 );
  this->m_VirtualDomainDirection.SetIdentity();

  this->m_PhysicallyConsistentOrigin1.Fill( 0.0 );
  this->m_PhysicallyConsistentOrigin2.Fill( 0.0 );

  this->m_PhysicallyConsistentDirection1.SetIdentity();
  this->m_PhysicallyConsistentDirection2.SetIdentity();

  this->m_UsePhysicalConsistency = true;
}

template<class TInputImage1, class TInputImage2>
void
TransformVirtualDomainCalculator<TInputImage1, TInputImage2>
::CalculateVirtualDomainParameters()
{
  OriginType origin1 = this->m_InputImage1->GetOrigin();
  OriginType origin2 = this->m_InputImage2->GetOrigin();

  DirectionType direction1 = this->m_InputImage1->GetDirection();
  DirectionType direction2 = this->m_InputImage2->GetDirection();

  if( this->m_UsePhysicalConsistency )
    {
    this->CalculatePhysicallyConsistentDomainParameters( this->m_InputImage1,
      origin1, direction1 );
    this->m_PhysicallyConsistentOrigin1 = origin1;
    this->m_PhysicallyConsistentDirection1 = direction1;

    typedef CastImageFilter<Image2Type, Image1Type> CasterType;
    typename CasterType::Pointer caster = CasterType::New();
    caster->SetInput( this->m_InputImage2 );
    caster->Update();

    this->CalculatePhysicallyConsistentDomainParameters( caster->GetOutput(),
      origin2, direction2 );

    this->m_PhysicallyConsistentOrigin2 = origin2;
    this->m_PhysicallyConsistentDirection2 = direction2;
    }

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    this->m_VirtualDomainOrigin[d] = 0.5 * ( origin1[d] + origin2[d] );
    }
  this->m_VirtualDomainDirection = this->CalculateUnbiasedDirectionCosineMatrix(
    direction1, direction2 );
}

template<class TInputImage1, class TInputImage2>
void
TransformVirtualDomainCalculator<TInputImage1, TInputImage2>
::CalculatePhysicallyConsistentDomainParameters( Image1Pointer image,
  OriginType &physicalDomainOrigin, DirectionType &physicalDomainDirection )
{
  if( !this->m_InputImage1 || !this->m_InputImage2 )
    {
    itkExceptionMacro( "One or both input images have not been set." );
    }

  // Determine the image corners.  We keep track of the relative location of
  // the corners using a binary labeling system.  For example, in a 3-D coordinate
  // system aligned with the x,y,z axes, we have 8 points labeled as follows:
  //
  //  1. 000  min_x, min_y, min_z
  //  2. 001  min_x, min_y, max_z
  //  3. 010  min_x, max_y, min_z
  //  4. 011  min_x, max_y, max_z
  //  5. 100  max_x, min_y, min_z
  //  6. 101  max_x, min_y, max_z
  //  7. 110  max_x, max_y, min_z
  //  8. 111  max_x, max_y, max_z
  //
  // We use this binary description of the corners in n-dimensions because it
  // allows us to know the adjacent neighbors of an arbitrary image corner. For
  // example, suppose we locate the transform domain origin at the corner 011
  // the adjacent neighbors which form the rotated coordinate system are
  // 111, 001, and 010.  Notice that we just change 1 bit at a time from the
  // origin to determine these axes.  Thus bitwise operators are used
  // throughout the code so that the initializer is generalized to n-dimensions.

  typedef typename Image1Type::PointType             ImagePointType;
  typedef typename ImagePointType::CoordRepType      CoordRepType;

  typedef PointSet<CoordRepType, ImageDimension>     PointSetType;
  typename PointSetType::Pointer cornerPoints = PointSetType::New();
  cornerPoints->Initialize();

  typedef typename PointSetType::PointType           PointType;
  typedef typename PointSetType::PointIdentifier     PointIdentifier;
  typedef typename PointType::RealType               RealType;
  typedef typename PointType::VectorType             VectorType;

  // We first convert the image corners into points which reside in physical
  // space and label them as indicated above.  We also store them using the
  // point set class which gives us easy access to the bounding box.

  IndexType startIndex = image->GetRequestedRegion().GetIndex();
  for( unsigned int d = 0; d < vcl_pow( 2.0, ImageDimension ); d++ )
    {
    IndexType whichIndex;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      whichIndex[i] = startIndex[i] + static_cast<unsigned int>( ( d >> i ) & 1 ) *
        ( image->GetRequestedRegion().GetSize()[i] - 1 );
      }
    ImagePointType point;
    image->TransformIndexToPhysicalPoint( whichIndex, point );
    PointType corner;
    corner.CastFrom( point );
    cornerPoints->SetPoint( d, corner );
    }

  // We next determine which corner is the transform domain origin by which
  // point is closest to the minimum of the bounding box.

  physicalDomainOrigin.Fill( 0 );
  PointIdentifier physicalDomainOriginId = 0;
  RealType minDistance = NumericTraits<RealType>::max();

  for( unsigned int d = 0; d < cornerPoints->GetNumberOfPoints(); d++ )
    {
    PointType corner;
    cornerPoints->GetPoint( d, &corner );

    RealType distance = corner.SquaredEuclideanDistanceTo(
      cornerPoints->GetBoundingBox()->GetMinimum() );
    if( distance < minDistance )
      {
      physicalDomainOrigin.CastFrom( corner );
      minDistance = distance;
      physicalDomainOriginId = static_cast<PointIdentifier>( d );
      }
    }

  // Now we need to find the transform direction matrix.  This is done
  // by using the domain origin and its adjacent neighbors to determine a new
  // rotated coordinate system.

  physicalDomainDirection.SetIdentity();

  // We first determine which image axis is the most aligned with each physical axis.

  PointIdentifier minCornerId[ImageDimension];
  double minAngle[ImageDimension];

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    minAngle[d] = NumericTraits<double>::max();

    VectorType vectorAxis( 0.0 );
    vectorAxis[d] = 1.0;

    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      PointIdentifier oppositeCornerId = ( 1 << i ) ^ physicalDomainOriginId;

      PointType corner;
      cornerPoints->GetPoint( oppositeCornerId, &corner );

      VectorType vector = corner - physicalDomainOrigin;
      vector.Normalize();

      double theta = angle( vectorAxis.GetVnlVector(), vector.GetVnlVector() );

      if( theta < minAngle[d] )
        {
        bool alreadyFound = false;
        for( unsigned int j = 0; j < d; j++ )
          {
          if( minCornerId[j] == oppositeCornerId )
            {
            alreadyFound = true;
            break;
            }
          }
        if( !alreadyFound )
          {
          minCornerId[d] = oppositeCornerId;
          minAngle[d] = theta;
          }
        }
      }
    }

  // Now that we know which image axes correspond to the unrotated coordinate
  // axes in physical space, we can easily construct the rotation matrix which
  // rotates a point from the unrotated coordinate system to the rotated
  // coordinate system.  This is done by placing the rotated axis vectors as
  // columns in the rotation matrix.

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    PointType corner;
    cornerPoints->GetPoint( minCornerId[d], &corner );

    VectorType vector = corner - physicalDomainOrigin;
    vector.Normalize();

    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      physicalDomainDirection[i][d] = vector[i];
      }
    }
}

template<class TInputImage1, class TInputImage2>
typename TransformVirtualDomainCalculator<TInputImage1, TInputImage2>
::DirectionType
TransformVirtualDomainCalculator<TInputImage1, TInputImage2>
::CalculateUnbiasedDirectionCosineMatrix( DirectionType direction1,
  DirectionType direction2 )
{
  DirectionType Aold = direction1;
  DirectionType Bold = DirectionType( direction2.GetInverse() );

  DirectionType A;
  DirectionType B;

  float epsilon = 1.0;
  unsigned int iterations = 0;
  while( epsilon > 1.0e-6 && iterations++ < 1000 )
    {
    DirectionType BoldInverse = DirectionType( Bold.GetInverse() );
    DirectionType AoldInverse = DirectionType( Aold.GetInverse() );

    A = ( Aold + BoldInverse ) / 2.0;
    B = ( Bold + AoldInverse ) / 2.0;

    epsilon = ( ( A - Aold ).GetVnlMatrix() ).frobenius_norm();

    Aold = A;
    Bold = B;
    }
  if( epsilon > 1.0e6 && iterations >= 1000 )
    {
    itkExceptionMacro( "Unbiased cosine matrix estimation failed to converge." );
    }

  return A;
}

template<class TInputImage1, class TInputImage2>
void
TransformVirtualDomainCalculator<TInputImage1, TInputImage2>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << "Virtual domain origin: " << this->m_VirtualDomainOrigin
    << std::endl;
  os << "Virtual domain direction: " << std::endl
    << this->m_VirtualDomainDirection << std::flush;
  if( this->m_UsePhysicalConsistency )
    {
    os << "Use physical consistency." << std::endl;
    os << indent << "Input image 1 physically consistent parameters:"
      << std::endl;
    os << indent << "   origin = "
      << this->m_PhysicallyConsistentOrigin1 << std::endl;
    os << indent << "   direction = "
      << this->m_PhysicallyConsistentDirection1 << std::endl;
    os << indent << "Input image 2 physically consistent parameters:"
      << std::endl;
    os << indent << "   origin = "
      << this->m_PhysicallyConsistentOrigin2 << std::endl;
    os << indent << "   direction = "
      << this->m_PhysicallyConsistentDirection2 << std::endl;
    }
  else
    {
    os << "Do not use physical consistency." << std::endl;
    }
}
}  // namespace itk

#endif
