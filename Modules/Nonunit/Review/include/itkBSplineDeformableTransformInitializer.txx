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
#ifndef __itkBSplineDeformableTransformInitializer_txx
#define __itkBSplineDeformableTransformInitializer_txx

#include "itkBSplineDeformableTransformInitializer.h"

#include "itkPointSet.h"
#include "itkBoundingBox.h"

namespace itk
{
template<class TTransform, class TImage>
BSplineDeformableTransformInitializer<TTransform, TImage>
::BSplineDeformableTransformInitializer() :
  m_Transform( NULL )
{
}

template<class TTransform, class TImage>
void
BSplineDeformableTransformInitializer<TTransform, TImage>
::InitializeTransform() const
{
  if( !this->m_Transform )
    {
    itkExceptionMacro( "Transform has not been set." );
    return;
    }
  if( !this->m_Image )
    {
    itkExceptionMacro( "Image has not been set." );
    return;
    }
  if( TImage::GetImageDimension() != SpaceDimension )
    {
    itkExceptionMacro( "Image dimensionality does not match the transform." );
    return;
    }

  OriginType                        transformDomainOrigin;
  PhysicalDimensionsType            transformDomainPhysicalDimensions;
  DirectionType                     transformDomainDirection;

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

  typedef typename ImageType::PointType              ImagePointType;
  typedef typename ImagePointType::CoordRepType      CoordRepType;

  typedef PointSet<CoordRepType, SpaceDimension>     PointSetType;
  typename PointSetType::Pointer cornerPoints = PointSetType::New();
  cornerPoints->Initialize();

  typedef typename PointSetType::PointType           PointType;
  typedef typename PointSetType::PointIdentifier     PointIdentifier;
  typedef typename PointType::RealType               RealType;
  typedef typename PointType::VectorType             VectorType;

  // We first convert the image corners into points which reside in physical
  // space and label them as indicated above.  We also store them using the
  // point set class which gives us easy access to the bounding box.

  IndexType startIndex = this->m_Image->GetRequestedRegion().GetIndex();
  for( unsigned int d = 0; d < vcl_pow( 2.0, SpaceDimension ); d++ )
    {
    IndexType whichIndex;
    for( unsigned int i = 0; i < SpaceDimension; i++ )
      {
      whichIndex[i] = startIndex[i] + static_cast<unsigned int>( ( d >> i ) & 1 ) *
        ( this->m_Image->GetRequestedRegion().GetSize()[i] - 1 );
      }
    ImagePointType point;
    this->m_Image->TransformIndexToPhysicalPoint( whichIndex, point );
    PointType corner;
    corner.CastFrom( point );
    cornerPoints->SetPoint( d, corner );
    }

  // We next determine which corner is the transform domain origin by which
  // point is closest to the minimum of the bounding box.

  typedef BoundingBox<unsigned int, SpaceDimension,
    typename PointSetType::CoordRepType,
    typename PointSetType::PointsContainer> BoundingBoxType;
  typename BoundingBoxType::Pointer bbox = BoundingBoxType::New();
  bbox->SetPoints( cornerPoints->GetPoints() );
  bbox->ComputeBoundingBox();

  transformDomainOrigin.Fill( 0 );
  PointIdentifier transformDomainOriginId = 0;
  RealType minDistance = NumericTraits<RealType>::max();

  for( unsigned int d = 0; d < cornerPoints->GetNumberOfPoints(); d++ )
    {
    PointType corner;
    cornerPoints->GetPoint( d, &corner );

    RealType distance = corner.SquaredEuclideanDistanceTo(
      bbox->GetMinimum() );
    if( distance < minDistance )
      {
      transformDomainOrigin.CastFrom( corner );
      minDistance = distance;
      transformDomainOriginId = static_cast<PointIdentifier>( d );
      }
    }

  // Now we need to find the transform direction matrix.  This is done
  // by using the domain origin and its adjacent neighbors to determine a new
  // rotated coordinate system.

  transformDomainDirection.SetIdentity();

  // We first determine which image axis is the most aligned with each physical axis.

  PointIdentifier minCornerId[SpaceDimension];
  double minAngle[SpaceDimension];

  for( unsigned int d = 0; d < SpaceDimension; d++ )
    {
    minAngle[d] = NumericTraits<double>::max();

    VectorType vectorAxis( 0.0 );
    vectorAxis[d] = 1.0;

    for( unsigned int i = 0; i < SpaceDimension; i++ )
      {
      PointIdentifier oppositeCornerId = ( 1 << i ) ^ transformDomainOriginId;

      PointType corner;
      cornerPoints->GetPoint( oppositeCornerId, &corner );

      VectorType vector = corner - transformDomainOrigin;
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

  // Now that we know which image axes corresponds to the unrotated coordinate
  // axes in physical space, we can easily construct the rotation matrix which
  // rotates a point from the unrotated coordinate system to the rotated
  // coordinate system.  This is done by placing the rotated axis vectors as
  // columns in the rotation matrix.

  for( unsigned int d = 0; d < SpaceDimension; d++ )
    {
    PointType corner;
    cornerPoints->GetPoint( minCornerId[d], &corner );

    VectorType vector = corner - transformDomainOrigin;

    // Note that specifying the size and spacing separately doesn't matter in
    // the case of the B-spline transform since the B-spline transform is a
    // continuous object over its finite domain.  However, their product matters
    // in that it defines the physical dimensions of the continuous B-spline
    // domain.  Here we specify an arbitrary size and force the spacing to
    // be a quantity such that ( ( size - 1 ) * spacing ) equals the dimensions
    // of the domain.

    transformDomainPhysicalDimensions[d] = vector.GetNorm();
    vector.Normalize();

    for( unsigned int i = 0; i < SpaceDimension; i++ )
      {
      transformDomainDirection[i][d] = vector[i];
      }
    }

  for( unsigned int d = 0; d < SpaceDimension; d++ )
    {
    PointType corner;
    cornerPoints->GetPoint( minCornerId[d], &corner );
    }

  this->m_Transform->SetTransformDomainOrigin( transformDomainOrigin );
  this->m_Transform->SetTransformDomainPhysicalDimensions(
    transformDomainPhysicalDimensions );
  this->m_Transform->SetTransformDomainDirection( transformDomainDirection );
}

template<class TTransform, class TImage>
void
BSplineDeformableTransformInitializer<TTransform, TImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Transform: " << std::endl;
  if ( this->m_Transform )
    {
    os << indent << this->m_Transform  << std::endl;
    }
  else
    {
    os << indent << "None" << std::endl;
    }
  os << indent << "Image: " << this->m_Image << std::endl;
}

}  // namespace itk

#endif
