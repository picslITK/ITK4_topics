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
#ifndef __itkOptimizerParameterEstimator_txx
#define __itkOptimizerParameterEstimator_txx

#include "itkOptimizerParameterEstimator.h"

namespace itk
{

template< class TFixedImage, class TMovingImage, class TTransform >
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::OptimizerParameterEstimator()
{
  // Euclidean distance
  m_LNorm = 2;
  // true means to transform from fixed domain to moving domain
  m_TransformForward = true;

  m_ScaleStrategy = ScalesFromShift;

}

/** Compute parameter scales */
template< class TFixedImage, class TMovingImage, class TTransform >
void
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::EstimateScales(ParametersType parameters,
                 ScalesType &parameterScales)
{

  switch (m_ScaleStrategy)
    {
    case ScalesFromShift:
      {
      this->EstimateScalesFromMaximumShift(parameters, parameterScales);
      break;
      }
    case ScalesFromJacobian:
      {
      this->EstimateScalesFromJacobian(parameters, parameterScales);
      break;
      }
    default:
      {
      itkExceptionMacro(" Undefined strategy to decide scales.");
      }
    }
}

/** Compute parameter scales */
template< class TFixedImage, class TMovingImage, class TTransform >
void
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::EstimateScalesFromMaximumShift(ParametersType parameters,
                                 ScalesType &parameterScales)
{
  TransformPointer transform = TransformType::New();
  transform->SetParameters(parameters);
  this->SetTransform(transform);

  double maxShift;
  unsigned int numPara = parameters.size();

  if (parameterScales.size() != numPara)
    {
    itkExceptionMacro(" The size of scales does not match that of transform parameters.");
    return;
    }

  this->SampleImageDomain();

  ParametersType deltaParameters(numPara);
  for (unsigned int j=0; j<numPara; j++)
    {
    deltaParameters[j] = 0;
    }

  // To avoid division-by-zero, assign a small value for parameters
  // whose impact on voxel shift is zero. This small value may be
  // the minimum non-zero shift.
  double minNonZeroShift = NumericTraits<double>::max();

  // compute variation for each transform parameter
  for (unsigned int i=0; i<numPara; i++)
    {
    deltaParameters[i] = 1;
    maxShift = this->ComputeMaximumVoxelShift(parameters, deltaParameters);
    deltaParameters[i] = 0;

    parameterScales[i] = maxShift;
    if ( maxShift != 0 && maxShift < minNonZeroShift )
      {
      minNonZeroShift = maxShift;
      }
    }

  for (unsigned int i=0; i<numPara; i++)
    {
    if (parameterScales[i] == 0)
      {
      parameterScales[i] = minNonZeroShift;
      }
    parameterScales[i] *= parameterScales[i];
    }

}

/** Compute parameter scales */
template< class TFixedImage, class TMovingImage, class TTransform >
void
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::EstimateScalesFromJacobian(ParametersType parameters,
                             ScalesType &parameterScales)
{
  TransformPointer transform = TransformType::New();
  transform->SetParameters(parameters);
  this->SetTransform(transform);

  this->SampleImageDomain();

  unsigned int numPara = parameterScales.size();
  double normSquare;
  double *norms = new double[numPara];
  unsigned int numSamples = 0;

  numSamples = m_ImageSamples.size();

  PointType point;
  typename TransformType::JacobianType jacobian;

  for (unsigned int p=0; p<numPara; p++)
    {
    norms[p] = 0;
    }

  // find max shift by checking each sample point
  for (unsigned int c=0; c<numSamples; c++)
    {
    point = m_ImageSamples[c];
    jacobian = transform->GetJacobian(point);

    for (unsigned int p=0; p<numPara; p++)
      {
      normSquare = 0;
      for (unsigned int d=0; d<ImageDimension; d++)
        {
        normSquare = normSquare + jacobian[d][p] * jacobian[d][p];
        }
      norms[p] = norms[p] + normSquare;
      }
    } //for numSamples

  for (unsigned int p=0; p<numPara; p++)
    {
    parameterScales[p] = norms[p] / numSamples;
    }

}

/** Set the sample points for computing pixel shifts */
template< class TFixedImage, class TMovingImage, class TTransform >
void
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::SampleImageDomain()
{
  switch (m_ScaleStrategy)
    {
    case ScalesFromShift:
      {
      if (m_TransformForward)
        {
        this->SampleWithCornerPoints<TFixedImage>();
        }
      else
        {
        this->SampleWithCornerPoints<TMovingImage>();
        }
      break;
      }
    case ScalesFromJacobian:
      {
      if (m_TransformForward)
        {
        this->SampleImageDomainRandomly<TFixedImage>();
        }
      else
        {
        this->SampleImageDomainRandomly<TMovingImage>();
        }
      break;
      }
    default:
      {
      itkExceptionMacro(" Undefined strategy to decide scales.");
      }
    }
}

/**
 * Sample the physical coordinates of image in uniform random
 */
template< class TFixedImage, class TMovingImage, class TTransform >
template< class TImage>
void
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::SampleImageDomainRandomly()
{
  typename TImage::ConstPointer image;

  if ( typeid(TImage) == typeid(TFixedImage) ) { image = m_FixedImage; }
  else if ( typeid(TImage) == typeid(TMovingImage) ) { image = m_MovingImage; }
  else itkExceptionMacro(" TImage must be either TFixedImage or TMovingImage.");

  int numSamples = 1000;
  m_ImageSamples.clear();

  // Set up a random interator within the user specified fixed image region.
  typedef ImageRandomConstIteratorWithIndex<TImage> RandomIterator;
  RandomIterator randIter( image, image->GetLargestPossibleRegion() );

  randIter.SetNumberOfSamples( numSamples );
  randIter.GoToBegin();

  PointType point;
  IndexType index;

  for (int i=0; i<numSamples; i++)
    {
    index = randIter.GetIndex();
    image->TransformIndexToPhysicalPoint( index, point );
    m_ImageSamples.push_back(point);
    ++randIter;
    }
}

/**
 * Get the physical coordinates of image corners
 */
template< class TFixedImage, class TMovingImage, class TTransform >
template< class TImage >
void
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::SampleWithCornerPoints()
{
  typename TImage::ConstPointer image;

  if ( typeid(TImage) == typeid(TFixedImage) ) { image = m_FixedImage; }
  else if ( typeid(TImage) == typeid(TMovingImage) ) { image = m_MovingImage; }
  else itkExceptionMacro(" TImage must be either TFixedImage or TMovingImage.");

  m_ImageSamples.clear();

  ImageRegionType region = image->GetLargestPossibleRegion();
  IndexType firstCorner = region.GetIndex();
  IndexType corner;
  PointType point;

  SizeType size = region.GetSize();
  int cornerNumber = 1 << ImageDimension; // 2^ImageDimension

  for(int i=0; i<cornerNumber; i++)
    {
    int bit;
    for (int d=0; d<ImageDimension; d++)
      {
      bit = (int) (( i & (1 << d) ) != 0); // 0 or 1
      corner[d] = firstCorner[d] + bit * (size[d] - 1);
      }

    image->TransformIndexToPhysicalPoint(corner, point);
    m_ImageSamples.push_back(point);
    }
}

/**
 * Compute the maximum shift when one transform is changed to another
 */
template < class TFixedImage, class TMovingImage, class TTransform >
double
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::ComputeMaximumVoxelShift(ParametersType parameters,
                    ParametersType deltaParameters)
{
  double voxelShift = 0.0;

  ParametersType oldParameters = parameters;
  ParametersType newParameters(oldParameters.size());
  for (unsigned int p=0; p<oldParameters.size(); p++)
    {
    newParameters[p] = oldParameters[p] + deltaParameters[p];
    }
  TransformPointer oldTransform = TransformType::New();
  TransformPointer newTransform = TransformType::New();
  oldTransform->SetParameters(oldParameters);
  newTransform->SetParameters(newParameters);

  double distance;
  unsigned int numSamples = 0;
  PointType point, oldMappedPoint, newMappedPoint;

  ContinuousIndex<double, ImageDimension> oldMappedIndex, newMappedIndex;
  ContinuousIndex<double, ImageDimension> diffIndex;

  numSamples = m_ImageSamples.size();

  // find max shift by checking each sample point
  for (unsigned int c=0; c<numSamples; c++)
    {
    point = this->m_ImageSamples[c];
    oldMappedPoint = oldTransform->TransformPoint(point);
    newMappedPoint = newTransform->TransformPoint(point);

    if (this->m_TransformForward)
      {
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex(oldMappedPoint, oldMappedIndex);
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex(newMappedPoint, newMappedIndex);
      }
    else
      {
      this->m_FixedImage->TransformPhysicalPointToContinuousIndex(oldMappedPoint, oldMappedIndex);
      this->m_FixedImage->TransformPhysicalPointToContinuousIndex(newMappedPoint, newMappedIndex);
      }
    for (unsigned int d=0; d<ImageDimension; d++)
      {
      diffIndex[d] = oldMappedIndex[d] - newMappedIndex[d];
      }

    distance = ComputeLNorm(diffIndex);
    if ( voxelShift < distance )
      {
      voxelShift = distance;
      }
  } // end for numSamples
  return voxelShift;
}

/**
 * Compute the L-norm of a point
 */
template< class TFixedImage, class TMovingImage, class TTransform >
double
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::ComputeLNorm(Point<double, ImageDimension> point)
{
  double distance = 0;

  if (m_LNorm == 2) // Euclidean distance
    {
    for (unsigned int d=0; d<ImageDimension; d++)
      {
      distance += point[d] * point[d];
      }
    distance = vcl_sqrt(distance);
    }
  else if (m_LNorm == 1)
    {
    for (unsigned int d=0; d<ImageDimension; d++)
      {
      distance += vcl_abs(point[d]);
      }
    }
  else if (m_LNorm == -1) //L-infinity norm
    {
    for (unsigned int d=0; d<ImageDimension; d++)
      {
      if (distance < vcl_abs(point[d]))
        {
        distance = vcl_abs(point[d]);
        }
      }
    }
  else
    {
    std::cerr << "OptimizerParameterEstimator: norm undefined" << std::endl;
    }

  return distance;
}

template< class TFixedImage, class TMovingImage, class TTransform >
void
OptimizerParameterEstimator< TFixedImage, TMovingImage, TTransform >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "TransformType   = " << std::endl;
  os << indent << typeid(TransformType).name()  << std::endl;

  os << indent << "m_LNorm   = " << this->m_LNorm << std::endl;
  os << indent << "TransformForward   = " << this->m_TransformForward << std::endl;

  os << indent << "FixedImage   = " << std::endl;
  if( this->m_FixedImage )
    {
    os << indent << this->m_FixedImage  << std::endl;
    }
  else
    {
    os << indent << "None" << std::endl;
    }

  os << indent << "MovingImage   = " << std::endl;
  if( this->m_MovingImage )
    {
    os << indent << this->m_MovingImage  << std::endl;
    }
  else
    {
    os << indent << "None" << std::endl;
    }

  os << indent << "ImageSamples Size   = " << std::endl;
  os << indent << this->m_ImageSamples.size()  << std::endl;

}

}  // namespace itk

#endif /* __itkOptimizerParameterEstimator_txx */
