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
#ifndef __itkOptimizerParameterEstimator_hxx
#define __itkOptimizerParameterEstimator_hxx

#include "itkOptimizerParameterEstimator.h"

namespace itk
{

template< class TMetric >
OptimizerParameterEstimator< TMetric >
::OptimizerParameterEstimator()
{
  // Euclidean distance
  m_LNorm = 2;

  // estimate paramter scales of the moving transform
  m_TransformForward = true;

  // using maximum voxel shift to estimate scales
  m_ScaleStrategy = ScalesFromShift;

  // the metric object must be set before EstimateScales()
  m_Metric = NULL;

}

/** Compute parameter scales */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::EstimateScales(ScalesType &parameterScales)
{
  if ( m_Metric == (MetricPointer)NULL )
    {
    itkExceptionMacro(" OptimizerParameterEstimator: m_Metric == NULL.");
    }

  if (parameterScales.size() != this->GetTransform()->GetNumberOfParameters())
    {
    itkExceptionMacro(" The size of scales does not match that of transform parameters.");
    }

  this->SampleImageDomain();

  switch (m_ScaleStrategy)
    {
    case ScalesFromShift:
      {
      this->EstimateScalesFromMaximumShift(parameterScales);
      break;
      }
    case ScalesFromJacobian:
      {
      this->EstimateScalesFromJacobian(parameterScales);
      break;
      }
    default:
      {
      itkExceptionMacro(" Undefined strategy to decide scales.");
      }
    }
}

/** Get the transform being estimated */
template< class TMetric >
TransformBase *
OptimizerParameterEstimator< TMetric >
::GetTransform()
{
  if (m_TransformForward)
    {
    return m_MovingTransform.GetPointer();
    }
  else
    {
    return m_FixedTransform.GetPointer();
    }
}

/** Transform a physical point to its continuous index */
template< class TMetric >
template< class TContinuousIndexType >
void
OptimizerParameterEstimator< TMetric >
::TransformPointToContinuousIndex(const VirtualPointType &point,
                                  TContinuousIndexType &mappedIndex)
{
  if (this->m_TransformForward)
    {
    MovingPointType mappedPoint;
    mappedPoint = m_MovingTransform->TransformPoint(point);
    this->m_MovingImage->TransformPhysicalPointToContinuousIndex(mappedPoint, mappedIndex);
    }
  else
    {
    FixedPointType mappedPoint;
    mappedPoint = m_FixedTransform->TransformPoint(point);
    this->m_FixedImage->TransformPhysicalPointToContinuousIndex(mappedPoint, mappedIndex);
    }
}

/** Get the transform Jacobian w.r.t parameters at a point */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::ComputeJacobianWithRespectToParameters( const VirtualPointType  & point,
                                          JacobianType & jacobian ) const
{
  if (m_TransformForward)
    {
    m_MovingTransform->ComputeJacobianWithRespectToParameters(point, jacobian);
    }
  else
    {
    m_FixedTransform->ComputeJacobianWithRespectToParameters(point, jacobian);
    }
}

/** Estimate parameter scales from maximum voxel shift */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::EstimateScalesFromMaximumShift(ScalesType &parameterScales)
{
  double maxShift;
  unsigned int numPara = parameterScales.size();

  ParametersType deltaParameters(numPara);
  deltaParameters.Fill(0.0);

  // minNonZeroShift: the minimum non-zero shift.
  double minNonZeroShift = NumericTraits<double>::max();

  // compute voxel shift generated from each transform parameter
  for (unsigned int i=0; i<numPara; i++)
    {
    deltaParameters[i] = 1;
    maxShift = this->ComputeMaximumVoxelShift(deltaParameters);
    deltaParameters[i] = 0;

    if (maxShift != maxShift)
      {
      itkExceptionMacro("OptimizerParameterEstimator: maximum voxel shift is undefined with current parameters.");
      }

    parameterScales[i] = maxShift;
    if ( maxShift != 0 && maxShift < minNonZeroShift )
      {
      minNonZeroShift = maxShift;
      }
    }

  if (minNonZeroShift == NumericTraits<double>::max())
    {
    std::cout << std::endl
      << "Warning: any change in the parameters yields zero voxel shift."
      << std::endl
      << "Warning: nothing could be optimized."
      << std::endl;
    parameterScales.fill(0.0);
    }
  else
    {
    for (unsigned int i=0; i<numPara; i++)
      {
      if (parameterScales[i] == 0)
        {
        // To avoid division-by-zero in optimizers, assign a small value for a zero scale.
        parameterScales[i] = minNonZeroShift;
        }
      parameterScales[i] *= parameterScales[i];
      }
    }
}

/** Compute parameter scales from average jacobian norms.
 *  For each parameter, compute the squared norm of its transform Jacobian,
 *  then average the squared norm over the sample points. This average is
 *  used as the scale of this parameter.
 */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::EstimateScalesFromJacobian(ScalesType &parameterScales)
{
  unsigned int numPara = parameterScales.size();
  double normSquare;
  double *norms = new double[numPara];

  unsigned int numSamples = 0;
  unsigned int dim;

  numSamples = m_ImageSamples.size();
  if (m_TransformForward)
    {
    dim = MovingImageDimension;
    }
  else
    {
    dim = FixedImageDimension;
    }

  VirtualPointType point;

  for (unsigned int p=0; p<numPara; p++)
    {
    norms[p] = 0;
    }

  // checking each sample point
  for (unsigned int c=0; c<numSamples; c++)
    {
    point = m_ImageSamples[c];

    JacobianType jacobian;
    this->ComputeJacobianWithRespectToParameters( point, jacobian );

    for (unsigned int p=0; p<numPara; p++)
      {
      normSquare = 0;
      for (unsigned int d=0; d<dim; d++)
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

/** Initialize the sample points in the virtual image domain */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::SampleImageDomain()
{
  switch (m_ScaleStrategy)
    {
    case ScalesFromShift:
      {
      this->SampleWithCornerPoints();
      break;
      }
    case ScalesFromJacobian:
      {
      this->SampleImageDomainRandomly();
      break;
      }
    default:
      {
      itkExceptionMacro(" Undefined strategy to decide scales.");
      }
    }
}

/**
 * Sample the physical points of the virtual image in a uniform random distribution
 */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::SampleImageDomainRandomly()
{

  VirtualImagePointer image = this->m_VirtualImage;

  int numSamples = 1000;
  m_ImageSamples.clear();

  // Set up a random interator within the user specified fixed image region.
  typedef ImageRandomConstIteratorWithIndex<VirtualImageType> RandomIterator;
  RandomIterator randIter( image, image->GetLargestPossibleRegion() );

  randIter.SetNumberOfSamples( numSamples );
  randIter.GoToBegin();

  VirtualPointType point;
  VirtualIndexType index;

  for (int i=0; i<numSamples; i++)
    {
    index = randIter.GetIndex();
    image->TransformIndexToPhysicalPoint( index, point );
    m_ImageSamples.push_back(point);
    ++randIter;
    }
}

/**
 * Sample the physical points with image corners
 */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::SampleWithCornerPoints()
{
  VirtualImagePointer image = this->m_VirtualImage;
  m_ImageSamples.clear();

  VirtualRegionType region = image->GetLargestPossibleRegion();
  VirtualIndexType firstCorner = region.GetIndex();
  VirtualIndexType corner;
  VirtualPointType point;

  VirtualSizeType size = region.GetSize();
  int cornerNumber = 1 << VirtualImageDimension; // 2^ImageDimension

  for(int i=0; i<cornerNumber; i++)
    {
    int bit;
    for (int d=0; d<VirtualImageDimension; d++)
      {
      bit = (int) (( i & (1 << d) ) != 0); // 0 or 1
      corner[d] = firstCorner[d] + bit * (size[d] - 1);
      }

    image->TransformIndexToPhysicalPoint(corner, point);
    m_ImageSamples.push_back(point);
    }
}

/**
 * Compute the maximum shift when a transform is changed with deltaParameters
 */
template< class TMetric >
double
OptimizerParameterEstimator< TMetric >
::ComputeMaximumVoxelShift(ParametersType deltaParameters)
{
  //when ComputeMaximumVoxelShift is called without EstimateScales being called,
  //we need to make sure SampleImageDomain is called
  this->SampleImageDomain();

  double shift;
  if (m_TransformForward)
    {
    shift = this->ComputeTemplatedMaximumVoxelShift<MovingTransformType>(deltaParameters);
    }
  else
    {
    shift = this->ComputeTemplatedMaximumVoxelShift<FixedTransformType>(deltaParameters);
    }
  return shift;
}

/** The templated version of EstimateScalesFromMaximumShift.
 *  The template argument TTransform may be either MovingTransformType
 *  or FixedTransformType.
 */
template< class TMetric >
template< class TTransform >
double
OptimizerParameterEstimator< TMetric >
::ComputeTemplatedMaximumVoxelShift(ParametersType deltaParameters)
{
  double voxelShift = 0.0;
  unsigned int dim;
  if (m_TransformForward)
    {
    dim = MovingImageDimension;
    }
  else
    {
    dim = FixedImageDimension;
    }

  const ParametersType oldParameters = this->GetTransform()->GetParameters();
  ParametersType newParameters(oldParameters.size());
  for (unsigned int p=0; p<oldParameters.size(); p++)
    {
    newParameters[p] = oldParameters[p] + deltaParameters[p];
    }

  double distance;
  unsigned int numSamples = m_ImageSamples.size();

  VirtualPointType point;

  typedef ContinuousIndex<double, TTransform::OutputSpaceDimension> ContinuousIndexType;

  ContinuousIndexType newMappedIndex;
  ContinuousIndexType diffIndex;

  //store the old mapped indices to reduce calls to Transform::SetParameters()
  ContinuousIndexType *oldMappedIndices = new ContinuousIndexType[numSamples];

  // compute the indices mapped by the old transform
  for (unsigned int c=0; c<numSamples; c++)
    {
    point = this->m_ImageSamples[c];
    this->TransformPointToContinuousIndex<ContinuousIndexType>(point, oldMappedIndices[c]);
    } // end for numSamples

  // set the new parameters in the transform
  this->GetTransform()->SetParameters(newParameters);

  // compute the indices mapped by the new transform
  for (unsigned int c=0; c<numSamples; c++)
    {
    point = this->m_ImageSamples[c];
    this->TransformPointToContinuousIndex<ContinuousIndexType>(point, newMappedIndex);

    for (unsigned int d=0; d<dim; d++)
      {
      diffIndex[d] = oldMappedIndices[c][d] - newMappedIndex[d];
      }

    // find max shift by checking each sample point
    distance = ComputeLNorm<ContinuousIndexType>(diffIndex);
    if ( voxelShift < distance )
      {
      voxelShift = distance;
      }
  } // end for numSamples

  // restore the parameters in the transform
  this->GetTransform()->SetParameters(oldParameters);

  delete[] oldMappedIndices;

  return voxelShift;
}

/**
 * Compute the L-norm of a point
 */
template< class TMetric >
template< class TContinuousIndexType >
double
OptimizerParameterEstimator< TMetric >
::ComputeLNorm(TContinuousIndexType point)
{
  double distance = 0;
  unsigned int dim = TContinuousIndexType::IndexDimension;

  if (m_LNorm == 2) // Euclidean distance
    {
    for (unsigned int d=0; d<dim; d++)
      {
      distance += point[d] * point[d];
      }
    distance = vcl_sqrt(distance);
    }
  else if (m_LNorm == 1)
    {
    for (unsigned int d=0; d<dim; d++)
      {
      distance += vcl_abs(point[d]);
      }
    }
  else if (m_LNorm == -1) //L-infinity norm
    {
    for (unsigned int d=0; d<dim; d++)
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

/** Print the information about this class */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "MetricType   = " << std::endl;
  os << indent << typeid(MetricType).name()  << std::endl;

  os << indent << "FixedTransformType   = " << std::endl;
  os << indent << typeid(FixedTransformType).name()  << std::endl;

  os << indent << "MovingTransformType   = " << std::endl;
  os << indent << typeid(MovingTransformType).name()  << std::endl;

  os << indent << "m_FixedImage   = " << std::endl;
  if( this->m_FixedImage )
    {
    os << indent << this->m_FixedImage  << std::endl;
    }
  else
    {
    os << indent << "None" << std::endl;
    }

  os << indent << "m_MovingImage   = " << std::endl;
  if( this->m_MovingImage )
    {
    os << indent << this->m_MovingImage  << std::endl;
    }
  else
    {
    os << indent << "None" << std::endl;
    }

  os << indent << "m_VirtualImage   = " << std::endl;
  if( this->m_VirtualImage )
    {
    os << indent << this->m_VirtualImage  << std::endl;
    }
  else
    {
    os << indent << "None" << std::endl;
    }

  os << indent << "m_ImageSamples.size = " << std::endl;
  os << indent << this->m_ImageSamples.size()  << std::endl;

  os << indent << "m_LNorm = " << this->m_LNorm << std::endl;
  os << indent << "m_TransformForward = " << this->m_TransformForward << std::endl;
  os << indent << "m_ScaleStrategy = " << this->m_ScaleStrategy << std::endl;

}

}  // namespace itk

#endif /* __itkOptimizerParameterEstimator_txx */
