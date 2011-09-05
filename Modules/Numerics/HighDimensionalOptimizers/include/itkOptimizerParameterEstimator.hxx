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
  // estimate paramter scales of the moving transform
  m_TransformForward = true;

  // using maximum voxel shift to estimate scales
  m_ScaleStrategy = ScalesFromShift;

  // the metric object must be set before EstimateScales()
  m_Metric = NULL;
  this->m_FixedImage   = NULL;
  this->m_MovingImage  = NULL;
  this->m_VirtualImage = NULL;

  this->m_MovingTransform = NULL;
  this->m_FixedTransform  = NULL;
}

/** SetMetric sets the metric used in the estimation process. SetMetric
 *  gets the images and transforms from the metric. Please make sure the metric
 *  has these members set when SetMetric(metric) is called.
 */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::SetMetric(MetricType *metric)
{
  if (metric == NULL)
    {
    itkExceptionMacro("OptimizerParameterEstimator: the metric argument is NULL");
    }

  this->m_Metric = metric;
  this->SetFixedImage(metric->GetFixedImage());
  this->SetMovingImage(metric->GetMovingImage());
  this->SetVirtualImage(metric->GetVirtualDomainImage());

  //the transforms will be modified with some delta parameters but be reset to its original values
  this->SetMovingTransform(const_cast<MovingTransformType *>(m_Metric->GetMovingTransform()));
  this->SetFixedTransform(const_cast<FixedTransformType *>(m_Metric->GetFixedTransform()));

  this->Modified();

  this->CheckInputs();
}

template< class TMetric >
bool
OptimizerParameterEstimator< TMetric >
::CheckInputs() const
{
if (m_Metric == (MetricPointer)NULL)
    {
    itkExceptionMacro("OptimizerParameterEstimator: the metric is NULL");
    }
  if (this->GetFixedImage() == NULL
    || this->GetMovingImage() == NULL
    || this->GetVirtualImage() == NULL)
    {
    itkExceptionMacro("OptimizerParameterEstimator: the image(s) in the metric is NULL");
    }
  if (this->GetMovingTransform() == NULL
    || this->GetFixedTransform() == NULL)
    {
    itkExceptionMacro("OptimizerParameterEstimator: the transform(s) in the metric is NULL.");
    }
  return true;
}

/** Compute parameter scales */
template< class TMetric >
void
OptimizerParameterEstimator< TMetric >
::EstimateScales(ScalesType &parameterScales)
{
  this->CheckInputs();

  parameterScales.SetSize(this->GetTransform()->GetNumberOfParameters());

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
      itkExceptionMacro("Undefined strategy to decide scales.");
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
  typedef ObjectToObjectMetric::InternalComputationValueType ValueType;
  ValueType maxShift;
  IndexValueType numPara = parameterScales.size();

  ParametersType deltaParameters(numPara);
  deltaParameters.Fill(0.0);

  // minNonZeroShift: the minimum non-zero shift.
  ValueType minNonZeroShift = NumericTraits<ValueType>::max();

  // compute voxel shift generated from each transform parameter
  for (IndexValueType i=0; i<numPara; i++)
    {
    deltaParameters[i] = 1;
    maxShift = this->ComputeMaximumVoxelShift(deltaParameters);
    deltaParameters[i] = 0;

    // check for NaN
    if (maxShift != maxShift)
      {
      itkExceptionMacro("OptimizerParameterEstimator: maximum voxel shift is undefined with current parameters.");
      }

    parameterScales[i] = maxShift;
    if ( maxShift > NumericTraits<ValueType>::epsilon() && maxShift < minNonZeroShift )
      {
      minNonZeroShift = maxShift;
      }
    }

  if (minNonZeroShift == NumericTraits<ValueType>::max())
    {
    std::cout << std::endl
      << "Warning: any change in the parameters yields zero voxel shift."
      << std::endl
      << "Warning: the default scales (1.0) are used to avoid division-by-zero."
      << std::endl;
    parameterScales.fill(1.0);
    }
  else
    {
    for (IndexValueType i=0; i<numPara; i++)
      {
      if (parameterScales[i] <= NumericTraits<ValueType>::epsilon())
        {
        // To avoid division-by-zero in optimizers, assign a small value for a zero scale.
        parameterScales[i] = minNonZeroShift * minNonZeroShift;
        }
      else
        {
        parameterScales[i] *= parameterScales[i];
        }
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
  IndexValueType numPara = parameterScales.size();

  typedef ObjectToObjectMetric::InternalComputationValueType ValueType;

  ValueType normSquare;
  itk::Array<ValueType> norms(numPara);

  IndexValueType numSamples = m_ImageSamples.size();
  IndexValueType dim;

  if (m_TransformForward)
    {
    dim = MovingImageDimension;
    }
  else
    {
    dim = FixedImageDimension;
    }

  for (IndexValueType p=0; p<numPara; p++)
    {
    norms[p] = 0;
    }

  // checking each sample point
  for (IndexValueType c=0; c<numSamples; c++)
    {
    VirtualPointType point = m_ImageSamples[c];

    JacobianType jacobian;
    this->ComputeJacobianWithRespectToParameters( point, jacobian );

    for (IndexValueType p=0; p<numPara; p++)
      {
      normSquare = 0;
      for (IndexValueType d=0; d<dim; d++)
        {
        normSquare = normSquare + jacobian[d][p] * jacobian[d][p];
        }
      norms[p] = norms[p] + normSquare;
      }
    } //for numSamples

  if (numSamples == 0)
    {
    parameterScales.Fill(1.0);
    }
  else
    {
    for (IndexValueType p=0; p<numPara; p++)
      {
      parameterScales[p] = norms[p] / numSamples;
      }
    }

}

/** Initialize the sample points in the virtual image domain.
 *  The results are stored into m_ImageSamples.
 */
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
      itkExceptionMacro("Undefined strategy to decide scales.");
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
  IndexValueType dim;
  if (m_TransformForward)
    {
    dim = MovingImageDimension;
    }
  else
    {
    dim = FixedImageDimension;
    }

  //We're purposely copying the parameters,
  //for temporary use of the transform to calculate the voxel shift.
  //After it is done, we will reset to the old parameters.
  const ParametersType oldParameters = this->GetTransform()->GetParameters();
  ParametersType newParameters(oldParameters.size());
  for (IndexValueType p=0; p<oldParameters.size(); p++)
    {
    newParameters[p] = oldParameters[p] + deltaParameters[p];
    }

  double distance;
  IndexValueType numSamples = m_ImageSamples.size();

  VirtualPointType point;

  typedef ContinuousIndex<double, TTransform::OutputSpaceDimension> ContinuousIndexType;

  ContinuousIndexType newMappedIndex;
  ContinuousIndexType diffIndex;

  //store the old mapped indices to reduce calls to Transform::SetParameters()
  std::vector<ContinuousIndexType> oldMappedIndices(numSamples);

  // compute the indices mapped by the old transform
  for (IndexValueType c=0; c<numSamples; c++)
    {
    point = this->m_ImageSamples[c];
    this->TransformPointToContinuousIndex<ContinuousIndexType>(point, oldMappedIndices[c]);
    } // end for numSamples

  // set the new parameters in the transform
  this->GetTransform()->SetParameters(newParameters);

  // compute the indices mapped by the new transform
  for (IndexValueType c=0; c<numSamples; c++)
    {
    point = this->m_ImageSamples[c];
    this->TransformPointToContinuousIndex<ContinuousIndexType>(point, newMappedIndex);

    // find max shift by checking each sample point
    distance = newMappedIndex.EuclideanDistanceTo(oldMappedIndices[c]);
    if ( voxelShift < distance )
      {
      voxelShift = distance;
      }
  } // end for numSamples

  // restore the parameters in the transform
  this->GetTransform()->SetParameters(oldParameters);

  return voxelShift;
}

/** Get the number of image dimensions */
template< class TMetric >
IndexValueType
OptimizerParameterEstimator< TMetric >
::GetImageDimension() const
{
  return Self::VirtualImageDimension;
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

  os << indent << "MovingTransformType  = " << std::endl;
  os << indent << typeid(MovingTransformType).name()  << std::endl;

  os << indent << "FixedImageType    = " << std::endl;
  os << indent << typeid(FixedImageType).name()  << std::endl;

  os << indent << "MovingImageType   = " << std::endl;
  os << indent << typeid(MovingImageType).name()  << std::endl;

  os << indent << "VirtualImageType  = " << std::endl;
  os << indent << typeid(VirtualImageType).name()  << std::endl;

  os << indent << "m_ImageSamples.size = " << std::endl;
  os << indent << this->m_ImageSamples.size()  << std::endl;

  os << indent << "m_TransformForward = " << this->m_TransformForward << std::endl;
  os << indent << "m_ScaleStrategy = " << this->m_ScaleStrategy << std::endl;

}

}  // namespace itk

#endif /* __itkOptimizerParameterEstimator_txx */
