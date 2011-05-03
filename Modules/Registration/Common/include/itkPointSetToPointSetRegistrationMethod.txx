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
#ifndef __itkPointSetToPointSetRegistrationMethod_txx
#define __itkPointSetToPointSetRegistrationMethod_txx

#include "itkPointSetToPointSetRegistrationMethod.h"

namespace itk
{
/**
 * Constructor
 */
template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::PointSetToPointSetRegistrationMethod()
{
  this->SetNumberOfRequiredOutputs(1);    // for the Transform

  m_FixedPointSet   = 0; // has to be provided by the user.
  m_MovingPointSet  = 0; // has to be provided by the user.
  m_Transform       = 0; // has to be provided by the user.
  m_Metric          = 0; // has to be provided by the user.
  m_Optimizer       = 0; // has to be provided by the user.

  m_InitialTransformParameters = ParametersType(1);
  m_LastTransformParameters = ParametersType(1);

  m_InitialTransformParameters.Fill(0.0f);
  m_LastTransformParameters.Fill(0.0f);

  TransformOutputPointer transformDecorator =
    static_cast< TransformOutputType * >(
      this->MakeOutput(0).GetPointer() );

  this->ProcessObject::SetNthOutput( 0, transformDecorator.GetPointer() );
}

/*
 * Set the initial transform parameters
 */
template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
void
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::SetInitialTransformParameters(const ParametersType & param)
{
  m_InitialTransformParameters = param;
  this->Modified();
}

/**
 * Initialize by setting the interconnects between components.
 */
template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
void
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::Initialize()
throw ( ExceptionObject )
{
  if ( !m_FixedPointSet )
    {
    itkExceptionMacro(<< "FixedPointSet is not present");
    }

  if ( !m_MovingPointSet )
    {
    itkExceptionMacro(<< "MovingPointSet is not present");
    }

  if ( !m_Metric )
    {
    itkExceptionMacro(<< "Metric is not present");
    }

  if ( !m_Optimizer )
    {
    itkExceptionMacro(<< "Optimizer is not present");
    }

  if ( !m_Transform )
    {
    itkExceptionMacro(<< "Transform is not present");
    }

  // Setup the metric
  m_Metric->SetMovingPointSet(m_MovingPointSet);
  m_Metric->SetFixedPointSet(m_FixedPointSet);
  m_Metric->SetTransform(m_Transform);

  m_Metric->Initialize();

  // Setup the optimizer
  m_Optimizer->SetCostFunction(m_Metric);

  // Validate initial transform parameters
  if ( m_InitialTransformParameters.Size() !=
       m_Transform->GetNumberOfParameters() )
    {
    itkExceptionMacro(<< "Size mismatch between initial parameter and transform");
    }

  m_Optimizer->SetInitialPosition(m_InitialTransformParameters);

  //
  // Connect the transform to the Decorator.
  //
  TransformOutputType *transformOutput =
    static_cast< TransformOutputType * >( this->ProcessObject::GetOutput(0) );

  transformOutput->Set( m_Transform.GetPointer() );
}

/*
 * Starts the Registration Process
 */
template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
void
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::StartRegistration(void)
{
  // StartRegistration is an old API from before
  // the RegistrationMethod was a subclass of ProcessObject.
  // Historically, one could call StartRegistration() instead of
  // calling Update().  However, when called directly by the user, the
  // inputs to the RegistrationMethod may not be up to date.  This
  // may cause an unexpected behavior.
  //
  // Since we cannot eliminate StartRegistration for backward
  // compability reasons, we check whether StartRegistration was
  // called directly or whether Update() (which in turn called
  // StartRegistration()).
  if ( !m_Updating )
    {
    this->Update();
    }
  else
    {
    try
      {
      // initialize the interconnects between components
      this->Initialize();
      }
    catch ( ExceptionObject & err )
      {
      m_LastTransformParameters = ParametersType(1);
      m_LastTransformParameters.Fill(0.0f);

      // pass exception to caller
      throw err;
      }

    try
      {
      // do the optimization
      m_Optimizer->StartOptimization();
      }
    catch ( ExceptionObject & err )
      {
      // An error has occurred in the optimization.
      // Update the parameters
      m_LastTransformParameters = m_Optimizer->GetCurrentPosition();

      // Pass exception to caller
      throw err;
      }

    // get the results
    m_LastTransformParameters = m_Optimizer->GetCurrentPosition();

    m_Transform->SetParameters(m_LastTransformParameters);
    }
}

/**
 * PrintSelf
 */
template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
void
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Metric: " << m_Metric.GetPointer() << std::endl;
  os << indent << "Optimizer: " << m_Optimizer.GetPointer() << std::endl;
  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "Fixed PointSet: " << m_FixedPointSet.GetPointer() << std::endl;
  os << indent << "Moving PointSet: " << m_MovingPointSet.GetPointer() << std::endl;
  os << indent << "Initial Transform Parameters: " << m_InitialTransformParameters << std::endl;
  os << indent << "Last    Transform Parameters: " << m_LastTransformParameters << std::endl;
}

template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
void
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::GenerateData()
{
  this->StartRegistration();
}

/**
 *  Get Output
 */
template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
const typename PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                                     TMovingPointSet,
                                                     TValueType >
::TransformOutputType *
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::GetOutput() const
{
  return static_cast< const TransformOutputType * >( this->ProcessObject::GetOutput(0) );
}

template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
DataObject::Pointer
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::MakeOutput(unsigned int output)
{
  switch ( output )
    {
    case 0:
      return static_cast< DataObject * >( TransformOutputType::New().GetPointer() );
      break;
    default:
      itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
      return 0;
    }
}

/**
 *
 */
template< typename TFixedPointSet,
          typename TMovingPointSet,
          typename TValueType >
unsigned long
PointSetToPointSetRegistrationMethod< TFixedPointSet,
                                      TMovingPointSet,
                                      TValueType >
::GetMTime() const
{
  unsigned long mtime = Superclass::GetMTime();
  unsigned long m;

  // Some of the following should be removed once ivars are put in the
  // input and output lists

  if ( m_Transform )
    {
    m = m_Transform->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_Metric )
    {
    m = m_Metric->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_Optimizer )
    {
    m = m_Optimizer->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_FixedPointSet )
    {
    m = m_FixedPointSet->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_MovingPointSet )
    {
    m = m_MovingPointSet->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  return mtime;
}
} // end namespace itk

#endif
