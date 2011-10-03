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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#include "itkProcessObject.h"
#include "itkCommand.h"

#include <functional>
#include <algorithm>

namespace itk
{
/**
 * Instantiate object with no start, end, or progress methods.
 */
ProcessObject
::ProcessObject() :
  m_Inputs(),
  m_Outputs(),
  m_CachedInputReleaseDataFlags()
{
  m_NumberOfRequiredInputs = 0;

  m_NumberOfRequiredOutputs = 0;

  m_AbortGenerateData = false;
  m_Progress = 0.0f;
  m_Updating = false;

  m_Threader = MultiThreader::New();
  m_NumberOfThreads = m_Threader->GetNumberOfThreads();

  m_ReleaseDataBeforeUpdateFlag = true;

  m_NumberOfIndexedInputs = 0;
  m_NumberOfIndexedOutputs = 0;
}

/**
 * This is a default implementation to make sure we have something.
 * Once all the subclasses of ProcessObject provide an appopriate
 * MakeOutput(), then ProcessObject::MakeOutput() can be made pure
 * virtual.
 */
DataObject::Pointer
ProcessObject
::MakeOutput( const DataObjectIdentifierType & name )
{
  itkDebugMacro("MakeOutput(" << name << ")");
  if( this->IsIndexedName(name) )
    {
    return this->MakeOutput( this->MakeIndexFromName(name) );
    }
  return static_cast< DataObject * >( DataObject::New().GetPointer() );
}

DataObject::Pointer
ProcessObject
::MakeOutput(unsigned int)
{
  return static_cast< DataObject * >( DataObject::New().GetPointer() );
}


/**
 * Destructor for the ProcessObject class. We've got to
 * UnRegister() the use of any input classes.
 */
ProcessObject
::~ProcessObject()
{
  // Tell each output that we are going away.  If other objects have a
  // reference to one of these outputs, the data object will not be deleted
  // when the process object is deleted.  However, the data object's source
  // will still point back to the now nonexistent process object if we do not
  // clean things up now.
  for ( DataObjectPointerMap::iterator it = m_Outputs.begin(); it != m_Outputs.end(); it++ )
    {
    if ( it->second )
      {
      // let the output know we no longer want to associate with the object
      it->second->DisconnectSource(this, it->first);
      // let go of our reference to the data object
      it->second = 0;
      }
    }
}

/**
 * Called by constructor to set up input array.
 */
void
ProcessObject
::SetNumberOfIndexedInputs(unsigned int num)
{
  if( this->GetNumberOfIndexedInputs() != num )
    {
    // remove the extra inputs
    for( unsigned int i=num; i<this->GetNumberOfIndexedInputs(); i++ )
      {
      this->RemoveInput( this->MakeNameFromIndex( i ) );
      }
    m_NumberOfIndexedInputs = num;
    this->Modified();
    }
}

/**
 *
 */
void
ProcessObject
::SetNumberOfInputs(unsigned int num)
{
  this->SetNumberOfIndexedInputs(num);
}

/**
 * Get the number of specified inputs
 */
ProcessObject::DataObjectPointerArraySizeType
ProcessObject
::GetNumberOfValidRequiredInputs() const
{
  DataObjectPointerArraySizeType count = 0;
  for( unsigned int i = 0; i < m_NumberOfRequiredInputs; i++ )
    {
    if( this->GetInput( i ) )
      {
      count++;
      }
    }
  return count;
}

/**
 * Adds an input to the first null position in the input list.
 * Expands the list memory if necessary
 */
void
ProcessObject
::AddInput(DataObject *input)
{
  for ( unsigned idx = 0; idx < this->GetNumberOfIndexedInputs(); ++idx )
    {
    if ( !this->GetInput( idx ) )
      {
      this->SetNthInput( idx, input );
      return;
      }
    }
  this->SetNthInput( this->GetNumberOfIndexedInputs(), input );
}

/**
 * Remove an input.
 */
void
ProcessObject
::RemoveInput(const DataObjectIdentifierType & key)
{
  DataObjectPointerMap::iterator it = m_Inputs.find(key);
  if ( it != m_Inputs.end() )
    {
    m_Inputs.erase( it );
    this->Modified();
    }
}

/**
 *
 */
void
ProcessObject
::RemoveInput(DataObject * input)
{
  if( !input )
    {
    return;
    }
  unsigned int nb = this->GetNumberOfIndexedInputs();
  for( unsigned int i = 0; i < nb; i++ )
    {
    if( this->GetInput(i) == input )
      {
      this->RemoveInput( i );
      return;
      }
    }
  itkDebugMacro("tried to remove an input that was not in the list");
}

/**
 * Remove an input.
 *
 */
void
ProcessObject
::RemoveInput(unsigned int idx)
{
  if( idx == this->GetNumberOfIndexedInputs() - 1 )
    {
    // just remove the last indexed input
    this->SetNumberOfInputs( this->GetNumberOfIndexedInputs() -1 );
    }
  else
    {
    this->RemoveInput( this->MakeNameFromIndex( idx ) );
    }
}

/**
 *
 */
void
ProcessObject
::SetInput(const DataObjectIdentifierType & key, DataObject * input)
{
  DataObjectPointerMap::iterator it = m_Inputs.find(key);
  if( it == m_Inputs.end() )
    {
    // this is a whole new entry
    m_Inputs[key] = input;
    this->Modified();
    }
  else if( it->second.GetPointer() != input )
    {
    // there is already an entry, but not with the right object
    it->second = input;
    this->Modified();
    }
  // the entry is already there - there is nothing to do
}

/** Set an Input of this filter. This method
 * does Register()/UnRegister() manually to
 * deal with the fact that smart pointers aren't
 * around to do the reference counting.
 */
void
ProcessObject
::SetNthInput(unsigned int idx, DataObject *input)
{
  if ( idx >= this->GetNumberOfIndexedInputs() )
    {
    this->SetNumberOfInputs(idx + 1);
    }
  this->SetInput( this->MakeNameFromIndex(idx), input );
}

/**
 * Model a queue on the input list by providing a push back
 */
void
ProcessObject
::PushBackInput(const DataObject *input)
{
  this->SetNthInput( this->GetNumberOfIndexedInputs(), const_cast< DataObject * >( input ) );
}

/**
 * Model a stack on the input list by providing a pop back
 */
void
ProcessObject
::PopBackInput()
{
  if ( this->GetNumberOfIndexedInputs() > 0 )
    {
    this->SetNumberOfInputs( this->GetNumberOfIndexedInputs() - 1 );
    }
}

/**
 *
 */
void
ProcessObject
::PushFrontInput(const DataObject *input)
{
  int nb = this->GetNumberOfIndexedInputs();
  for( int i = nb-1; i >= 0; i-- )
    {
    this->SetNthInput( i+1, this->GetInput(i) );
    }
  this->SetNthInput( 0, const_cast< DataObject * >( input ) );
}

/**
 *
 */
void
ProcessObject
::PopFrontInput()
{
  unsigned int nb = this->GetNumberOfIndexedInputs();
  if( nb > 0 )
    {
    for( unsigned int i = 1; i < nb; i++ )
      {
      this->SetNthInput( i-1, this->GetInput( i ) );
      }
    this->SetNumberOfInputs( nb - 1 );
    }
}

/**
 * Remove an output.
 *
 */
void
ProcessObject
::RemoveOutput(const DataObjectIdentifierType & key)
{
  // find the input in the list of inputs
  DataObjectPointerMap::iterator it = m_Outputs.find(key);
  if ( it != m_Outputs.end() )
    {
    // let the output know we no longer want to associate with the object
    it->second->DisconnectSource( this, it->first );
    m_Outputs.erase( it );
    // let go of our reference to the data object
    this->Modified();
    }
}

/**
 *
 */
void
ProcessObject
::RemoveOutput(DataObject * output)
{
  if( !output )
    {
    return;
    }
  unsigned int nb = this->GetNumberOfIndexedOutputs();
  for( unsigned int i = 0; i < nb; i++ )
    {
    if( this->GetOutput(i) == output )
      {
      this->RemoveOutput( i );
      return;
      }
    }
  itkDebugMacro("tried to remove an output that was not in the list");
}

void
ProcessObject
::RemoveOutput(unsigned int idx)
{
  if( idx == this->GetNumberOfIndexedOutputs() - 1 )
    {
    // just remove the last indexed output
    this->SetNumberOfOutputs( this->GetNumberOfIndexedOutputs() -1 );
    }
  else
    {
    this->RemoveOutput( this->MakeNameFromIndex( idx ) );
    }
}

/**
  Set an Output of this filter. This method
 * does Register()/UnRegister() manually to
 * deal with the fact that smart pointers aren't
 * around to do the reference counting.
 */
void
ProcessObject
::SetOutput(const DataObjectIdentifierType & name, DataObject * output)
{
  // copy the key, because it might be destroyed in that method, so a reference
  // is not enough.
  DataObjectIdentifierType key = name;

  // does this change anything?
  DataObjectPointerMap::const_iterator it = m_Outputs.find(key);
  if ( it != m_Outputs.end() && it->second.GetPointer() == output)
    {
    return;
    }

  // Keep a handle to the original output and disconnect the old output from
  // the pipeline
  DataObjectPointer oldOutput;
  if ( m_Outputs[key] )
    {
    oldOutput = m_Outputs[key];
    m_Outputs[key]->DisconnectSource(this, key);
    }

  if ( output )
    {
    output->ConnectSource(this, key);
    }
  // save the current reference (which releases the previous reference)
  m_Outputs[key] = output;

  // if we are clearing an output, we need to create a new blank output
  // so we are prepared for the next Update(). this copies the requested
  // region ivar
  if ( !m_Outputs[key] )
    {
    itkDebugMacro(" creating new output object.");
    DataObjectPointer newOutput = this->MakeOutput(key);
    this->SetOutput(key, newOutput);

    // If we had an output object before, copy the requested region
    // ivars and release data flag to the the new output
    if ( oldOutput )
      {
      newOutput->SetRequestedRegion(oldOutput);
      newOutput->SetReleaseDataFlag( oldOutput->GetReleaseDataFlag() );
      }
    }

  this->Modified();
}

/**
 *
 */
void
ProcessObject
::SetNthOutput(unsigned int idx, DataObject *output)
{
  if ( idx >= this->GetNumberOfIndexedOutputs() )
    {
    this->SetNumberOfOutputs(idx + 1);
    }
  this->SetOutput( this->MakeNameFromIndex(idx), output );
}

/**
 * Adds an output to the first null position in the output list.
 * Expands the list memory if necessary
 */
void
ProcessObject
::AddOutput(DataObject *output)
{
  for ( unsigned int idx = 0; idx < this->GetNumberOfIndexedOutputs(); ++idx )
    {
    if ( !this->GetOutput( idx ) )
      {
      this->SetNthOutput( idx, output );
      return;
      }
    }
  this->SetNthOutput( this->GetNumberOfIndexedOutputs(), output );
}

/**
 * Called by constructor to set up output array.
 */
void
ProcessObject
::SetNumberOfIndexedOutputs(unsigned int num)
{
  if( this->GetNumberOfIndexedOutputs() != num )
    {
    // remove the extra outputs
    for( unsigned int i=num; i<this->GetNumberOfIndexedOutputs(); i++ )
      {
      this->RemoveOutput( this->MakeNameFromIndex( i ) );
      }
    m_NumberOfIndexedOutputs = num;
    this->Modified();
    }
}

/**
 *
 */
void
ProcessObject
::SetNumberOfOutputs(unsigned int num)
{
  this->SetNumberOfIndexedOutputs(num);
}

/**
 *
 */
DataObject *
ProcessObject
::GetOutput(const DataObjectIdentifierType & key)
{
  DataObjectPointerMap::iterator it = m_Outputs.find(key);
  if ( it == m_Outputs.end() )
    {
    return NULL;
    }
  return it->second.GetPointer();
}

const DataObject *
ProcessObject
::GetOutput(const DataObjectIdentifierType & key) const
{
  DataObjectPointerMap::const_iterator it = m_Outputs.find(key);
  if ( it == m_Outputs.end() )
    {
    return NULL;
    }
  return it->second.GetPointer();
}

/**
 *
 */
DataObject *
ProcessObject
::GetOutput(unsigned int i)
{
  return this->GetOutput( this->MakeNameFromIndex(i) );
}

const DataObject *
ProcessObject
::GetOutput(unsigned int i) const
{
  return this->GetOutput( this->MakeNameFromIndex(i) );
}

const DataObject *
ProcessObject
::GetPrimaryOutput() const
{
  return this->GetOutput("Primary");
}

DataObject *
ProcessObject
::GetPrimaryOutput()
{
  return this->GetOutput("Primary");
}

void
ProcessObject
::SetPrimaryOutput(DataObject * object)
{
  this->SetOutput("Primary", object);
}

bool
ProcessObject
::HasOutput( const DataObjectIdentifierType & key ) const
{
  DataObjectPointerMap::const_iterator it = m_Outputs.find(key);
  return it != m_Outputs.end();
}

ProcessObject::NameArray
ProcessObject
::GetOutputNames() const
{
  NameArray res(m_Outputs.size());
  for ( DataObjectPointerMap::const_iterator it = m_Outputs.begin(); it != m_Outputs.end(); it++ )
    {
    res.push_back( it->first );
    }

  return res;
}

// ProcessObject::ConstDataObjectPointerArray
// ProcessObject
// ::GetOutputs() const
// {
//   ConstDataObjectPointerArray res;
//   res.reserve(m_Outputs.size());
//   for ( DataObjectPointerMap::const_iterator it = m_Outputs.begin(); it != m_Outputs.end(); it++ )
//     {
//     res.push_back( it->second.GetPointer() );
//     }
//   return res;
// }

ProcessObject::DataObjectPointerArray
ProcessObject
::GetOutputs()
{
  DataObjectPointerArray res;
  res.reserve(m_Outputs.size());
  for ( DataObjectPointerMap::iterator it = m_Outputs.begin(); it != m_Outputs.end(); it++ )
    {
    res.push_back( it->second.GetPointer() );
    }
  return res;
}

ProcessObject::DataObjectPointerArraySizeType
ProcessObject
::GetNumberOfIndexedOutputs() const
{
  if( m_NumberOfIndexedOutputs > 0 )
    {
    return m_NumberOfIndexedOutputs;
    }
  if( this->GetPrimaryOutput() )
    {
    return 1;
    }
  return 0;
}

// ProcessObject::ConstDataObjectPointerArray
// ProcessObject
// ::GetIndexedOutputs() const
// {
//   ConstDataObjectPointerArray res(this->GetNumberOfIndexedOutputs());
//   for ( unsigned int i = 0; i < this->GetNumberOfIndexedOutputs(); i++ )
//     {
//     res[i] = this->GetOutput(i);
//     }
//   return res;
// }

ProcessObject::DataObjectPointerArray
ProcessObject
::GetIndexedOutputs()
{
  DataObjectPointerArray res(this->GetNumberOfIndexedOutputs());
  for ( unsigned int i = 0; i < this->GetNumberOfIndexedOutputs(); i++ )
    {
    res[i] = this->GetOutput(i);
    }
  return res;
}

/**
 *
 */
DataObject *
ProcessObject
::GetInput(const DataObjectIdentifierType & key)
{
  DataObjectPointerMap::iterator it = m_Inputs.find(key);
  if ( it == m_Inputs.end() )
    {
    return NULL;
    }
  return it->second.GetPointer();
}

const DataObject *
ProcessObject
::GetInput(const DataObjectIdentifierType & key) const
{
  DataObjectPointerMap::const_iterator it = m_Inputs.find(key);
  if ( it == m_Inputs.end() )
    {
    return NULL;
    }
  return it->second.GetPointer();
}

/**
 *
 */
DataObject *
ProcessObject
::GetInput(unsigned int i)
{
  return this->GetInput( this->MakeNameFromIndex(i) );
}

const DataObject *
ProcessObject
::GetInput(unsigned int i) const
{
  return this->GetInput( this->MakeNameFromIndex(i) );
}

const DataObject *
ProcessObject
::GetPrimaryInput() const
{
  return this->GetInput("Primary");
}

DataObject *
ProcessObject
::GetPrimaryInput()
{
  return this->GetInput("Primary");
}

void
ProcessObject
::SetPrimaryInput(DataObject * object)
{
  this->SetInput("Primary", object);
}

bool
ProcessObject
::HasInput( const DataObjectIdentifierType & key ) const
{
  DataObjectPointerMap::const_iterator it = m_Inputs.find(key);
  return it != m_Inputs.end();
}

ProcessObject::NameArray
ProcessObject
::GetInputNames() const
{
  NameArray res(m_Inputs.size());
  for ( DataObjectPointerMap::const_iterator it = m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    res.push_back( it->first );
    }
  return res;
}

// ProcessObject::ConstDataObjectPointerArray
// ProcessObject
// ::GetInputs() const
// {
//   ConstDataObjectPointerArray res;
//   res.reserve(m_Inputs.size());
//   for ( DataObjectPointerMap::const_iterator it = m_Inputs.begin(); it != m_Inputs.end(); it++ )
//     {
//     res.push_back( it->second.GetPointer() );
//     }
//   return res;
// }

ProcessObject::DataObjectPointerArray
ProcessObject
::GetInputs()
{
  DataObjectPointerArray res;
  res.reserve(m_Inputs.size());
  for ( DataObjectPointerMap::iterator it = m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    res.push_back( it->second.GetPointer() );
    }
  return res;
}

ProcessObject::DataObjectPointerArraySizeType
ProcessObject
::GetNumberOfIndexedInputs() const
{
  if( m_NumberOfIndexedInputs > 0 )
    {
    return m_NumberOfIndexedInputs;
    }
  if( this->GetPrimaryInput() )
    {
    return 1;
    }
  return 0;
}

// ProcessObject::ConstDataObjectPointerArray
// ProcessObject
// ::GetIndexedInputs() const
// {
//   ConstDataObjectPointerArray res(this->GetNumberOfIndexedInputs());
//   for ( unsigned int i = 0; i < this->GetNumberOfIndexedInputs(); i++ )
//     {
//     res[i] = this->GetInput(i);
//     }
//   return res;
// }

ProcessObject::DataObjectPointerArray
ProcessObject
::GetIndexedInputs()
{
  DataObjectPointerArray res(this->GetNumberOfIndexedInputs());
  for ( unsigned int i = 0; i < this->GetNumberOfIndexedInputs(); i++ )
    {
    res[i] = this->GetInput(i);
    }
  return res;
}

ProcessObject::DataObjectIdentifierType
ProcessObject
::MakeNameFromIndex(unsigned int idx) const
{
  if( idx == 0 )
    {
    return "Primary";
    }
  DataObjectIdentifierType baseName = "IndexedDataObject";
  std::ostringstream oss;
  oss << idx;
  return baseName + oss.str();
}

unsigned int
ProcessObject
::MakeIndexFromName(const DataObjectIdentifierType & name) const
{
  if( name == "Primary" )
    {
    itkDebugMacro("MakeIndexFromName("<<name<<") -> 0");
    return 0;
    }
  DataObjectIdentifierType baseName = "IndexedDataObject";
  unsigned int baseSize = baseName.size();
  if( name.size() <= baseSize || name.substr(0, baseSize) != baseName )
    {
    itkDebugMacro("MakeIndexFromName("<<name<<") -> exception bad base name");
    itkExceptionMacro(<<"Not an indexed data object: " << name);
    }
  DataObjectIdentifierType idxStr = name.substr(baseSize);
  unsigned int idx;
  if( !(std::istringstream(idxStr) >> idx) )
    {
    itkDebugMacro("MakeIndexFromName("<<name<<") -> exception not an index");
    itkExceptionMacro(<<"Not an indexed data object: " << name);
    }
  itkDebugMacro("MakeIndexFromName("<<name<<") -> "<< idx);
  return idx;
}

bool
ProcessObject
::IsIndexedName(const DataObjectIdentifierType & name) const
{
  if( name == "Primary" )
    {
    return true;
    }
  DataObjectIdentifierType baseName = "IndexedDataObject";
  unsigned int baseSize = baseName.size();
  if( name.size() <= baseSize || name.substr(0, baseSize) != baseName )
    {
    return false;
    }
  DataObjectIdentifierType idxStr = name.substr(baseSize);
  unsigned int idx;
  return (std::istringstream(idxStr) >> idx);
}

/**
 * Update the progress of the process object. If a ProgressMethod exists,
 * execute it. Then set the Progress ivar to amount. The parameter amount
 * should range between (0,1).
 */
void
ProcessObject
::UpdateProgress(float amount)
{
  m_Progress = amount;
  this->InvokeEvent( ProgressEvent() );
}

/**
 *
 */
bool
ProcessObject
::GetReleaseDataFlag() const
{
  if ( this->GetPrimaryOutput() )
    {
    return this->GetPrimaryOutput()->GetReleaseDataFlag();
    }
  itkWarningMacro(<< "Output doesn't exist!");
  return false;
}

/**
 *
 */
void
ProcessObject
::SetReleaseDataFlag(bool val)
{
  for ( DataObjectPointerMap::iterator it=m_Outputs.begin(); it != m_Outputs.end(); it++ )
    {
    if ( it->second )
      {
      it->second->SetReleaseDataFlag(val);
      }
    }
}

/**
 *
 */
void
ProcessObject
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  Indent indent2 = indent.GetNextIndent();
  if ( !m_Inputs.empty() )
    {
    os << indent << "Inputs: " << std::endl;
    for ( DataObjectPointerMap::const_iterator it = m_Inputs.begin(); it != m_Inputs.end(); it++ )
      {
      os << indent2 << it->first << ": (" << it->second.GetPointer() << ")" << std::endl;
      }
    }
  else
    {
    os << indent << "No Inputs\n";
    }

  if ( !m_Outputs.empty() )
    {
    os << indent << "Outputs: " << std::endl;
    for ( DataObjectPointerMap::const_iterator it = m_Outputs.begin(); it != m_Outputs.end(); it++ )
      {
      os << indent2 << it->first << ": (" << it->second.GetPointer() << ")" << std::endl;
      }
    }
  else
    {
    os << indent << "No Outputs\n";
    }

  os << indent << "NumberOfIndexedInputs: "
     << m_NumberOfIndexedInputs << std::endl;
  os << indent << "NumberOfRequiredInputs: "
     << m_NumberOfRequiredInputs << std::endl;
  os << indent << "NumberOfIndexedOutputs: "
     << m_NumberOfIndexedOutputs << std::endl;
  os << indent << "NumberOfRequiredOutputs: "
     << m_NumberOfRequiredOutputs << std::endl;

  os << indent << "Number Of Threads: "
     << m_NumberOfThreads << std::endl;

  os << indent << "ReleaseDataFlag: "
     << ( this->GetReleaseDataFlag() ? "On" : "Off" ) << std::endl;

  os << indent << "ReleaseDataBeforeUpdateFlag: "
     << ( m_ReleaseDataBeforeUpdateFlag ? "On" : "Off" ) << std::endl;

  os << indent << "AbortGenerateData: " << ( m_AbortGenerateData ? "On" : "Off" ) << std::endl;
  os << indent << "Progress: " << m_Progress << std::endl;

  os << indent << "Multithreader: " << std::endl;
  m_Threader->PrintSelf( os, indent.GetNextIndent() );
}

/**
 * The following methods are used to coordinate the execution of the
 * data processing pipeline.
 */

/**
 *
 */
void
ProcessObject
::Update()
{
  if ( this->GetPrimaryOutput() )
    {
    this->GetPrimaryOutput()->Update();
    }
}

void
ProcessObject
::ResetPipeline()
{
  if ( this->GetPrimaryOutput() )
    {
    this->GetPrimaryOutput()->ResetPipeline();
    }
  else
    {
    // important to make it work on process objects without outputs
    this->PropagateResetPipeline();
    }
}

void
ProcessObject
::PropagateResetPipeline()
{
  //
  // Reset this object.
  //
  // Clear the updating flag.
  m_Updating = 0;

  //
  // Loop through the inputs
  //
  for ( DataObjectPointerMap::iterator it=m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    if ( it->second )
      {
      it->second->PropagateResetPipeline();
      }
    }
}

/**
 *
 */
void
ProcessObject
::VerifyPreconditions()
{

  /**
    * Count the number of required inputs which have been assigned
    */
  DataObjectPointerArraySizeType ninputs = this->GetNumberOfValidRequiredInputs();

  if ( ninputs < m_NumberOfRequiredInputs )
    {
    itkExceptionMacro(<< "At least " << m_NumberOfRequiredInputs
                      << " inputs are required but only " << ninputs
                      << " are specified.");
    }
}


/**
 *
 */
void
ProcessObject
::VerifyInputInformation()
{
}


/**
 *
 */
void
ProcessObject
::UpdateOutputInformation()
{
  unsigned long                  t1, t2;
  DataObject *                   input;
  DataObject *                   output;

  /**
   * Watch out for loops in the pipeline
   */
  if ( m_Updating )
    {
    /**
     * Since we are in a loop, we will want to update. But if
     * we don't modify this filter, then we will not execute
     * because our OutputInformationMTime will be more recent than
     * the MTime of our output.
     */
    this->Modified();
    return;
    }

  /**
   * Verify that the process object has been configured correctly,
   * that all required inputs are set, and needed parameters are set
   * appropriately before we continue the pipeline, i.e. is the filter
   * in a state that it can be run.
   */
  this->VerifyPreconditions();

  /**
   * We now wish to set the PipelineMTime of each output DataObject to
   * the largest of this ProcessObject's MTime, all input DataObject's
   * PipelineMTime, and all input's MTime.  We begin with the MTime of
   * this ProcessObject.
   */
  t1 = this->GetMTime();

  /**
   * Loop through the inputs
   */
  for ( DataObjectPointerMap::iterator it=m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    if ( it->second )
      {
      input = it->second;

      /**
       * Propagate the UpdateOutputInformation call
       */
      m_Updating = true;
      input->UpdateOutputInformation();
      m_Updating = false;

      /**
       * What is the PipelineMTime of this input? Compare this against
       * our current computation to find the largest one.
       */
      t2 = input->GetPipelineMTime();

      if ( t2 > t1 )
        {
        t1 = t2;
        }

      /**
       * Pipeline MTime of the input does not include the MTime of the
       * data object itself. Factor these mtimes into the next PipelineMTime
       */
      t2 = input->GetMTime();
      if ( t2 > t1 )
        {
        t1 = t2;
        }
      }
    }

  /**
   * Call GenerateOutputInformation for subclass specific information.
   * Since UpdateOutputInformation propagates all the way up the pipeline,
   * we need to be careful here to call GenerateOutputInformation only if
   * necessary. Otherwise, we may cause this source to be modified which
   * will cause it to execute again on the next update.
   */
  if ( t1 > m_OutputInformationMTime.GetMTime() )
    {
    for ( DataObjectPointerMap::iterator it=m_Outputs.begin(); it != m_Outputs.end(); it++ )
      {
      output = it->second;
      if ( output )
        {
        output->SetPipelineMTime(t1);
        }
      }

    /**
     * Verify that all the inputs are consistent with the
     * requirements of the filter. For example, subclasses might want
     * to ensure all the inputs are in the same coordinate frame.
     */
    this->VerifyInputInformation();

    /**
     * Finally, generate the output information.
     */
    this->GenerateOutputInformation();

    /**
     * Keep track of the last time GenerateOutputInformation() was called
     */
    m_OutputInformationMTime.Modified();
    }
}

/**
 *
 */
void
ProcessObject
::PropagateRequestedRegion(DataObject *output)
{
  /**
   * check flag to avoid executing forever if there is a loop
   */
  if ( m_Updating )
    {
    return;
    }

  /**
   * Give the subclass a chance to indicate that it will provide
   * more data then required for the output. This can happen, for
   * example, when a source can only produce the whole output.
   * Although this is being called for a specific output, the source
   * may need to enlarge all outputs.
   */
  this->EnlargeOutputRequestedRegion(output);

  /**
   * Give the subclass a chance to define how to set the requested
   * regions for each of its outputs, given this output's requested
   * region.  The default implementation is to make all the output
   * requested regions the same.  A subclass may need to override this
   * method if each output is a different resolution.
   */
  this->GenerateOutputRequestedRegion(output);

  /**
   * Give the subclass a chance to request a larger requested region on
   * the inputs. This is necessary when, for example, a filter
   * requires more data at the "internal" boundaries to
   * produce the boundary values - such as an image filter that
   * derives a new pixel value by applying some operation to a
   * neighborhood of surrounding original values.
   */
  this->GenerateInputRequestedRegion();

  /**
   * Now that we know the input requested region, propagate this
   * through all the inputs.
   */
  m_Updating = true;
  for ( DataObjectPointerMap::iterator it=m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    if ( it->second )
      {
      it->second->PropagateRequestedRegion();
      }
    }
  m_Updating = false;
}

/**
 * By default we require all the input to produce the output. This is
 * overridden in the subclasses since we can often produce the output with
 * just a portion of the input data.
 */
void
ProcessObject
::GenerateInputRequestedRegion()
{
  for ( DataObjectPointerMap::iterator it=m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    if ( it->second )
      {
      it->second->SetRequestedRegionToLargestPossibleRegion();
      }
    }
}

/**
 * By default we set all the output requested regions to be the same.
 */
void
ProcessObject
::GenerateOutputRequestedRegion(DataObject *output)
{
  for ( DataObjectPointerMap::iterator it=m_Outputs.begin(); it != m_Outputs.end(); it++ )
    {
    if ( it->second && it->second != output )
      {
      it->second->SetRequestedRegion(output);
      }
    }
}

/**
 *
 */
void
ProcessObject
::PrepareOutputs()
{
  if ( this->GetReleaseDataBeforeUpdateFlag() )
    {
    for ( DataObjectPointerMap::iterator it=m_Outputs.begin(); it != m_Outputs.end(); it++ )
      {
      if ( it->second )
        {
        it->second->PrepareForNewData();
        }
      }
    }
}

/**
 *
 */
void
ProcessObject
::ReleaseInputs()
{
  for ( DataObjectPointerMap::iterator it=m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    if ( it->second )
      {
      if ( it->second->ShouldIReleaseData() )
        {
        it->second->ReleaseData();
        }
      }
    }
}

/**
 *
 */
void
ProcessObject
::UpdateOutputData( DataObject * itkNotUsed(output) )
{
  /**
   * prevent chasing our tail
   */
  if ( m_Updating )
    {
    return;
    }

  /**
   * Prepare all the outputs. This may deallocate previous bulk data.
   */
  this->PrepareOutputs();

  /**
   * Propagate the update call - make sure everything we
   * might rely on is up-to-date
   * Must call PropagateRequestedRegion before UpdateOutputData if multiple
   * inputs since they may lead back to the same data object.
   */
  m_Updating = true;
  if ( m_Inputs.size() == 1 )
    {
    if ( this->GetPrimaryInput() )
      {
      this->GetPrimaryInput()->UpdateOutputData();
      }
    }
  else
    {
    for ( DataObjectPointerMap::iterator it=m_Inputs.begin(); it != m_Inputs.end(); it++ )
      {
      if ( it->second )
        {
        it->second->PropagateRequestedRegion();
        it->second->UpdateOutputData();
        }
      }
    }

  /**
   * Cache the state of any ReleaseDataFlag's on the inputs. While the
   * filter is executing, we need to set the ReleaseDataFlag's on the
   * inputs to false in case the current filter is implemented using a
   * mini-pipeline (which will try to release the inputs).  After the
   * filter finishes, we restore the state of the ReleaseDataFlag's
   * before the call to ReleaseInputs().
   */
  this->CacheInputReleaseDataFlags();

  /**
   * Tell all Observers that the filter is starting
   */
  this->InvokeEvent( StartEvent() );

  /**
   * GenerateData this object - we have not aborted yet, and our progress
   * before we start to execute is 0.0.
   */
  m_AbortGenerateData = false;
  m_Progress = 0.0f;

  try
    {
    this->GenerateData();
    }
  catch ( ProcessAborted & excp )
    {
    this->InvokeEvent( AbortEvent() );
    this->ResetPipeline();
    this->RestoreInputReleaseDataFlags();
    throw excp;
    }
  catch (...)
    {
    this->ResetPipeline();
    this->RestoreInputReleaseDataFlags();
    throw;
    }

  /**
   * If we ended due to aborting, push the progress up to 1.0 (since
   * it probably didn't end there)
   *
   */
  if ( m_AbortGenerateData )
    {
    this->UpdateProgress(1.0f);
    }

  /**
   * Notify end event observers
   */
  this->InvokeEvent( EndEvent() );

  /**
   * Now we have to mark the data as up to date.
   */
  for ( DataObjectPointerMap::iterator it=m_Outputs.begin(); it != m_Outputs.end(); it++ )
    {
    if ( it->second )
      {
      it->second->DataHasBeenGenerated();
      }
    }

  /**
   * Restore the state of any input ReleaseDataFlags
   */
  this->RestoreInputReleaseDataFlags();

  /**
   * Release any inputs if marked for release
   */
  this->ReleaseInputs();

  // Mark that we are no longer updating the data in this filter
  m_Updating = false;
}

/**
 *
 */
void
ProcessObject
::CacheInputReleaseDataFlags()
{
  m_CachedInputReleaseDataFlags.clear();
  for ( DataObjectPointerMap::iterator it = m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    if ( it->second )
      {
      m_CachedInputReleaseDataFlags[it->first] = it->second->GetReleaseDataFlag();
      it->second->ReleaseDataFlagOff();
      }
    else
      {
      m_CachedInputReleaseDataFlags[it->first] = false;
      }
    }
}

/**
 *
 */
void
ProcessObject
::RestoreInputReleaseDataFlags()
{
  for ( DataObjectPointerMap::iterator it = m_Inputs.begin(); it != m_Inputs.end(); it++ )
    {
    if ( it->second )
      {
      it->second->SetReleaseDataFlag(m_CachedInputReleaseDataFlags[it->first]);
      }
    }
  m_CachedInputReleaseDataFlags.clear();
}

/**
 * Default implementation - copy information from first input to all outputs
 */
void
ProcessObject
::GenerateOutputInformation()
{
  DataObject * input = this->GetPrimaryInput();

  if ( input )
    {
    for ( DataObjectPointerMap::iterator it=m_Outputs.begin(); it != m_Outputs.end(); it++ )
      {
      if ( it->second )
        {
        it->second->CopyInformation(input);
        }
      }
    }
}

/**
 *
 */
void
ProcessObject
::UpdateLargestPossibleRegion()
{
  this->UpdateOutputInformation();

  if ( this->GetPrimaryOutput() )
    {
    this->GetPrimaryOutput()->SetRequestedRegionToLargestPossibleRegion();
    this->GetPrimaryOutput()->Update();
    }
}
} // end namespace itk
