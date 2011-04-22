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
#ifndef __itkStdStreamLogOutput_h
#define __itkStdStreamLogOutput_h

#include <iostream>
#include <string>

#include "itkSimpleFastMutexLock.h"
#include "itkLogOutput.h"

namespace itk
{
/** \class StdStreamLogOutput
 *  \brief Class StdStreamLogOutput represents a standard stream output stream.
 *  This class provides thread safety for the standard stream output stream.
 *
 * \author Hee-Su Kim, Compute Science Dept. Kyungpook National University,
 *                     ISIS Center, Georgetown University.
 *
 *
 *  \ingroup OSSystemObjects LoggingObjects
 * \ingroup ITK-Common
 */

class ITKCommon_EXPORT StdStreamLogOutput:public LogOutput
{
public:

  typedef StdStreamLogOutput         Self;
  typedef LogOutput                  Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  typedef std::ostream StreamType;

  typedef std::ostream *StreamPointerType;

  itkTypeMacro(StdStreamLogOutput, LogOutput);

  itkNewMacro(StdStreamLogOutput);

  itkGetConstMacro(Stream, StreamPointerType);

  /** Set a standard stream pointer */
  void SetStream(StreamType & Stream);

  /** flush a buffer */
  virtual void Flush();

  /** Write to multiple outputs */
  virtual void Write(double timestamp);

  /** Write to a buffer */
  virtual void Write(std::string const & content);

  /** Write to a buffer */
  virtual void Write(std::string const & content, double timestamp);

protected:
  /** Constructor */
  StdStreamLogOutput();

  /** Destructor */
  virtual ~StdStreamLogOutput();

  void PrintSelf(std::ostream & os, Indent indent) const;

private:

  StreamPointerType m_Stream;

  SimpleFastMutexLock m_Mutex;
};
}

#endif //__itkStdStreamLogOutput_h
