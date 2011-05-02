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
#ifndef __itkTransformIOBase_h
#define __itkTransformIOBase_h
#include "itkLightProcessObject.h"
#include "itkTransformBase.h"
#include <list>
#include <iostream>
#include <fstream>
#include <string>

namespace itk
{
/** \class TransformIOBase
 * \brief Abstract superclass defining the Transform IO interface.
 *
 * TransformIOBase is a pure virtual base class for dervied classes that
 * read/write transform data. It's used by the TransformFileReader and
 * TransformFileWriter classes.  End users don't directly manipulate classes
 * derived from TransformIOBase; the TransformIOFactory is used by
 * the Reader/Writer to pick a concrete derived class to do the actual
 * reading/writing of transforms.
 * \ingroup ITK-Transform
 */
class ITK_EXPORT TransformIOBase:public LightProcessObject
{
public:
  /** Standard class typedefs */
  typedef TransformIOBase      Self;
  typedef LightProcessObject   Superclass;
  typedef SmartPointer< Self > Pointer;
  /** Run-time type information (and related methods). */
  itkTypeMacro(TransformIOBase, Superclass);

  /** Transform types */
  typedef TransformBase TransformType;
  /** For writing, a const transform list gets passed in, for
   * reading, a non-const transform list is created from the file.
   */
  typedef TransformType::Pointer             TransformPointer;
  typedef std::list< TransformPointer >      TransformListType;
  typedef TransformType::ConstPointer        ConstTransformPointer;
  typedef std::list< ConstTransformPointer > ConstTransformListType;

  /** Set/Get the name of the file to be read. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

  /** Reads the data from disk into the memory buffer provided. */
  virtual void Read() = 0;

  /** Writes the transform list to disk. */
  virtual void Write() = 0;

  /** Determine the file type. Returns true if this TransformIO can read the
   * file specified. */
  virtual bool CanReadFile(const char *) = 0;

  /** Determine the file type. Returns true if this TransformIO can read the
   * file specified. */
  virtual bool CanWriteFile(const char *)  = 0;

  /** Get the list of transforms resulting from a file read */
  TransformListType & GetTransformList() { return m_ReadTransformList; }
  TransformListType & GetReadTransformList() { return m_ReadTransformList; }
  ConstTransformListType & GetWriteTransformList() { return m_WriteTransformList; }

  /** Set the list of transforms before writing */
  void SetTransformList(ConstTransformListType & transformList);

  /** Set the writer to append to the specified file */
  itkSetMacro(AppendMode, bool);
  itkGetConstMacro(AppendMode, bool);
  itkBooleanMacro(AppendMode);
protected:
  TransformIOBase();
  virtual ~TransformIOBase();
  void PrintSelf(std::ostream & os, Indent indent) const;

  void OpenStream(std::ofstream & out, bool binary);

  void CreateTransform(TransformPointer & ptr, const std::string & ClassName);

private:
  std::string            m_FileName;
  TransformListType      m_ReadTransformList;
  ConstTransformListType m_WriteTransformList;
  bool                   m_AppendMode;
};
} // end namespace itk
#endif // __itkTransformIOBase
