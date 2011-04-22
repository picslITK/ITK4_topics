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
#ifndef __itkVTKImageIO_h
#define __itkVTKImageIO_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <fstream>
#include "itkStreamingImageIOBase.h"

namespace itk
{
/** \class VTKImageIO
 *
 *  \brief ImageIO class for reading VTK images
 *
 * This implementation was taken fron the Insight Joural:
 * http://hdl.handle.net/10380/3171
 *
 * \ingroup IOFilters
 *
 * \ingroup ITK-IO-VTK
 */
class ITK_EXPORT VTKImageIO:
  public StreamingImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef VTKImageIO                 Self;
  typedef StreamingImageIOBase       Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VTKImageIO, StreamingImageIOBase);

// see super class for documentation
  //
  // overidden to return true only when supported
  virtual bool CanStreamWrite(void);

  // see super class for documentation
  //
  // overidden to return true only when supported
  virtual bool CanStreamRead(void);


  /*-------- This part of the interface deals with reading data. ------ */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanReadFile(const char *);

  /** Set the spacing and dimesion information for the current filename. */
  virtual void ReadImageInformation();

  /** Reads the data from disk into the memory buffer provided. */
  virtual void Read(void *buffer);

  /*-------- This part of the interfaces deals with writing data. ----- */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanWriteFile(const char *);

  /** Writes the spacing and dimentions of the image.
   * Assumes SetFileName has been called with a valid file name. */
  virtual void WriteImageInformation() {}

  /** Writes the data to disk from the memory buffer provided. Make sure
   * that the IORegion has been set properly. */
  virtual void Write(const void *buffer);

  /** returns the header size, if it is unknown it will return 0 */
  virtual SizeType GetHeaderSize() const { return this->m_HeaderSize; }

protected:
  VTKImageIO();
  ~VTKImageIO();

  void PrintSelf(std::ostream & os, Indent indent) const;

  void InternalReadImageInformation(std::ifstream & file);

  void WriteImageInformation(const void *buffer);

  void ReadHeaderSize(std::ifstream & file);

  /** Convenient method to read a buffer as ASCII text. */
  virtual void ReadBufferAsASCII(std::istream & os, void *buffer,
                         IOComponentType ctype,
                         SizeType numberOfBytesToBeRead);

  /** Convenient method to write a buffer as ASCII text. */
  virtual void WriteBufferAsASCII(std::ostream & os, const void *buffer,
                          IOComponentType ctype,
                          SizeType numberOfBytesToWrite);

  /** We have a special method to read symmetric second rank tensors because
   * the VTK file format expands the symmetry and only supports 3D tensors. */
  virtual void ReadSymmetricTensorBufferAsBinary(std::istream& os,
    void *buffer,
    StreamingImageIOBase::SizeType num);

  /** We have a special method to write symmetric second rank tensors because
   * the VTK file format expands the symmetry and only supports 3D tensors. */
  virtual void WriteSymmetricTensorBufferAsBinary(std::ostream& os,
    const void *buffer,
    StreamingImageIOBase::SizeType num);

private:
  VTKImageIO(const Self &);    //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  void SetPixelTypeFromString(const std::string & pixelType);

  SizeType m_HeaderSize;
};
} // end namespace itk

#endif // __itkVTKImageIO_h
