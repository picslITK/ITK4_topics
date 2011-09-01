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

#ifndef __itkOBJMeshIO_h
#define __itkOBJMeshIO_h

#include "itkMeshIOBase.h"

#include <fstream>

namespace itk
{
/** \class OBJMeshIO
 * \brief This class defines how to read and write Object file format.
 * \ingroup IOFilters
 * \ingroup ITKIOMesh
 */

class ITK_EXPORT OBJMeshIO:public MeshIOBase
{
public:
  /** Standard class typedefs. */
  typedef OBJMeshIO                  Self;
  typedef MeshIOBase                 Superclass;
  typedef SmartPointer< const Self > ConstPointer;
  typedef SmartPointer< Self >       Pointer;

  typedef Superclass::SizeValueType    SizeValueType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OBJMeshIO, MeshIOBase);

  /*-------- This part of the interfaces deals with reading data. ----- */

  /** Determine if the file can be read with this MeshIO implementation.
  * \param FileNameToRead The name of the file to test for reading.
  * \post Sets classes MeshIOBase::m_FileName variable to be FileNameToWrite
  * \return Returns true if this MeshIO can read the file specified.
  */
  virtual bool CanReadFile(const char *FileNameToRead);

  /** Set the spacing and dimension information for the set filename. */
  virtual void ReadMeshInformation();

  /** Reads the data from disk into the memory buffer provided. */
  virtual void ReadPoints(void *buffer);

  virtual void ReadCells(void *buffer);

  virtual void ReadPointData(void *buffer);

  virtual void ReadCellData(void *buffer);

  /*-------- This part of the interfaces deals with writing data. ----- */

  /** Determine if the file can be written with this MeshIO implementation.
   * \param FileNameToWrite The name of the file to test for writing.
   * \post Sets classes MeshIOBase::m_FileName variable to be FileNameToWrite
   * \return Returns true if this MeshIO can write the file specified.
   */
  virtual bool CanWriteFile(const char *FileNameToWrite);

  /** Set the spacing and dimension information for the set filename. */
  virtual void WriteMeshInformation();

  /** Writes the data to disk from the memory buffer provided. Make sure
   * that the IORegions has been set properly. */
  virtual void WritePoints(void *buffer);

  virtual void WriteCells(void *buffer);

  virtual void WritePointData(void *buffer);

  virtual void WriteCellData(void *buffer);

  virtual void Write();

protected:
  /** Write points to output stream */
  template< typename T >
  void WritePoints(T *buffer, std::ofstream & outputFile)
  {
    SizeValueType index = itk::NumericTraits< SizeValueType >::Zero;

    for ( SizeValueType ii = 0; ii < this->m_NumberOfPoints; ii++ )
      {
      outputFile << "v ";
      for ( unsigned int jj = 0; jj < this->m_PointDimension; jj++ )
        {
        outputFile << buffer[index++] << "  ";
        }
      outputFile << '\n';
      }
  }

  template< typename T >
  void WriteCells(T *buffer, std::ofstream & outputFile)
  {
    SizeValueType index = itk::NumericTraits< SizeValueType >::Zero;

    for ( SizeValueType ii = 0; ii < this->m_NumberOfCells; ii++ )
      {
      outputFile << "f ";
      index++;
      unsigned int numberOfCellPoints = static_cast< unsigned int >( buffer[index++] );

      for ( unsigned int jj = 0; jj < numberOfCellPoints; jj++ )
        {
        outputFile << buffer[index++] + 1 << "  ";
        }
      outputFile << '\n';
      }
  }

  /** Write point data to output stream */
  template< typename T >
  void WritePointData(T *buffer, std::ofstream & outputFile)
  {
    SizeValueType index = itk::NumericTraits< SizeValueType >::Zero;

    for ( SizeValueType ii = 0; ii < this->m_NumberOfPointPixels; ii++ )
      {
      outputFile << "vn ";
      for ( unsigned int jj = 0; jj < this->m_PointDimension; jj++ )
        {
        outputFile << buffer[index++] << "  ";
        }

      outputFile << '\n';
      }
  }

protected:
  OBJMeshIO();
  virtual ~OBJMeshIO(){}

  void PrintSelf(std::ostream & os, Indent indent) const;

  void OpenFile();

  void CloseFile();

private:
  OBJMeshIO(const Self &);      // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  std::ifstream  m_InputFile;
  std::streampos m_PointsStartPosition;  // file position for points rlative to
                                         // std::ios::beg
};
} // end namespace itk

#endif
