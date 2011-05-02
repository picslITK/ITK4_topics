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
#ifndef __itkPolygonGroupSpatialObjectXMLFile_h
#define __itkPolygonGroupSpatialObjectXMLFile_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkPolygonGroupSpatialObject.h"
#include "itkXMLFile.h"
namespace itk
{
/* 3D Polygon Groups only ones that make sense for this data type */
typedef PolygonGroupSpatialObject< 3 > PGroupSpatialObjectType;

/** \class PolygonGroupSpatialObjectXMLFileReader
 *
 * Reads an XML-format file containing a list of polygons, and
 * creates a corresponding PolygonGroupSpatialObject
 * \ingroup ITK-IO-SpatialObjects
 */
class ITK_EXPORT PolygonGroupSpatialObjectXMLFileReader:
  public XMLReader< PGroupSpatialObjectType >
{
public:
  /** Standard typedefs */
  typedef PolygonGroupSpatialObjectXMLFileReader Self;
  typedef XMLReader< PGroupSpatialObjectType >   Superclass;
  typedef SmartPointer< Self >                   Pointer;

  typedef PGroupSpatialObjectType   PolygonGroupType;
  typedef PolygonSpatialObject< 3 > PolygonSpatialObjectType;
  typedef SpatialObjectPoint< 3 >   PointType;
  typedef std::vector< PointType >  PointListType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(PolygonGroupSpatialObjectXMLFileReader, XMLReader);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
public:
  /** Determine if a file can be read */
  virtual int CanReadFile(const char *name);

protected:
  PolygonGroupSpatialObjectXMLFileReader() {}
  virtual ~PolygonGroupSpatialObjectXMLFileReader() {}

  virtual void StartElement(const char *name, const char **atts);

  virtual void EndElement(const char *name);

  virtual void CharacterDataHandler(const char *inData, int inLength);

private:
  PolygonGroupSpatialObjectXMLFileReader(const Self &); //purposely not
                                                        // implemented
  void operator=(const Self &);                         //purposely not
                                                        // implemented

  PGroupSpatialObjectType::Pointer  m_PGroup;
  PolygonSpatialObjectType::Pointer m_CurPoly;
  PointListType                     m_CurPointList;
  std::string                       m_CurCharacterData;
};

/** \class PolygonGroupSpatialObjectXMLFileWriter
 *
 * Writes an XML-format file containing a list of polygons,
 * based on a PolygonGroupSpatialObject.
 * \ingroup ITK-IO-SpatialObjects
 */
class ITK_EXPORT PolygonGroupSpatialObjectXMLFileWriter:
  public XMLWriterBase< PGroupSpatialObjectType >
{
public:
  /** standard typedefs */
  typedef XMLWriterBase< PGroupSpatialObjectType > Superclass;
  typedef PolygonGroupSpatialObjectXMLFileWriter   Self;
  typedef SmartPointer< Self >                     Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PolygonGroupSpatialObjectXMLFileWriter,
               XMLWriterBase< PGroupSpatialObjectType > );
  typedef PGroupSpatialObjectType   PolygonGroupType;
  typedef PolygonSpatialObject< 3 > PolygonSpatialObjectType;
  /** Test whether a file is writable. */
  virtual int CanWriteFile(const char *name);

  /** Actually write out the file in question */
  virtual int WriteFile();

protected:
  PolygonGroupSpatialObjectXMLFileWriter() {}
  virtual ~PolygonGroupSpatialObjectXMLFileWriter() {}
private:
  PolygonGroupSpatialObjectXMLFileWriter(const Self &); //purposely not
                                                        // implemented
  void operator=(const Self &);                         //purposely not
                                                        // implemented
};
}
#endif
