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
#ifndef __itkMetaVesselTubeConverter_h
#define __itkMetaVesselTubeConverter_h

#include "metaVesselTube.h"
#include "itkMetaConverterBase.h"
#include "itkVesselTubeSpatialObject.h"

/** \class MetaVesselTubeConverter
 *  \brief converts between MetaObject<->SpatialObject
 *  \sa MetaConverterBase
 *  \ingroup ITK-SpatialObjects
 */
namespace itk
{
/** \class MetaVesselTubeConverter
 * \brief This is the MetaVesselTubeConverter class.
 *
 * \ingroup ITK-SpatialObjects
 */
template< unsigned int NDimensions = 3 >
class ITK_EXPORT MetaVesselTubeConverter :
    public MetaConverterBase< NDimensions >
{
public:
  /** Standard class typedefs */
  typedef MetaVesselTubeConverter          Self;
  typedef MetaConverterBase< NDimensions > Superclass;
  typedef SmartPointer< Self >             Pointer;
  typedef SmartPointer< const Self >       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MetaVesselTubeConverter, MetaConverterBase);

  typedef typename Superclass::SpatialObjectType SpatialObjectType;
  typedef typename SpatialObjectType::Pointer    SpatialObjectPointer;
  typedef typename Superclass::MetaObjectType    MetaObjectType;

  /** Specific class types for conversion */
  typedef VesselTubeSpatialObject<NDimensions>               VesselTubeSpatialObjectType;
  typedef typename VesselTubeSpatialObjectType::Pointer      VesselTubeSpatialObjectPointer;
  typedef typename VesselTubeSpatialObjectType::ConstPointer VesselTubeSpatialObjectConstPointer;
  typedef MetaVesselTube                                     VesselTubeMetaObjectType;

  /** Convert the MetaObject to Spatial Object */
  virtual SpatialObjectPointer MetaObjectToSpatialObject(const MetaObjectType *mo);

  /** Convert the SpatialObject to MetaObject */
  virtual MetaObjectType *SpatialObjectToMetaObject(const SpatialObjectType *spatialObject);

protected:
  /** Create the specific MetaObject for this class */
  virtual MetaObjectType *CreateMetaObject();

  MetaVesselTubeConverter();
  ~MetaVesselTubeConverter() {}

private:
  MetaVesselTubeConverter(const Self &);   //purposely not implemented
  void operator=(const Self &);       //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
  #include "itkMetaVesselTubeConverter.txx"
#endif

#endif
