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
#ifndef __itkMetaDataObject_hxx
#define __itkMetaDataObject_hxx

#include "itkMetaDataObject.h"

namespace itk
{
template< class MetaDataObjectType >
MetaDataObject< MetaDataObjectType >
::MetaDataObject(void)
{
  //Nothing to do, m_MetaDataObjectValue takes this types default value.
}

template< class MetaDataObjectType >
MetaDataObject< MetaDataObjectType >
::~MetaDataObject(void)
{
  //Nothing to do here.
}

template< class MetaDataObjectType >
MetaDataObject< MetaDataObjectType >
::MetaDataObject(const MetaDataObjectType InitializerValue):
  m_MetaDataObjectValue(InitializerValue)
{
  //Nothing to be done here
}

template< class MetaDataObjectType >
MetaDataObject< MetaDataObjectType >
::MetaDataObject(const MetaDataObject< MetaDataObjectType > & TemplateObject):
  Superclass(), m_MetaDataObjectValue(TemplateObject.m_MetaDataObjectValue)
{
  //Nothing to be done here
}

template< class MetaDataObjectType >
const char *
MetaDataObject< MetaDataObjectType >
::GetMetaDataObjectTypeName(void) const
{
  return typeid( MetaDataObjectType ).name();
}

template< class MetaDataObjectType >
const std::type_info &
MetaDataObject< MetaDataObjectType >
::GetMetaDataObjectTypeInfo(void) const
{
  return typeid( MetaDataObjectType );
}

template< class MetaDataObjectType >
const MetaDataObjectType &
MetaDataObject< MetaDataObjectType >
::GetMetaDataObjectValue(void) const
{
  return m_MetaDataObjectValue;
}

template< class MetaDataObjectType >
void
MetaDataObject< MetaDataObjectType >
::SetMetaDataObjectValue(const MetaDataObjectType & NewValue)
{
  m_MetaDataObjectValue = NewValue;
}

template< class MetaDataObjectType >
void
MetaDataObject< MetaDataObjectType >
::Print(std::ostream & os) const
{
  Superclass::Print(os);
}
} // end namespace itk

#endif
