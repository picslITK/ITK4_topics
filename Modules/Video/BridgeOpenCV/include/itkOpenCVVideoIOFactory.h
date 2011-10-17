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
#ifndef __itkOpenCVVideoIOFactory_h
#define __itkOpenCVVideoIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkVideoIOBase.h"

namespace itk
{
/** \class OpenCVVideoIOFactory
 * \brief Create instances of OpenCVVideoIO objects using an object factory.
 *
 * \ingroup Video-IO-OpenCV
 */
class ITK_EXPORT OpenCVVideoIOFactory:public ObjectFactoryBase
{
public:
  /** Standard class typedefs. */
  typedef OpenCVVideoIOFactory       Self;
  typedef ObjectFactoryBase          Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Class methods used to interface with the registered factories. */
  virtual const char * GetITKSourceVersion(void) const;

  virtual const char * GetDescription(void) const;

  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OpenCVVideoIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
  {
    OpenCVVideoIOFactory::Pointer OpenCVFactory = OpenCVVideoIOFactory::New();

    ObjectFactoryBase::RegisterFactoryInternal(OpenCVFactory);
  }

protected:
  OpenCVVideoIOFactory();
  ~OpenCVVideoIOFactory();
private:
  OpenCVVideoIOFactory(const Self &); //purposely not implemented
  void operator=(const Self &);    //purposely not implemented
};
} // end namespace itk

#endif
