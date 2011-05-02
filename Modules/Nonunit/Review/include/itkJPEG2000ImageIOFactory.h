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
#ifndef __itkJPEG2000ImageIOFactory_h
#define __itkJPEG2000ImageIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
/** \class JPEG2000ImageIOFactory
 * \brief Supports for the JPEG2000 file format based on openjpeg
 *
 *  JPEG2000 offers a large collection of interesting features including:
 *  compression (lossless and lossy), streaming, multi-channel images.
 *
 * \ingroup ITK-Review
 */
class ITK_EXPORT JPEG2000ImageIOFactory:public ObjectFactoryBase
{
public:
  /** Standard class typedefs. */
  typedef JPEG2000ImageIOFactory     Self;
  typedef ObjectFactoryBase          Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Class methods used to interface with the registered factories. */
  virtual const char * GetITKSourceVersion() const;

  virtual const char * GetDescription() const;

  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);
  static JPEG2000ImageIOFactory * FactoryNew() { return new JPEG2000ImageIOFactory; }

  /** Run-time type information (and related methods). */
  itkTypeMacro(JPEG2000ImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory()
  {
    JPEG2000ImageIOFactory::Pointer metaFactory = JPEG2000ImageIOFactory::New();

    ObjectFactoryBase::RegisterFactory(metaFactory);
  }

protected:
  JPEG2000ImageIOFactory();
  ~JPEG2000ImageIOFactory();
private:
  JPEG2000ImageIOFactory(const Self &); //purposely not implemented
  void operator=(const Self &);         //purposely not implemented
};
} // end namespace itk

#endif
