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
#ifndef __itkPhilipsRECImageIOFactory_h
#define __itkPhilipsRECImageIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
/** \class PhilipsRECImageIOFactory
 * \brief Create instances of PhilipsRECImageIO objects using an object factory.
 *
 * \author Don C. Bigler
 *         The Pennsylvania State University 2005
 *
 * This implementation was contributed as a paper to the Insight Journal
 * http://insight-journal.org/midas/handle.php?handle=1926/1381
 *
 * \ingroup ITK-Review
 */
class ITK_EXPORT PhilipsRECImageIOFactory:public ObjectFactoryBase
{
public:
  /** Standard class typedefs. */
  typedef PhilipsRECImageIOFactory   Self;
  typedef ObjectFactoryBase          Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Class methods used to interface with the registered factories. */
  virtual const char * GetITKSourceVersion(void) const;

  virtual const char * GetDescription(void) const;

  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);
  static PhilipsRECImageIOFactory * FactoryNew()
  {
    return new PhilipsRECImageIOFactory;
  }

  /** Run-time type information (and related methods). */
  itkTypeMacro(PhilipsRECImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
  {
    PhilipsRECImageIOFactory::Pointer factory =
      PhilipsRECImageIOFactory::New();

    ObjectFactoryBase::RegisterFactory(factory);
  }

protected:
  PhilipsRECImageIOFactory();
  ~PhilipsRECImageIOFactory();
private:
  PhilipsRECImageIOFactory(const Self &); //purposely not implemented
  void operator=(const Self &);           //purposely not implemented
};
} // end namespace itk

#endif
