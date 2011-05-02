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
#ifndef __itkFEMObjectFactory_h
#define __itkFEMObjectFactory_h

#include "itkFastMutexLock.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <stdlib.h>

namespace itk {
namespace fem {

/**
 * \class FEMObjectFactory
 * \brief Create objects of derived classes by specifying a class ID.
 *
 * ObjectFactory is used to create various objects of derived classes by
 * specifying an ID of a derived class. Before the objects can be created by
 * object factory, you should first call the Register method for each class:
 *
 * \code
 * int ID_Derived=FEMObjectFactory<BaseClass>.Register( NewDerivedClass,
 *                                                      "NewDerivedClassName"
 *                                                    );
 * \endcode
 *
 *  - BaseClass: Base class from which the registered classes are derived
 *
 *  - ID_Derived: Integer returned by the Register function that
 *                specifies the derived class' ID. Class ID's
 *                are assigned automatically.
 *
 *  - NewDerivedClass: Function that creates a new instance of a derived
 *                     class and returns a pointer to it as a pointer
 *                     to the base class. Normally you should define this
 *                     function as:
 *
 * \code
 *                        - BaseClass* NewDerivedClass()
 *                              { return new DerivedClass; }
 * \endcode
 *
 *  - NewDerivedClassname: String with a name of a derived class.
 *
 * You should also put the NewDerivedClass function in an anonymous namespace
 * and make it static. So a complete registration of a derived class
 * looks like:
 *
 * \code
 *  - namespace { static Base* NewDerivedClass() { return new DerivedClass; } }
 *    const int OF_Derived=FEMObjectFactory<BaseClass>::Register( NewDerivedClass, "NewDerivedObjectName" );
 * \endcode
 * \ingroup ITK-FEM
 */
template<class T>
class FEMObjectFactory
{

  /**
   * Type that holds a pointer to function which can create a new object of a derived class.
   */
  typedef typename T::Pointer (*COF)();

  /**
   * Type that holds class name as a string.
   */
  typedef std::string StrClassName;

  /**
   * Type that holds an array of pairs of the COF pointers to functions and class names.
   */
  typedef std::vector<std::pair<COF,StrClassName> > COF_Array;
  typedef typename COF_Array::value_type COF_Array_value_type;

public:

  /**
   * Create a new object based on class identifier id and return a pointer to it.
   */
  static typename T::Pointer Create(int id)
    {
    return (Instance().cofs_[id].first)();
    }

  /**
   * Register the class with the factory. A pointer to a 'create'
   * function and class name as a string must be provided. Function
   * returns the newly assigned ID of the class, which can be later
   * used to create objects of that class.
   */
  static int Register(COF f, const char *str)
    {
    int clid=-1;
    Instance().m_MutexLock.Lock();
    Instance().cofs_.push_back( COF_Array_value_type(f,str) );
    clid = static_cast<int>( Instance().cofs_.size()-1 );
    Instance().m_MutexLock.Unlock();
    return clid;
    }

  /**
   * Return the name of the class (as a string) for the given ID.
   */
  static StrClassName ID2ClassName(int id)
    {
    return Instance().cofs_[id].second;
    }

  /**
   * Find the ID of the class with specified name (this is a slow function).
   * If you have to create many objects of the same class, obtain the class ID
   * with this function first and then use that ID with the Create member
   * function.
   */
  static int ClassName2ID(StrClassName str)
    {
    int j=0;
    for(typename COF_Array::const_iterator i=Instance().cofs_.begin(); i != Instance().cofs_.end(); i++)
      {
      if (i->second==str) return j;
      j++;
      }
    return -1;
    }

private:

  /**
   * Array that stores pairs of create functions and class names.
   */
  COF_Array cofs_;

  /**
   * Mutex lock to protect modification to the cofs_ array during
   * class registration.
   */
  mutable SimpleFastMutexLock m_MutexLock;

  /**
   * Private constructor. This class is implemented as a singleton, so we
   * don't allow anybody from outside to construct it.
   */
  FEMObjectFactory();

  /**
   * Private copy constructor.
   */
  FEMObjectFactory(const FEMObjectFactory&);

  /**
   * Private destructor.
   */
  ~FEMObjectFactory();

  /**
   * Access to the only instance of the FEMObjectFactory object.
   */
  inline static FEMObjectFactory& Instance();

  /**
   * Deletes the object in obj member. This function is
   * called when application finishes (atexit() function).
   */
  static void CleanUP();

  /**
   * Pointer to the only instance of the FEMObjectFactory class.
   */
  static FEMObjectFactory* obj;

private:
  /**
   * \class Dummy
   * \brief This class is defined in FEMObjectFactory just to get rid of some
            warnings about destructor being private in gcc.
   * \ingroup ITK-FEM
   */
  class Dummy {};

  /**
   * By defining a Dummy friend class, some warnings in gcc about destructor
   * being private are eliminated.
   */
  friend class Dummy;

};


template<class T>
FEMObjectFactory<T>* FEMObjectFactory<T>::obj = 0;

template<class T>
FEMObjectFactory<T>::FEMObjectFactory() {}

template<class T>
FEMObjectFactory<T>::FEMObjectFactory(const FEMObjectFactory<T>&) {}

template<class T>
FEMObjectFactory<T>::~FEMObjectFactory() {}

extern "C"
{
  typedef void(*c_void_cast)();
}
template<class T>
FEMObjectFactory<T>& FEMObjectFactory<T>::Instance()
{
  if (!obj)
    {
    /**
     * Create a new FEMObjectFactory object if we don't have it already.
     */
    obj=new FEMObjectFactory;

    /**
     * Make sure that the object that ws just created is also destroyed
     * when program finishes.
     */
    atexit(reinterpret_cast<c_void_cast>(&CleanUP));

    }

  /**
   * Return the actual FEMObjectFactory object
   */
  return *obj;
}

template<class T>
void FEMObjectFactory<T>::CleanUP() { delete obj; }

}} // end namespace itk::fem

#endif // #ifndef __itkFEMFEMObjectFactory_h
