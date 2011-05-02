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
/**
 * Copyright (c) 1996
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Moscow Center for SPARC Technology makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */

#ifndef __itk_hash_map_h
#define __itk_hash_map_h

#if defined( _MSC_VER )
#pragma warning ( disable : 4786 )
#endif

#if ( defined( __GNUC__ ) && ( ( ( __GNUC__ == 3 ) && ( __GNUC_MINOR__ >= 1 ) || ( __GNUC__ > 3 ) ) || ( ( __GNUC__ == \
                                                                                                           4 )         \
  && defined( __INTEL_COMPILER )  )  ) ) || ( defined( __IBMCPP__ ) && __IBMCPP__ >= 600 )
// Use this hash_map for GNU_C versions >= 3.1, IBMCPP >=600, or Intel compilers
// with GCCv4

// #if (defined(__GNUC__) && ((__GNUC__==4) && (__GNUC_MINOR__>=3) ) )
// in gcc 4.3 the header <ext/hash_map> was deprecated and replaced with
// <unordered_map>.
// note that this is still experimental in gcc 4.3
// #include <tr1/unordered_map>
// #else
#include <ext/hash_map>
// #endif

namespace itk
{
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
using __gnu_cxx::hash_multimap;
}

#else

#include "itk_hashtable.h"

// The following is required for CodeWarrior
#include "vcl_compiler.h"

namespace itk
{
#define VCL_IMPORT_CONTAINER_TYPEDEFS(super)             \
  typedef typename super::value_type      value_type;      \
  typedef typename super::reference       reference;       \
  typedef typename super::size_type       size_type;       \
  typedef typename super::const_reference const_reference; \
  typedef typename super::difference_type difference_type;

#define VCL_IMPORT_ITERATORS(super)              \
  typedef typename super::iterator       iterator; \
  typedef typename super::const_iterator const_iterator;

#define VCL_IMPORT_REVERSE_ITERATORS(super)                            \
  typedef typename super::const_reverse_iterator const_reverse_iterator; \
  typedef typename super::reverse_iterator       reverse_iterator;

template< class Key, class T, class HashFcn, class EqualKey, class Alloc >
class hash_map;
template< class Key, class T, class HashFcn, class EqualKey, class Alloc >
class hash_multimap;

template< class Key, class T, class HashFcn, class EqualKey, class Alloc >
bool operator==(const hash_map< Key, T, HashFcn, EqualKey, Alloc > &,
                const hash_map< Key, T, HashFcn, EqualKey, Alloc > &);

template< class Key, class T, class HashFcn, class EqualKey, class Alloc >
bool operator==(const hash_multimap< Key, T, HashFcn, EqualKey, Alloc > &,
                const hash_multimap< Key, T, HashFcn, EqualKey, Alloc > &);

/** \class hash_map
 * \brief Replacement for STL hash map because some systems do not support it,
 * or support it incorrectly.
 * \ingroup ITK-Common
 */
template< class Key, class T,
          VCL_DFL_TMPL_PARAM_STLDECL(HashFcn, hash< Key > ),
          VCL_DFL_TMPL_PARAM_STLDECL(EqualKey, std::equal_to< Key > ),
          VCL_DFL_TYPE_PARAM_STLDECL(Alloc, std::allocator< char > ) >
class hash_map
{
private:
  typedef std::select1st< std::pair< const Key, T > >                                   sel1st;
  typedef hashtable< std::pair< const Key, T >, Key, HashFcn, sel1st, EqualKey, Alloc > ht;
  typedef hash_map< Key, T, HashFcn, EqualKey, Alloc >                                  self;
public:
  VCL_IMPORT_CONTAINER_TYPEDEFS(ht)
  VCL_IMPORT_ITERATORS(ht)
  typedef typename ht::key_type      key_type;
  typedef typename ht::hasher        hasher;
  typedef typename ht::key_equal     key_equal;
  typedef T                          data_type;
  typedef typename ht::pointer       pointer;
  typedef typename ht::const_pointer const_pointer;
private:
  ht rep;
public:
  hasher hash_funct() const { return rep.hash_funct(); }
  key_equal key_eq() const { return rep.key_eq(); }
public:
  hash_map():rep( 100, hasher(), key_equal() ) {}
  hash_map(size_type n):rep( n, hasher(), key_equal() ) {}
  hash_map(size_type n, const hasher & hf):rep( n, hf, key_equal() ) {}
  hash_map(size_type n, const hasher & hf, const key_equal & eql):
    rep(n, hf, eql) {}

  hash_map(const value_type *f, const value_type *l):
    rep( 100, hasher(), key_equal() ) { rep.insert_unique(f, l); }
  hash_map(const value_type *f, const value_type *l, size_type n):
    rep( n, hasher(), key_equal() ) { rep.insert_unique(f, l); }
  hash_map(const value_type *f, const value_type *l, size_type n,
           const hasher & hf):
    rep( n, hf, key_equal() ) { rep.insert_unique(f, l); }
  hash_map(const value_type *f, const value_type *l, size_type n,
           const hasher & hf, const key_equal & eql):
    rep(n, hf, eql) { rep.insert_unique(f, l); }

  hash_map(const_iterator f, const_iterator l):
    rep( 100, hasher(), key_equal() ) { rep.insert_unique(f, l); }
  hash_map(const_iterator f, const_iterator l, size_type n):
    rep( n, hasher(), key_equal() ) { rep.insert_unique(f, l); }
  hash_map(const_iterator f, const_iterator l, size_type n,
           const hasher & hf):
    rep( n, hf, key_equal() ) { rep.insert_unique(f, l); }
  hash_map(const_iterator f, const_iterator l, size_type n,
           const hasher & hf, const key_equal & eql):
    rep(n, hf, eql) { rep.insert_unique(f, l); }
public:
  size_type size() const { return rep.size(); }
  size_type max_size() const { return rep.max_size(); }
  bool empty() const { return rep.empty(); }
  void swap(self & hs) { rep.swap(hs.rep); }

  friend bool operator== ITK_FRIEND_TEMPLATE_FUNCTION_ARGUMENT(self) (const self &, const self &);

  iterator begin() { return rep.begin(); }
  iterator end() { return rep.end(); }
  const_iterator begin() const { return rep.begin(); }
  const_iterator end() const { return rep.end(); }
public:
  std::pair< iterator, bool > insert(const value_type & obj)
  { return rep.insert_unique(obj); }
  void insert(const value_type *f, const value_type *l) { rep.insert_unique(f, l); }
  void insert(const_iterator f, const_iterator l) { rep.insert_unique(f, l); }
  std::pair< iterator, bool > insert_noresize(const value_type & obj)
  { return rep.insert_unique_noresize(obj); }

  iterator find(const key_type & key) { return rep.find(key); }
  const_iterator find(const key_type & key) const { return rep.find(key); }

  T & operator[](const key_type & key)
  {
    value_type val( key, T() );

    return rep.find_or_insert(val).second;
  }

  size_type count(const key_type & key) const { return rep.count(key); }

  std::pair< iterator, iterator > equal_range(const key_type & key)
  { return rep.equal_range(key); }
  std::pair< const_iterator, const_iterator > equal_range(const key_type & key) const
  { return rep.equal_range(key); }

  size_type erase(const key_type & key) { return rep.erase(key); }
  void erase(iterator it) { rep.erase(it); }
  void erase(iterator f, iterator l) { rep.erase(f, l); }
  void clear() { rep.clear(); }
public:
  void resize(size_type hint) { rep.resize(hint); }
  size_type bucket_count() const { return rep.bucket_count(); }
  size_type max_bucket_count() const { return rep.max_bucket_count(); }
  size_type elems_in_bucket(size_type n) const
  { return rep.elems_in_bucket(n); }
};

template< class Key, class T, VCL_DFL_TMPL_PARAM_STLDECL(HashFcn, hash< Key > ),
          VCL_DFL_TMPL_PARAM_STLDECL(EqualKey, std::equal_to< Key > ),
          VCL_DFL_TYPE_PARAM_STLDECL(Alloc, std::allocator< char > ) >
class hash_multimap
{
private:
  typedef hashtable< std::pair< const Key, T >, Key, HashFcn,
                     std::select1st< std::pair< const Key, T > >, EqualKey, Alloc > ht;
  typedef hash_multimap< Key, T, HashFcn, EqualKey, Alloc > self;
public:
  VCL_IMPORT_CONTAINER_TYPEDEFS(ht)
  VCL_IMPORT_ITERATORS(ht)
  typedef typename ht::key_type      key_type;
  typedef typename ht::hasher        hasher;
  typedef typename ht::key_equal     key_equal;
  typedef T                          data_type;
  typedef typename ht::pointer       pointer;
  typedef typename ht::const_pointer const_pointer;

  hasher hash_funct() const { return rep.hash_funct(); }
  key_equal key_eq() const { return rep.key_eq(); }
private:
  ht rep;
public:
  hash_multimap():rep( 100, hasher(), key_equal() ) {}
  hash_multimap(size_type n):rep( n, hasher(), key_equal() ) {}
  hash_multimap(size_type n, const hasher & hf):rep( n, hf, key_equal() ) {}
  hash_multimap(size_type n, const hasher & hf, const key_equal & eql):
    rep(n, hf, eql) {}

  hash_multimap(const value_type *f, const value_type *l):
    rep( 100, hasher(), key_equal() ) { rep.insert_equal(f, l); }
  hash_multimap(const value_type *f, const value_type *l, size_type n):
    rep( n, hasher(), key_equal() ) { rep.insert_equal(f, l); }
  hash_multimap(const value_type *f, const value_type *l, size_type n,
                const hasher & hf):
    rep( n, hf, key_equal() ) { rep.insert_equal(f, l); }
  hash_multimap(const value_type *f, const value_type *l, size_type n,
                const hasher & hf, const key_equal & eql):
    rep(n, hf, eql) { rep.insert_equal(f, l); }

  hash_multimap(const_iterator f, const_iterator l):
    rep( 100, hasher(), key_equal() ) { rep.insert_equal(f, l); }
  hash_multimap(const_iterator f, const_iterator l, size_type n):
    rep( n, hasher(), key_equal() ) { rep.insert_equal(f, l); }
  hash_multimap(const_iterator f, const_iterator l, size_type n,
                const hasher & hf):
    rep( n, hf, key_equal() ) { rep.insert_equal(f, l); }
  hash_multimap(const_iterator f, const_iterator l, size_type n,
                const hasher & hf, const key_equal & eql):
    rep(n, hf, eql) { rep.insert_equal(f, l); }
public:
  size_type size() const { return rep.size(); }
  size_type max_size() const { return rep.max_size(); }
  bool empty() const { return rep.empty(); }
  void swap(self & hs) { rep.swap(hs.rep); }

  friend bool operator== ITK_FRIEND_TEMPLATE_FUNCTION_ARGUMENT(self) (const self &, const self &);

  iterator begin() { return rep.begin(); }
  iterator end() { return rep.end(); }
  const_iterator begin() const { return rep.begin(); }
  const_iterator end() const { return rep.end(); }
public:
  iterator insert(const value_type & obj) { return rep.insert_equal(obj); }
  void insert(const value_type *f, const value_type *l) { rep.insert_equal(f, l); }
  void insert(const_iterator f, const_iterator l) { rep.insert_equal(f, l); }
  iterator insert_noresize(const value_type & obj)
  { return rep.insert_equal_noresize(obj); }

  iterator find(const key_type & key) { return rep.find(key); }
  const_iterator find(const key_type & key) const { return rep.find(key); }

  size_type count(const key_type & key) const { return rep.count(key); }

  std::pair< iterator, iterator > equal_range(const key_type & key)
  { return rep.equal_range(key); }
  std::pair< const_iterator, const_iterator > equal_range(const key_type & key) const
  { return rep.equal_range(key); }

  size_type erase(const key_type & key) { return rep.erase(key); }
  void erase(iterator it) { rep.erase(it); }
  void erase(iterator f, iterator l) { rep.erase(f, l); }
  void clear() { rep.clear(); }
public:
  void resize(size_type hint) { rep.resize(hint); }
  size_type bucket_count() const { return rep.bucket_count(); }
  size_type max_bucket_count() const { return rep.max_bucket_count(); }
  size_type elems_in_bucket(size_type n) const
  { return rep.elems_in_bucket(n); }
};

/** This method MUST NOT be declared "inline" because it a specialization of its template is
    declared as friend of a class. The hash_map class, in this case */
template< class Key, class T, class HashFcn, class EqualKey, class Alloc >
bool operator==(const hash_map< Key, T, HashFcn, EqualKey, Alloc > & hm1,
                const hash_map< Key, T, HashFcn, EqualKey, Alloc > & hm2)
{
  return hm1.rep == hm2.rep;
}

/** This method MUST NOT be declared "inline" because it a specialization of its template is
    declared as friend of a class. The hash_map class, in this case */
template< class Key, class T, class HashFcn, class EqualKey, class Alloc >
bool operator==(const hash_multimap< Key, T, HashFcn, EqualKey, Alloc > & hm1,
                const hash_multimap< Key, T, HashFcn, EqualKey, Alloc > & hm2)
{
  return hm1.rep == hm2.rep;
}

#define HASH_MAP_INSTANTIATE \
  extern "please include emulation/hash_map.txx instead"
} // end namespace itk

#endif

#endif // itk_emulation_hash_map_h
