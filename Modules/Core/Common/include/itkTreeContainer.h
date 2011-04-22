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
#ifndef __itkTreeContainer_h
#define __itkTreeContainer_h

#include "itkTreeContainerBase.h"
#include "itkPreOrderTreeIterator.h"

namespace itk
{
/** \class TreeContainer
 *  \brief TreeContainer class
 *
 * This class derives from the TreeContainerBase class.
 *
 * The class is templated over the type of the elements.
 *
 * Template parameters for class TreeContainer:
 *
 * - TValueType = Element type stored at each location in the Tree.
 *
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< class TValueType >
class TreeContainer:public TreeContainerBase< TValueType >
{
public:

  /** Standard typedefs */
  typedef TreeContainerBase< TValueType > Superclass;
  typedef TreeContainer< TValueType >     Self;
  typedef SmartPointer< Self >            Pointer;
  typedef SmartPointer< const Self >      ConstPointer;
  typedef TValueType                      ValueType;
  typedef TreeNode< ValueType >           TreeNodeType;

  /** Iterators typedef */
  typedef TreeIteratorBase< Self >     IteratorType;
  typedef PreOrderTreeIterator< Self > PreOrderIteratorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TreeContainer, TreeContainerBase);

  /** Constructor */
  TreeContainer(int defaultChildrenCount);

  /** Constructor */
  TreeContainer(TreeContainer< TValueType > & tree);

  /** Set the root as an element */
  virtual bool SetRoot(const TValueType element);

  /** The the root as an iterator position */
  bool SetRoot(IteratorType & pos);

  /** Set the root as a tree node */
  virtual bool SetRoot(TreeNode< TValueType > *node);

  /** Return true if the element is in the tree */
  bool Contains(const TValueType element);

  /** Return the number of elements in the tree */
  int Count() const;

  /** Return true if the element is a leaf */
  bool IsLeaf(const TValueType element);

  /** Return true if the element is a root */
  bool IsRoot(const TValueType element);

  /** Clear the tree */
  bool Clear();

  /** operator equal */
  bool operator==(TreeContainer< TValueType > & tree);

  /** Swap the iterators */
  bool Swap(IteratorType & v, IteratorType & w);

  /** Get the root */
  const TreeNodeType * GetRoot() const { return m_Root.GetPointer(); }

  /** Add a child to a given parent using values */
  bool Add(const TValueType child, const TValueType parent);

  /** Get node given a value */
  const TreeNodeType * GetNode(TValueType val) const;

protected:

  TreeContainer();
  virtual ~TreeContainer();

  typename TreeNodeType::Pointer m_Root;

  int m_DefaultChildrenCount;

  void PrintSelf(std::ostream & os, Indent indent) const;
};
} // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_TreeContainer(_, EXPORT, TypeX, TypeY)     \
  namespace itk                                                 \
  {                                                             \
  _( 1 ( class EXPORT TreeContainer< ITK_TEMPLATE_1 TypeX > ) ) \
  namespace Templates                                           \
  {                                                             \
  typedef TreeContainer< ITK_TEMPLATE_1 TypeX >                 \
  TreeContainer##TypeY;                                       \
  }                                                             \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkTreeContainer+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkTreeContainer.txx"
#endif

#endif
