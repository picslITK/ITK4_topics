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
#ifndef __itkPostOrderTreeIterator_h
#define __itkPostOrderTreeIterator_h

#include "itkTreeIteratorBase.h"

namespace itk
{
template< class TTreeType >
class PostOrderTreeIterator:public TreeIteratorBase< TTreeType >
{
public:

  /** Typedefs */
  typedef PostOrderTreeIterator             Self;
  typedef TreeIteratorBase< TTreeType >     Superclass;
  typedef TTreeType                         TreeType;
  typedef typename TTreeType::ValueType     ValueType;
  typedef typename Superclass::TreeNodeType TreeNodeType;
  typedef typename Superclass::NodeType     NodeType;

  /** Constructor */
  PostOrderTreeIterator(TreeType *tree);

  /** Get the type of the iterator */
  NodeType GetType() const;

  /** Clone function */
  TreeIteratorBase< TTreeType > * Clone();

protected:
  /** Return the next node */
  const ValueType & Next();

  /** Return true if the next node exists */
  bool HasNext() const;

protected:

  const TreeNodeType * FindNextNode() const;

  const TreeNodeType * FindMostRightLeaf(TreeNodeType *node) const;

  const TreeNodeType * FindSister(TreeNodeType *node) const;
};

/** Constructor */
template< class TTreeType >
PostOrderTreeIterator< TTreeType >::PostOrderTreeIterator(TTreeType *tree):
  TreeIteratorBase< TTreeType >(tree, NULL)
{
  if ( tree->GetRoot() == 0 )
    {
    this->m_Begin = 0;
    }
  else
    {
    const TreeNodeType *root = dynamic_cast<const TreeNodeType *>(tree->GetRoot());
    if(root == 0)
      {
      itkGenericExceptionMacro(<< "Can't downcast root node to TreeNodeType *");
      }
    this->m_Position = const_cast<TreeNodeType *>(root);
    this->m_Position = const_cast< TreeNodeType * >( FindMostRightLeaf(this->m_Position) );
    this->m_Begin = this->m_Position;
    }
}

/** Return the type of the iterator */
template< class TTreeType >
typename PostOrderTreeIterator< TTreeType >::NodeType
PostOrderTreeIterator< TTreeType >::GetType() const
{
  return TreeIteratorBase< TTreeType >::POSTORDER;
}

/** Return true if the next node exists */
template< class TTreeType >
bool
PostOrderTreeIterator< TTreeType >::HasNext() const
{
  if ( const_cast< TreeNodeType * >( FindNextNode() ) != NULL )
    {
    return true;
    }
  return false;
}

/** Go to the next node */
template< class TTreeType >
const typename PostOrderTreeIterator< TTreeType >::ValueType &
PostOrderTreeIterator< TTreeType >::Next()
{
  this->m_Position = const_cast< TreeNodeType * >( FindNextNode() );
  return this->m_Position->Get();
}

/** Find the next node */
template< class TTreeType >
const typename PostOrderTreeIterator< TTreeType >::TreeNodeType *
PostOrderTreeIterator< TTreeType >::FindNextNode() const
{
  if ( this->m_Position == NULL || this->m_Position == this->m_Root )
    {
    return NULL;
    }
  TreeNodeType *sister = const_cast< TreeNodeType * >( FindSister(this->m_Position) );

  if ( sister != NULL )
    {
    return FindMostRightLeaf(sister);
    }
  if(this->m_Position->GetParent() == 0)
    {
    return 0;
    }
  TreeNodeType *rval = dynamic_cast<TreeNodeType *>(this->m_Position->GetParent());
  if(rval == 0)
    {
      itkGenericExceptionMacro(<< "Can't downcast to TreeNodeType *");
    }
  return rval;
}

/** Find the sister node */
template< class TTreeType >
const typename PostOrderTreeIterator< TTreeType >::TreeNodeType *
PostOrderTreeIterator< TTreeType >::FindSister(TreeNodeType *node) const
{
  if ( !node->HasParent() )
    {
    return NULL;
    }

  TreeNodeType *parent = dynamic_cast<TreeNodeType *>(node->GetParent());
  if(parent == 0)
    {
    itkGenericExceptionMacro(<< "Can't downcast to TreeNodeType *");
    }

  int           childPosition = parent->ChildPosition(node);
  int           lastChildPosition = parent->CountChildren() - 1;

  while ( childPosition < lastChildPosition )
    {
    if(parent->GetChild(childPosition + 1) == 0)
      {
      childPosition++;
      }
    else
      {
      TreeNodeType *sister = dynamic_cast<TreeNodeType *>(parent->GetChild(childPosition + 1));
      if ( sister == 0)
      {
      itkGenericExceptionMacro(<< "Can't downcast to TreeNodeType *");
      }
      return sister;
      }
    }
  return NULL;
}

/** Find the most right leaf */
template< class TTreeType >
const typename PostOrderTreeIterator< TTreeType >::TreeNodeType *
PostOrderTreeIterator< TTreeType >::FindMostRightLeaf(TreeNodeType *node) const
{
  while ( node->HasChildren() )
    {
    TreeNodeType *helpNode;
    int           childCount = node->CountChildren();
    int           i = 0;

    do
      {
      if(node->GetChild(i) == 0)
        {
        helpNode = 0;
        }
      else
        {
        helpNode = dynamic_cast<TreeNodeType *>(node->GetChild(i));
       if(helpNode == 0)
          {
          itkGenericExceptionMacro(<< "Can't downcast to TreeNodeType *");
          }
        }
      i++;
      }
    while ( helpNode == NULL && i < childCount );

    if ( helpNode == NULL )
      {
      return node;
      }
    node = helpNode;
    }
  return node;
}

/** Clone function */
template< class TTreeType >
TreeIteratorBase< TTreeType > *PostOrderTreeIterator< TTreeType >::Clone()
{
  PostOrderTreeIterator< TTreeType > *clone =
    new PostOrderTreeIterator< TTreeType >( const_cast< TTreeType * >( this->m_Tree ) );
  *clone = *this;
  return clone;
}
} // end namespace itk

#endif
