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
#ifndef __itkCorrespondingMedialNodeClique_h
#define __itkCorrespondingMedialNodeClique_h

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkFixedArray.h"
#include "itkBloxCoreAtomPixel.h"

#include <list>

namespace itk
{
/**
 * \class CorrespondingMedialNodeClique
 * \brief CorrespondingMedialNodeClique is an item stored
 * in CorrespondingNodeList. Specifically it is stored in
 * corresponding node lists and contain pointers to a
 * set of medial nodes (cliques).
 *
 *
 * \ingroup ITK-Blox
 */

template< unsigned int VImageDimension, unsigned int VCliqueSize >
class CorrespondingMedialNodeClique
{
public:

  /** Medial node typedef. */
  typedef BloxCoreAtomPixel< VImageDimension > ItemType;

  /** A vector of pointers to medial nodes. */
  std::vector< ItemType * > m_ItemPointer;

  /** Set the pointer to medial nodes. */
  void SetNodePointer(ItemType *itemPointer, unsigned int index)
  { m_ItemPointer[index] = itemPointer; }

  /** Coordinate of node in clique in physical space. */
  typedef FixedArray< vnl_vector_fixed< double, VImageDimension >, VCliqueSize >
  CoordinateType;

  /** Center mass of node clique in physical space. */
  typedef vnl_vector_fixed< double, VCliqueSize > CenterOfMassType;

  /** Transform matrix. */
  typedef vnl_matrix_fixed< double, VImageDimension + 1, VImageDimension + 1 >
  TransformMatrixType;

  /** Set and get the coordinates of the nodes in the clique. */
  void SetNodeCoordinates(CoordinateType *coordinates)
  { m_NodeCoordinates = coordinates; }
  CoordinateType * GetNodeCoordinates() { return m_NodeCoordinates; }

  /** Set and get the center of mass of the clique. */
  void SetCenterOfMass(CenterOfMassType *centerOfMass)
  { m_CenterOfMass = centerOfMass; }
  CenterOfMassType * GetCenterOfMass() { return m_CenterOfMass; }

  /** Set and get the transform matrix. */
  void SetTransformMatrix(TransformMatrixType *transformMatrix)
  { m_TransformMatrix = transformMatrix; }
  TransformMatrixType * GetTransformMatrix() { return m_TransformMatrix; }

  /** Set and get the node index. */
  void SetNodeIndex(int index, int nodeIndex)
  { m_NodeIndex[index] = nodeIndex; }
  int GetNodeIndex(int index) { return m_NodeIndex[index]; }

  /** Set and get the correspondence value. */
  void SetCorrespondenceValue(int index, float correspondenceValue)
  { m_CorrespondenceValue[index] = correspondenceValue; }
  float GetCorrespondenceValue(int index)
  { return m_CorrespondenceValue[index]; }

  CorrespondingMedialNodeClique();
  ~CorrespondingMedialNodeClique();
private:

  /** Coordinate of the nodes of the clique. */
  CoordinateType *m_NodeCoordinates;

  /** Center of mass of the node clique. */
  CenterOfMassType *m_CenterOfMass;

  /** Transform matrix. */
  TransformMatrixType *m_TransformMatrix;

  /** Index of medial nodes in this clique. */
  int m_NodeIndex[VCliqueSize];

  /** Store the correspondence value. */
  float m_CorrespondenceValue[VCliqueSize];

  /** Average distance between nodes of clique in physical space. */
  double m_AverageDistance;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCorrespondingMedialNodeClique.txx"
#endif

#endif
