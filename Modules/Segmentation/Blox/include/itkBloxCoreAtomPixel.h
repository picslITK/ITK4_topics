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
#ifndef __itkBloxCoreAtomPixel_h
#define __itkBloxCoreAtomPixel_h

#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/algo/vnl_generalized_eigensystem.h"

#include "itkObject.h"
#include "itkBloxCoreAtomItem.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkBloxPixel.h"

namespace itk
{
/**
 * \class BloxCoreAtomPixel
 * \brief Holds a linked list of itk::BloxCoreAtomItem's
 *
 * \ingroup ImageObjects
 *
 * \ingroup ITK-Blox
 */

template< unsigned int NDimensions >
class ITK_EXPORT BloxCoreAtomPixel:public BloxPixel<
    BloxCoreAtomItem< NDimensions > >
{
public:

  /** Standard class typedefs. */
  typedef BloxCoreAtomPixel                            Self;
  typedef BloxPixel< BloxCoreAtomItem< NDimensions > > Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** The type of core atom item we process. */
  typedef BloxCoreAtomItem< NDimensions > CoreAtomItemType;

  /** The type of boundary point item we process. */
  typedef BloxBoundaryPointItem< NDimensions > BPItemType;

  /** The type used to store the position of the BoundaryPointItem. */
  typedef Point< double, NDimensions > PositionType;

  /** The type of vector used to store the gradient of the BoundaryPointItem. */
  typedef CovariantVector< double, NDimensions > GradientType;

  /** VNL type used in eigenanalysis. */
  typedef vnl_vector_fixed< double, NDimensions > VectorType;

  /** Vector type used to store eigenvalues. */
  typedef vnl_vector_fixed< double, NDimensions > EigenvalueType;

  /** Matrix type used to store eigenvectors. */
  typedef vnl_matrix_fixed< double, NDimensions, NDimensions > EigenvectorType;

  /** Generalized matrix type used for several different tasks. */
  typedef vnl_matrix_fixed< double, NDimensions, NDimensions > MatrixType;

  /** Calculate and store the mean of core atom diameters. */
  double CalcMeanCoreAtomDiameter();

  /** Perform eigenanalysis on the population of core atoms
   *  stored in this pixel. */
  bool DoCoreAtomEigenanalysis();

  /** Perform eigenanalysis on the voted CMatrix */
  void DoVotedEigenanalysis();

  /** Get statements */
  double GetMeanCoreAtomDiameter() { return m_MeanCoreAtomDiameter; }
  double GetMeanCoreAtomIntensity()
  {
    return m_MeanCoreAtomIntensity; \
  }

  EigenvalueType GetEigenvalues() { return m_Eigenvalues; }
  EigenvalueType GetVotedEigenvalues() { return m_VotedEigenvalues; }
  EigenvectorType GetEigenvectors() { return m_Eigenvectors; }
  EigenvectorType GetVotedEigenvectors() { return m_VotedEigenvectors; }
  PositionType GetLocationSums() { return m_LocationSums; }
  double GetWeightSum() { return m_WeightSum; }

  /** Get the raw CMatrix (prior to voting) */
  MatrixType * GetRawCMatrixPointer() { return &m_RawCMatrix; }

  /** Collect a vote and update m_VotedCMatrix */
  void CollectVote(MatrixType *pMatrix, double strength, double count);

  /** Re-normalizes the voted CMatrix after all votes are cast */
  void NormalizeVotedCMatrix();

  /** Calculate location of the pixel based on core atoms voting for
   * it.
   */
  void CalcWeightedCoreAtomLocation(double weight_factor, Self *votingPixel);

  /** Returns the calculated voted location */
  PositionType GetVotedLocation();

  BloxCoreAtomPixel();
  ~BloxCoreAtomPixel();
private:

  /** Average (arithmetic mean) of core atom diameters stored in this pixel. */
  double m_MeanCoreAtomDiameter;

  /** The raw CMatrix - this is the matrix that we do eigen analysis on. */
  MatrixType m_RawCMatrix;

  /** The eigenvalues of the core atom population in this pixel
   * These are stored in increasing order of value (not absolute value) from
   * indices 0 to n, where n is the number of dimensions in the source image */
  EigenvalueType m_Eigenvalues;

  /** The eigenvectors of the core atom population in this pixel
   * Each eigenvector is a column of this matrix */
  EigenvectorType m_Eigenvectors;

  /** The CMatrix that collects votes cast by other blox. */
  MatrixType m_VotedCMatrix;

  /** Same as above, but calculated from the voted CMatrix */
  EigenvalueType m_VotedEigenvalues;

  /** Same as above, but calculated from the voted CMatrix */
  EigenvectorType m_VotedEigenvectors;

  /** The number of core atoms in all of the blox's that have voted for
   * this blox (its constituency) */
  double m_ConstituencySize;

  /** Used to compute the voted location of the core atom population */
  PositionType m_LocationSums;

  /** Used to compute the voted location of the core atom population */
  PositionType m_VotedLocation;

  /** Total weights used to compute voted location */
  double m_WeightSum;

  /** Total weights used to compute voted location */
  double m_MeanCoreAtomIntensity;
};
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_BloxCoreAtomPixel(_, EXPORT, TypeX, TypeY)     \
  namespace itk                                                     \
  {                                                                 \
  _( 1 ( class EXPORT BloxCoreAtomPixel< ITK_TEMPLATE_1 TypeX > ) ) \
  namespace Templates                                               \
  {                                                                 \
  typedef BloxCoreAtomPixel< ITK_TEMPLATE_1 TypeX >                 \
  BloxCoreAtomPixel##TypeY;                                       \
  }                                                                 \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkBloxCoreAtomPixel+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkBloxCoreAtomPixel.txx"
#endif

#endif
