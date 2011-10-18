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
#ifndef __itkExpectationBasedPointSetMetric_h
#define __itkExpectationBasedPointSetMetric_h

#include "itkPointSetToPointSetMetric.h"
#include "itkCovariantVector.h"
#include "itkPointSet.h"
#include "itkImage.h"

namespace itk
{
/**
 * \class ExpectationBasedPointSetMetric
 * \brief Computes an expectation-based metric between two point sets.
 *
 * This information-theoretic point set measure models each point set
 * as a sum of Gaussians.  To speed up computation, evaluation of the local
 * value/derivative is done in a user-specified neighborhood using the k-d
 * tree constructed in the superclass.
 *
 *  Reference:
 *    Pluta J, Avants BB, Glynn S, Awate S, Gee JC, Detre JA,
 *    "Appearance and incomplete label matching for diffeomorphic template
 *     "based hippocampus segmentation", Hippocampus, 2009 Jun; 19(6):565-71.
 *
 * \ingroup ITKRegistrationRefactoring
 */
template<class TFixedPointSet, class TMovingPointSet = TFixedPointSet>
class ITK_EXPORT ExpectationBasedPointSetMetric:
  public PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
{
public:

  /** Standard class typedefs. */
  typedef ExpectationBasedPointSetMetric                       Self;
  typedef PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ExpectationBasedPointSetMetric, PointSetToPointSetMetric );

   /** Types transferred from the base class */
  itkSuperclassTraitMacro( MeasureType );
  itkSuperclassTraitMacro( DerivativeType );
  itkSuperclassTraitMacro( LocalDerivativeType );
  itkSuperclassTraitMacro( PointType );
  itkSuperclassTraitMacro( CoordRepType );
  itkSuperclassTraitMacro( PointIdentifier );
  itkSuperclassTraitMacro( NeighborsIdentifierType );

  /**
   * Calculates the local metric value for a single point.
   */
  virtual MeasureType GetLocalNeighborhoodValue( const PointType ) const;

  /**
   * Calculates the local value and derivative for a single point.
   */
  virtual void GetLocalNeighborhoodValueAndDerivative( const PointType,
    MeasureType &, LocalDerivativeType & ) const;

protected:
  ExpectationBasedPointSetMetric();
  virtual ~ExpectationBasedPointSetMetric() {}

  /** PrintSelf funtion */
  void PrintSelf( std::ostream & os, Indent indent ) const;

private:
  ExpectationBasedPointSetMetric( const Self & ); //purposely not implemented
  void operator=( const Self & );               //purposely not implemented

  CoordRepType                                     m_Sigma;
  unsigned int                                     m_EvaluationKNeighborhood;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkExpectationBasedPointSetMetric.hxx"
#endif

#endif
