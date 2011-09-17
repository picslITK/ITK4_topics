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


#ifndef __itkLevelSetEvolutionBase_h
#define __itkLevelSetEvolutionBase_h

#include <list>

#include "itkImage.h"
#include "itkDiscreteLevelSetImageBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkNumericTraits.h"
#include "itkLevelSetEvolutionStoppingCriterionBase.h"

namespace itk
{
/**
 *  \class LevelSetEvolutionBase
 *  \brief Class for iterating and evolving the dense level-set function
 *
 *  \tparam TEquationContainer Container holding the system of level set of equations
 *   \ingroup ITKLevelSetsv4
 */
template< class TEquationContainer >
class LevelSetEvolutionBase : public Object
{
public:
  typedef LevelSetEvolutionBase      Self;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object                     Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEvolutionBase, Object );

  typedef TEquationContainer                      EquationContainerType;
  typedef typename EquationContainerType::Pointer EquationContainerPointer;
  typedef typename EquationContainerType::TermContainerType
                                                  TermContainerType;
  typedef typename TermContainerType::Pointer     TermContainerPointer;

  typedef typename TermContainerType::TermType TermType;
  typedef typename TermType::Pointer           TermPointer;

  typedef typename TermContainerType::InputImageType InputImageType;
  typedef typename InputImageType::PixelType         InputImagePixelType;
  typedef typename InputImageType::ConstPointer      InputImageConstPointer;
  typedef typename InputImageType::RegionType        InputImageRegionType;
  typedef typename NumericTraits< InputImagePixelType >::RealType
                                                     InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int, InputImageType::ImageDimension );

  typedef typename TermContainerType::LevelSetContainerType LevelSetContainerType;
  typedef typename LevelSetContainerType::LevelSetIdentifierType
                                                            LevelSetIdentifierType;
  typedef typename LevelSetContainerType::Pointer           LevelSetContainerPointer;
  typedef typename LevelSetContainerType::LevelSetContainerConstIteratorType
                                                            LevelSetContainerConstIteratorType;
  typedef typename LevelSetContainerType::LevelSetContainerIteratorType
                                                            LevelSetContainerIteratorType;

  typedef typename LevelSetContainerType::LevelSetType  LevelSetType;
  typedef typename LevelSetType::Pointer                LevelSetPointer;
  typedef typename LevelSetType::ImageType              LevelSetImageType;
  typedef typename LevelSetType::OutputType             LevelSetOutputType;
  typedef typename LevelSetType::OutputRealType         LevelSetOutputRealType;
  typedef typename LevelSetImageType::Pointer           LevelSetImagePointer;
  typedef typename LevelSetType::LevelSetDataType       LevelSetDataType;


  typedef BinaryThresholdImageFilter< LevelSetImageType, LevelSetImageType >  ThresholdFilterType;
  typedef typename ThresholdFilterType::Pointer                               ThresholdFilterPointer;

  typedef SignedMaurerDistanceMapImageFilter< LevelSetImageType, LevelSetImageType >  MaurerType;
  typedef typename MaurerType::Pointer                                                MaurerPointer;

  typedef ImageRegionIteratorWithIndex< LevelSetImageType > LevelSetImageIteratorType;

  typedef ImageRegionConstIteratorWithIndex< LevelSetImageType > LevelSetImageConstIteratorType;

  typedef ImageRegionConstIteratorWithIndex< InputImageType > InputImageConstIteratorType;

  typedef std::list< IdentifierType >                    IdListType;
  typedef typename IdListType::iterator                  IdListIterator;
  typedef Image< IdListType, ImageDimension >            IdListImageType;
  typedef Image< short, ImageDimension >                 CacheImageType;
  typedef typename LevelSetContainerType::DomainMapImageFilterType
                                                         DomainMapImageFilterType;
  typedef typename DomainMapImageFilterType::Pointer     DomainMapImageFilterPointer;
  typedef typename DomainMapImageFilterType::DomainIteratorType
                                                         DomainIteratorType;


  typedef LevelSetEvolutionStoppingCriterionBase< LevelSetContainerType >
                                                  StoppingCriterionType;
  typedef typename StoppingCriterionType::Pointer StoppingCriterionPointer;

  itkSetObjectMacro( LevelSetContainer, LevelSetContainerType );
  itkGetObjectMacro( LevelSetContainer, LevelSetContainerType );

  /** Update the filter by computing the output level function
   * by calling RunOneIteration() once the instantiation of necessary variables
   * is verified */
  void Update();

  /** Set/Get the value of alpha for computing the time-step using CFL conditions */
  itkSetMacro( Alpha, LevelSetOutputRealType );
  itkGetMacro( Alpha, LevelSetOutputRealType );

  /** Set a user-specified value of the time-step */
  void SetTimeStep( const LevelSetOutputRealType& iDt );

  /** Set/Get the equation container for updating all the level sets */
  itkSetObjectMacro( EquationContainer, EquationContainerType );
  itkGetObjectMacro( EquationContainer, EquationContainerType );

  /** Set/Get the Stopping Criterion */
  itkGetObjectMacro( StoppingCriterion, StoppingCriterionType );
  itkSetObjectMacro( StoppingCriterion, StoppingCriterionType );


protected:
  LevelSetEvolutionBase();

  virtual ~LevelSetEvolutionBase();

  StoppingCriterionPointer  m_StoppingCriterion;

  EquationContainerPointer    m_EquationContainer;
  LevelSetContainerPointer    m_LevelSetContainer;
  LevelSetContainerPointer    m_UpdateBuffer;

  LevelSetOutputRealType      m_Alpha;
  LevelSetOutputRealType      m_Dt;
  LevelSetOutputRealType      m_RMSChangeAccumulator;
  bool                        m_UserGloballyDefinedTimeStep;

  /** Initialize the update buffers for all level sets to hold the updates of
   *  equations in each iteration */
  void AllocateUpdateBuffer();

  /** Run the iterative loops of calculating levelset function updates until
   *  the stopping criterion is satisfied */
  void RunOneIteration();

  /** Initialize the iteration by computing parameters in the terms of the level set equation */
  void InitializeIteration();

  /** Computer the update at each pixel and store in the update buffer */
  void ComputeIteration();

  /** Compute the time-step for the next iteration */
  void ComputeTimeStepForNextIteration();

  /** Update the levelset by 1 iteration from the computed updates */
  virtual void UpdateLevelSets();

  /** Update the equations at the end of 1 iteration */
  void UpdateEquations();

  /** Reinitialize the level set functions to a signed distance function */
  void Reinitialize();
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLevelSetEvolutionBase.hxx"
#endif

#endif // __itkLevelSetEvolutionBase_h
