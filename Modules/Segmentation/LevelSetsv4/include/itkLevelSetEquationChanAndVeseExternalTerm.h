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

#ifndef __itkLevelSetEquationChanAndVeseExternalTerm_h
#define __itkLevelSetEquationChanAndVeseExternalTerm_h

#include "itkLevelSetEquationChanAndVeseInternalTerm.h"
#include "itkNumericTraits.h"

namespace itk
{
/**
 *  \class LevelSetEquationChanAndVeseExternalTerm
 *  \brief Derived class to represent the external energy Chan And Vese term
 *
 *  \tparam TInput Input Image Type
 *  \tparam TLevelSetContainer Level set function container type
 *  \ingroup ITKLevelSetsv4
 */
template< class TInput, // Input image or mesh
          class TLevelSetContainer >
class LevelSetEquationChanAndVeseExternalTerm :
    public LevelSetEquationChanAndVeseInternalTerm< TInput, TLevelSetContainer >
{
public:
  typedef LevelSetEquationChanAndVeseExternalTerm         Self;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  typedef LevelSetEquationChanAndVeseInternalTerm< TInput,
                                    TLevelSetContainer >  Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationChanAndVeseExternalTerm,
                LevelSetEquationChanAndVeseInternalTerm );

  typedef typename Superclass::InputImageType     InputImageType;
  typedef typename Superclass::InputImagePointer  InputImagePointer;
  typedef typename Superclass::InputPixelType     InputPixelType;
  typedef typename Superclass::InputPixelRealType InputPixelRealType;

  typedef typename Superclass::LevelSetContainerType      LevelSetContainerType;
  typedef typename Superclass::LevelSetContainerPointer   LevelSetContainerPointer;
  typedef typename Superclass::LevelSetType               LevelSetType;
  typedef typename Superclass::LevelSetPointer            LevelSetPointer;
  typedef typename Superclass::LevelSetOutputPixelType    LevelSetOutputPixelType;
  typedef typename Superclass::LevelSetOutputRealType     LevelSetOutputRealType;
  typedef typename Superclass::LevelSetInputIndexType     LevelSetInputIndexType;
  typedef typename Superclass::LevelSetGradientType       LevelSetGradientType;
  typedef typename Superclass::LevelSetHessianType        LevelSetHessianType;
  typedef typename Superclass::LevelSetIdentifierType     LevelSetIdentifierType;

  typedef std::list< LevelSetIdentifierType >            IdListType;
  typedef typename IdListType::iterator                  IdListIterator;

  typedef typename Superclass::HeavisideType    HeavisideType;
  typedef typename Superclass::HeavisidePointer HeavisidePointer;

  /** Compute the product of Heaviside functions in the multi-levelset cases */
  virtual void ComputeProduct( const LevelSetInputIndexType& iP,
                               LevelSetOutputRealType& prod );

  /** Compute the product of Heaviside functions in the multi-levelset cases
   *  except the current levelset */
  virtual void ComputeProductTerm( const LevelSetInputIndexType& iP,
                                  LevelSetOutputRealType& prod );

  /** Supply updates at pixels to keep the term parameters always updated */
  virtual void UpdatePixel( const LevelSetInputIndexType& iP,
                           const LevelSetOutputRealType & oldValue,
                           const LevelSetOutputRealType & newValue );

protected:
  LevelSetEquationChanAndVeseExternalTerm();
  virtual ~LevelSetEquationChanAndVeseExternalTerm();

private:
  LevelSetEquationChanAndVeseExternalTerm( const Self& ); // purposely not implemented
  void operator = ( const Self& ); // purposely not implemented
};

}
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLevelSetEquationChanAndVeseExternalTerm.hxx"
#endif

#endif
