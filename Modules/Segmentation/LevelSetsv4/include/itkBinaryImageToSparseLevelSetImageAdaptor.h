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

#ifndef __itkBinaryImageToSparseLevelSetImageAdaptor_h
#define __itkBinaryImageToSparseLevelSetImageAdaptor_h

#include "itkBinaryImageToLevelSetImageAdaptorBase.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"

#include "itkImage.h"

#include "itkWhitakerSparseLevelSetImage.h"
#include "itkShiSparseLevelSetImage.h"
#include "itkMalcolmSparseLevelSetImage.h"

namespace itk
{
template< class TInput, class TOutput >
class BinaryImageToSparseLevelSetImageAdaptorBase :
    public BinaryImageToLevelSetImageAdaptorBase< TInput, TOutput >
{
public:
  typedef BinaryImageToSparseLevelSetImageAdaptorBase Self;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;
  typedef BinaryImageToLevelSetImageAdaptorBase< TInput, TOutput >
    Superclass;

  /** Run-time type information */
  itkTypeMacro( BinaryImageToSparseLevelSetImageAdaptor,
                BinaryImageToLevelSetImageAdaptorBase );

  typedef typename Superclass::InputImageType       InputImageType;
  typedef typename Superclass::InputImagePixelType  InputImagePixelType;
  typedef typename Superclass::InputImageIndexType  InputImageIndexType;
  typedef typename Superclass::InputImagePointer    InputImagePointer;
  typedef typename Superclass::InputImageRegionType InputImageRegionType;
  typedef typename Superclass::InputPixelRealType   InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

  typedef typename Superclass::LevelSetType             LevelSetType;
  typedef typename Superclass::LevelSetPointer          LevelSetPointer;

  typedef typename LevelSetType::InputType              LevelSetInputType;
  typedef typename LevelSetType::OutputType             LevelSetOutputType;

  typedef typename LevelSetType::LabelObjectType        LevelSetLabelObjectType;
  typedef typename LevelSetLabelObjectType::LabelType   LayerIdType;
  typedef typename LevelSetType::LabelObjectPointer     LevelSetLabelObjectPointer;
  typedef typename LevelSetType::LabelObjectLengthType  LevelSetLabelObjectLengthType;
  typedef typename LevelSetType::LabelObjectLineType    LevelSetLabelObjectLineType;

  typedef typename LevelSetType::LabelMapType           LevelSetLabelMapType;
  typedef typename LevelSetType::LabelMapPointer        LevelSetLabelMapPointer;

  typedef typename LevelSetType::LayerType              LevelSetLayerType;
  typedef typename LevelSetType::LayerIterator          LevelSetLayerIterator;
  typedef typename LevelSetType::LayerConstIterator     LevelSetLayerConstIterator;

  typedef Image< char, ImageDimension >         InternalImageType;
  typedef typename InternalImageType::Pointer   InternalImagePointer;

  typedef std::pair< LevelSetInputType, LevelSetOutputType >  LayerPairType;

  typedef ImageRegionIteratorWithIndex< InputImageType >      InputIteratorType;
  typedef ImageRegionIteratorWithIndex< InternalImageType >   InternalIteratorType;

  typedef ShapedNeighborhoodIterator< InternalImageType > NeighborhoodIteratorType;

protected:
  BinaryImageToSparseLevelSetImageAdaptorBase() : Superclass() {}
  virtual ~BinaryImageToSparseLevelSetImageAdaptorBase() {}

  LevelSetLabelMapPointer m_LabelMap;

  InternalImagePointer m_InternalImage;

private:
  BinaryImageToSparseLevelSetImageAdaptorBase( const Self& );
  void operator = ( const Self& );
};

////////////////////////////////////////////////////////////////////////////////
template< class TInput, class TOutput >
class BinaryImageToSparseLevelSetImageAdaptor
{
};

////////////////////////////////////////////////////////////////////////////////
template< class TInput, typename TOutput >
class BinaryImageToSparseLevelSetImageAdaptor<
    TInput,
    WhitakerSparseLevelSetImage< TOutput, TInput::ImageDimension > > :
  public BinaryImageToSparseLevelSetImageAdaptorBase<
      TInput,
      WhitakerSparseLevelSetImage< TOutput, TInput::ImageDimension > >
  {
public:
  typedef WhitakerSparseLevelSetImage< TOutput, TInput::ImageDimension >
    LevelSetType;

  typedef BinaryImageToSparseLevelSetImageAdaptor Self;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;
  typedef BinaryImageToSparseLevelSetImageAdaptorBase< TInput, LevelSetType >
    Superclass;


  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( BinaryImageToSparseLevelSetImageAdaptor,
                BinaryImageToLevelSetImageAdaptorBase );

  typedef typename Superclass::InputImageType       InputImageType;
  typedef typename Superclass::InputImagePixelType  InputImagePixelType;
  typedef typename Superclass::InputImageIndexType  InputImageIndexType;
  typedef typename Superclass::InputImagePointer    InputImagePointer;
  typedef typename Superclass::InputImageRegionType InputImageRegionType;
  typedef typename Superclass::InputPixelRealType   InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

//  typedef typename Superclass::LevelSetType             LevelSetType;
  typedef typename Superclass::LevelSetPointer                LevelSetPointer;

  typedef typename Superclass::LevelSetInputType              LevelSetInputType;
  typedef typename Superclass::LevelSetOutputType             LevelSetOutputType;

  typedef typename Superclass::LevelSetLabelObjectType        LevelSetLabelObjectType;
  typedef typename Superclass::LayerIdType                    LayerIdType;
  typedef typename Superclass::LevelSetLabelObjectPointer     LevelSetLabelObjectPointer;
  typedef typename Superclass::LevelSetLabelObjectLengthType  LevelSetLabelObjectLengthType;
  typedef typename Superclass::LevelSetLabelObjectLineType    LevelSetLabelObjectLineType;

  typedef typename Superclass::LevelSetLabelMapType           LevelSetLabelMapType;
  typedef typename Superclass::LevelSetLabelMapPointer        LevelSetLabelMapPointer;

  typedef typename Superclass::LevelSetLayerType              LevelSetLayerType;
  typedef typename Superclass::LevelSetLayerIterator          LevelSetLayerIterator;
  typedef typename Superclass::LevelSetLayerConstIterator     LevelSetLayerConstIterator;

  typedef typename Superclass::InternalImageType        InternalImageType;
  typedef typename Superclass::InternalImagePointer     InternalImagePointer;

  typedef typename Superclass::LayerPairType            LayerPairType;

  typedef typename Superclass::InputIteratorType        InputIteratorType;
  typedef typename Superclass::InternalIteratorType     InternalIteratorType;

  typedef typename Superclass::NeighborhoodIteratorType NeighborhoodIteratorType;

  void Initialize();

protected:
  /** Constructor */
  BinaryImageToSparseLevelSetImageAdaptor();

  /** Destructor */
  virtual ~BinaryImageToSparseLevelSetImageAdaptor();

private:

  BinaryImageToSparseLevelSetImageAdaptor( const Self& ); // purposely not implemented
  void operator = ( const Self& );  // purposely not implemented

  /** Fill layer adjacent (OutputLayer) to the layer (LayerToBeScanned) */
  void PropagateToOuterLayers( LayerIdType LayerToBeScanned, LayerIdType OutputLayer, LayerIdType TestValue );

  /** Fill the layer corresponding to zero level set */
  void FindActiveLayer();

  /** Fill layers adjacent to the zero level set (i.e. layer -1 and +1 )*/
  void FindPlusOneMinusOneLayer();

};

////////////////////////////////////////////////////////////////////////////////
template< class TInput >
class BinaryImageToSparseLevelSetImageAdaptor<
    TInput,
    ShiSparseLevelSetImage< TInput::ImageDimension > > :
    public BinaryImageToSparseLevelSetImageAdaptorBase< TInput, ShiSparseLevelSetImage< TInput::ImageDimension > >
{
public:
  typedef ShiSparseLevelSetImage< TInput::ImageDimension > LevelSetType;

  typedef BinaryImageToSparseLevelSetImageAdaptor Self;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;
  typedef BinaryImageToSparseLevelSetImageAdaptorBase< TInput, LevelSetType >
    Superclass;


  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( BinaryImageToSparseLevelSetImageAdaptor,
                BinaryImageToLevelSetImageAdaptorBase );

  typedef typename Superclass::InputImageType       InputImageType;
  typedef typename Superclass::InputImagePixelType  InputImagePixelType;
  typedef typename Superclass::InputImageIndexType  InputImageIndexType;
  typedef typename Superclass::InputImagePointer    InputImagePointer;
  typedef typename Superclass::InputImageRegionType InputImageRegionType;
  typedef typename Superclass::InputPixelRealType   InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

//  typedef typename Superclass::LevelSetType             LevelSetType;
  typedef typename Superclass::LevelSetPointer                LevelSetPointer;

  typedef typename Superclass::LevelSetInputType              LevelSetInputType;
  typedef typename Superclass::LevelSetOutputType             LevelSetOutputType;

  typedef typename Superclass::LevelSetLabelObjectType        LevelSetLabelObjectType;
  typedef typename Superclass::LayerIdType                    LayerIdType;
  typedef typename Superclass::LevelSetLabelObjectPointer     LevelSetLabelObjectPointer;
  typedef typename Superclass::LevelSetLabelObjectLengthType  LevelSetLabelObjectLengthType;
  typedef typename Superclass::LevelSetLabelObjectLineType    LevelSetLabelObjectLineType;

  typedef typename Superclass::LevelSetLabelMapType           LevelSetLabelMapType;
  typedef typename Superclass::LevelSetLabelMapPointer        LevelSetLabelMapPointer;

  typedef typename Superclass::LevelSetLayerType              LevelSetLayerType;
  typedef typename Superclass::LevelSetLayerIterator          LevelSetLayerIterator;
  typedef typename Superclass::LevelSetLayerConstIterator     LevelSetLayerConstIterator;

  typedef typename Superclass::InternalImageType        InternalImageType;
  typedef typename Superclass::InternalImagePointer     InternalImagePointer;

  typedef typename Superclass::LayerPairType            LayerPairType;

  typedef typename Superclass::InputIteratorType        InputIteratorType;
  typedef typename Superclass::InternalIteratorType     InternalIteratorType;

  typedef typename Superclass::NeighborhoodIteratorType NeighborhoodIteratorType;


  void Initialize();

protected:
  /** Constructor */
  BinaryImageToSparseLevelSetImageAdaptor();

  /** Destructor */
  ~BinaryImageToSparseLevelSetImageAdaptor();

  /** Find the active layer separating the foreground and background regions */
  void FindActiveLayer();

private:

  BinaryImageToSparseLevelSetImageAdaptor( const Self& ); // purposely not implemented
  void operator = ( const Self& );  // purposely not implemented
};


////////////////////////////////////////////////////////////////////////////////
template< class TInput >
class BinaryImageToSparseLevelSetImageAdaptor<
    TInput,
    MalcolmSparseLevelSetImage< TInput::ImageDimension > > :
  public BinaryImageToSparseLevelSetImageAdaptorBase< TInput, MalcolmSparseLevelSetImage< TInput::ImageDimension > >
{
public:
  typedef MalcolmSparseLevelSetImage< TInput::ImageDimension > LevelSetType;

  typedef BinaryImageToSparseLevelSetImageAdaptor Self;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;
  typedef BinaryImageToSparseLevelSetImageAdaptorBase< TInput, LevelSetType >
    Superclass;


  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( BinaryImageToSparseLevelSetImageAdaptor,
                BinaryImageToLevelSetImageAdaptorBase );

  typedef typename Superclass::InputImageType       InputImageType;
  typedef typename Superclass::InputImagePixelType  InputImagePixelType;
  typedef typename Superclass::InputImageIndexType  InputImageIndexType;
  typedef typename Superclass::InputImagePointer    InputImagePointer;
  typedef typename Superclass::InputImageRegionType InputImageRegionType;
  typedef typename Superclass::InputPixelRealType   InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

//  typedef typename Superclass::LevelSetType             LevelSetType;
  typedef typename Superclass::LevelSetPointer                LevelSetPointer;

  typedef typename Superclass::LevelSetInputType              LevelSetInputType;
  typedef typename Superclass::LevelSetOutputType             LevelSetOutputType;

  typedef typename Superclass::LevelSetLabelObjectType        LevelSetLabelObjectType;
  typedef typename Superclass::LayerIdType                    LayerIdType;
  typedef typename Superclass::LevelSetLabelObjectPointer     LevelSetLabelObjectPointer;
  typedef typename Superclass::LevelSetLabelObjectLengthType  LevelSetLabelObjectLengthType;
  typedef typename Superclass::LevelSetLabelObjectLineType    LevelSetLabelObjectLineType;

  typedef typename Superclass::LevelSetLabelMapType           LevelSetLabelMapType;
  typedef typename Superclass::LevelSetLabelMapPointer        LevelSetLabelMapPointer;

  typedef typename Superclass::LevelSetLayerType              LevelSetLayerType;
  typedef typename Superclass::LevelSetLayerIterator          LevelSetLayerIterator;
  typedef typename Superclass::LevelSetLayerConstIterator     LevelSetLayerConstIterator;

  typedef typename Superclass::InternalImageType        InternalImageType;
  typedef typename Superclass::InternalImagePointer     InternalImagePointer;

  typedef typename Superclass::LayerPairType            LayerPairType;

  typedef typename Superclass::InputIteratorType        InputIteratorType;
  typedef typename Superclass::InternalIteratorType     InternalIteratorType;

  typedef typename Superclass::NeighborhoodIteratorType NeighborhoodIteratorType;

  void Initialize();

protected:
  /** Constructor */
  BinaryImageToSparseLevelSetImageAdaptor();

  /** Destructor */
  virtual ~BinaryImageToSparseLevelSetImageAdaptor();

  /** Find the active layer separating the foreground and background regions */
  void FindActiveLayer();

  /** Ensure that the 0 level set layer is only of single pixel thickness */
  void CreateMinimalInterface();

private:

  BinaryImageToSparseLevelSetImageAdaptor( const Self& ); // purposely not implemented
  void operator = ( const Self& );  // purposely not implemented
};

}
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryImageToSparseLevelSetImageAdaptor.hxx"
#endif
#endif
