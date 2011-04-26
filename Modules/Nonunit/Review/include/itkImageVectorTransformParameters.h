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
#ifndef __itkImageVectorTransformParameters_h
#define __itkImageVectorTransformParameters_h

#include "itkTransformParameters.h"
#include "itkImage.h"

namespace itk
{
/** \class ImageVectorTransformParameters
 *  \brief Class to hold and manage parameters of Image<Vector> type used
 *  in Transforms.
 *
 */

/* Can we template of Image type instead, but require that Image be of type
 * Image< Vector< TValueType, NVectorDimension >, VImageDimension > ? */
template< typename TValueType,
          unsigned int NVectorDimension,
          unsigned int VImageDimension >
class ImageVectorTransformParameters : public TransformParameters< TValueType >
{
public:

  /** The element type stored at each location in the Array. */
  typedef TValueType                          ValueType;
  typedef ImageVectorTransformParameters      Self;
  typedef TransformParameters< TValueType >   Superclass;

  /** Image type that this class expects. */
  typedef Image< Vector<TValueType, NVectorDimension>,
                 VImageDimension >
                                                ParameterImageType;
  typedef typename ParameterImageType::Pointer  ParameterImagePointer;

  /** Default constructor. It is created with an empty array
   *  it has to be allocated later by assignment              */
  ImageVectorTransformParameters();

  /** Copy constructor.  Uses VNL copy construtor with correct
   *  setting for memory management.
   *  The vnl vector copy constructor creates new memory
   *  no matter the setting of let array manage memory of rhs.
   *
   * TODO Determine behavior when copying from obj pointing to image parameters.
   *  By default should copy image param data into Array portion of new object,
   *  i.e. into data_block. Is that what we want?
   */
  ImageVectorTransformParameters(const ImageVectorTransformParameters& rhs);

  /** Constructor with size. Size can only be changed by assignment */
  explicit ImageVectorTransformParameters(unsigned int dimension);

  /** Set a new data pointer for *both* the Array and parameter image,
   * pointing both to a different memory block.
   * The size of the new memroy block must be the same as current size of
   * Array and the parameter image's buffer, in elements of TValueType.
   * Memory must be managed by caller afterwards. */
  virtual void MoveDataPointer( TValueType * pointer );

  /** Set an image that holds the parameter data. The Array will be pointed
   * to the image data buffer, and set not to manage memory, so the image
   * still manages its memory. */
  void SetParameterImage( ParameterImageType * );

  /** Copy opertor */
  const Self & operator=(const Self & rhs);

  const Self & operator=(const Superclass & rhs);

  ~ImageVectorTransformParameters(){}

private:
  /** The parameter image used by the class */
  ParameterImagePointer           m_ParameterImage;

};

}//namespace itk

#if ITK_TEMPLATE_TXX
#include "itkImageVectorTransformParameters.txx"
#endif

#endif
