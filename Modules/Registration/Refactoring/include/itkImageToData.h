/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageSource.h,v $
  Language:  C++
  Date:      $Date: 2009-03-12 01:11:08 $
  Version:   $Revision: 1.59 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToData_h
#define __itkImageToData_h

#include "itkObjectToDataBase.h"
#include "itkImage.h"

namespace itk
{

/** \class ImageToData
 *  \brief Class for custom threading setup and dispatch of Image data.
 *
 * This class provides threading over an ImageRegion.
 * The SplitRequestedObject method splits the provided ImageRegion into
 * subregions, along the z-axis (?) ... details on motivation?
 *
 * Call SetOverallRegion to define the ImageRegion over which to thread.
 * Call SetThreadedGenerateData to define the worker callback function,
 *  which is called from each thread with a unique region to process.
 *\warning This callback function must be \c static if it is a class method.
 * Call SetHolder to provide a class instance...
 *
 * This class is templated over image dimension. The second template
 * parameter should always be left as default.
 *
 * \sa ObjectToDataBase
 * \ingroup DataSources
 */


template <unsigned int VDimension, class TDataHolder,
          typename TInputObject = ImageRegion<VDimension> >
class ITK_EXPORT ImageToData
  : public ObjectToDataBase<TInputObject, TDataHolder>
{
public:
  /** Standard class typedefs. */
  typedef ImageToData                                 Self;
  typedef ObjectToDataBase<TInputObject, TDataHolder> Superclass;
  typedef SmartPointer<Self>                          Pointer;
  typedef SmartPointer<const Self>                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToData,ObjectToDataBase);

  /** Type of the object being threaded over */
  typedef typename Superclass::InputObjectType  InputObjectType;

  /** Some convenient typedefs. */
  // typedef TImageRegion ImageRegionType;
  itkStaticConstMacro(ImageDimension, unsigned int, VDimension);

  typedef ImageRegion<VDimension>   ImageRegionType;
  typedef Size<VDimension>          SizeType;
  typedef Index<VDimension>         IndexType;

  /** Set the overall image region over which to operate.
   * This is equivalent to SetOverallObject, but named more intuitively
   * for this derived class. */
  virtual void SetOverallRegion(  ImageRegionType& region )
  { this->SetOverallObject( region ); }

protected:
  ImageToData();
  virtual ~ImageToData() {}

  /** Split the ImageRegion \c overallRegion into \c requestedTotal subregions,
   * returning subregion \c i as \c splitRegion.
   * This method is called \c requestedTotal times. The
   * pieces must not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedObject,
   * i.e. return value is less than or equal to \c requestedTotal. */
  virtual
  ThreadIdType SplitRequestedObject(ThreadIdType i,
                           ThreadIdType requestedTotal,
                           InputObjectType& overallRegion,
                           InputObjectType& splitRegion) const;

private:

  ImageToData(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkImageToData.txx"
#endif

#endif
