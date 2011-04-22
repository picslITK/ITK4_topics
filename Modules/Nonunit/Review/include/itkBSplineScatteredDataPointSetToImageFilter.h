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
#ifndef __itkBSplineScatteredDataPointSetToImageFilter_h
#define __itkBSplineScatteredDataPointSetToImageFilter_h

#include "itkPointSetToImageFilter.h"
#include "itkBSplineKernelFunction.h"
#include "itkCoxDeBoorBSplineKernelFunction.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorContainer.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

namespace itk
{
/** \class BSplineScatteredDataPointSetToImageFilter
 * \brief Image filter which provides a B-spline output approximation.
 *
 * Given an n-D image with scattered data, this filter finds
 * a fast approximation to that irregularly spaced data using uniform
 * B-splines.  The traditional method of inverting the observation
 * matrix to find a least-squares fit is made obsolete.  Therefore,
 * memory issues are not a concern and inverting large matrices is
 * not applicable.  In addition, this allows fitting to be multi-threaded.
 * This class generalizes from Lee's original paper to encompass
 * n-D data in m-D parametric space and any *feasible* B-spline order as well
 * as the option of specifying a confidence value for each point.
 *
 * In addition to specifying the input point set, one must specify the number
 * of control points.  The specified number of control points must be
 * > m_SplineOrder.  If one wishes to use the multilevel component of
 * this algorithm, one must also specify the number of levels in the
 * hierarchy.  If this is desired, the number of control points becomes
 * the number of control points for the coarsest level.  The algorithm
 * then increases the number of control points at each level so that
 * the B-spline n-D grid is refined to twice the previous level.
 *
 * \author Nicholas J. Tustison
 *
 * Contributed by Nicholas J. Tustison, James C. Gee
 * in the Insight Journal paper:
 * http://hdl.handle.net/1926/140
 *
 * \par REFERENCE
 * S. Lee, G. Wolberg, and S. Y. Shin, "Scattered Data Interpolation
 * with Multilevel B-Splines", IEEE Transactions on Visualization and
 * Computer Graphics, 3(3):228-244, 1997.
 *
 * \par REFERENCE
 * N.J. Tustison and J.C. Gee, "Generalized n-D C^k Scattered Data Approximation
 * with COnfidence Values", Proceedings of the MIAR conference, August 2006.
 * \ingroup ITK-Review
 */

template< class TInputPointSet, class TOutputImage >
class BSplineScatteredDataPointSetToImageFilter:
  public PointSetToImageFilter< TInputPointSet, TOutputImage >
{
public:
  typedef BSplineScatteredDataPointSetToImageFilter             Self;
  typedef PointSetToImageFilter<TInputPointSet, TOutputImage>   Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Extract dimension from the output image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TOutputImage::ImageDimension );

  typedef TOutputImage                              ImageType;
  typedef TInputPointSet                            PointSetType;

  /** Image typedef support. */
  typedef typename ImageType::PixelType             PixelType;
  typedef typename ImageType::RegionType            RegionType;
  typedef typename ImageType::SizeType              SizeType;
  typedef typename ImageType::IndexType             IndexType;

  /** PointSet typedef support. */
  typedef typename PointSetType::PointType          PointType;
  typedef typename PointSetType::Pointer            PointSetPointer;
  typedef typename PointSetType::PixelType          PointDataType;
  typedef typename PointSetType::PointDataContainer PointDataContainerType;

  /** Other typedef */
  typedef float                                     RealType;
  typedef VectorContainer<unsigned, RealType>       WeightsContainerType;

  /** Image types */
  typedef Image<PointDataType,
    itkGetStaticConstMacro( ImageDimension )>       PointDataImageType;
  typedef Image<RealType,
    itkGetStaticConstMacro( ImageDimension )>       RealImageType;
  typedef typename RealImageType::Pointer           RealImagePointer;
  typedef typename PointDataImageType::Pointer      PointDataImagePointer;
  typedef FixedArray<unsigned,
    itkGetStaticConstMacro(ImageDimension) >        ArrayType;

  /**
   * Interpolation kernel type (default spline order = 3)
   */
  typedef CoxDeBoorBSplineKernelFunction<3>         KernelType;
  typedef BSplineKernelFunction<0>                  KernelOrder0Type;
  typedef BSplineKernelFunction<1>                  KernelOrder1Type;
  typedef BSplineKernelFunction<2>                  KernelOrder2Type;
  typedef BSplineKernelFunction<3>                  KernelOrder3Type;

  // Helper functions

  /**
   * Set the spline order assuming it is the same in all parametric dimensions.
   * The spline order determines the continuity between B-spline elements and
   * the degree of polynomial used to construct the B-spline elements.  Default
   * = 3.
   */
  void SetSplineOrder( unsigned int );

  /**
   * Set the spline order for each parametric dimension separately.  The spline
   * order determines the continuity between B-spline elements and the degree of
   * polynomial used to construct the B-spline elements.  Default = 3.
   */
  void SetSplineOrder( const ArrayType & );

  /**
   * Get the spline order for all parametric dimensions.  The spline order
   * determines the continuity between B-spline elements and the degree of
   * polynomial used to construct the B-spline elements.  Default = 3.
   */
  itkGetConstReferenceMacro( SplineOrder, ArrayType );

  /**
   * Set the number of control points for each parametric dimension at the
   * initial fitting level.  The B-spline mesh size is equal to the number
   * of control points minus the spline order.  Default = 4 in each dimension.
   */
  itkSetMacro( NumberOfControlPoints, ArrayType );

  /**
   * Set the number of control points for each parametric dimension at the
   * initial fitting level.  The B-spline mesh size is equal to the number
   * of control points minus the spline order.  Default = 4 in each dimension.
   */
  itkGetConstReferenceMacro( NumberOfControlPoints, ArrayType );

  /**
   * Get the number of current control points for each parametric dimension at
   * the current fitting level.  The B-spline mesh size is equal to the number
   * of control points minus the spline order.  Default = 4 in each dimension.
   */
  itkGetConstReferenceMacro( CurrentNumberOfControlPoints, ArrayType );

  /**
   * Set the number of fitting levels assuming the number of fitting levels is
   * the same for each parametric dimension.  Starting with the mesh size
   * implied by setting the number of control points, the mesh size is doubled
   * at each fitting level.  Default = 1 in all parametric dimensions.
   */
  void SetNumberOfLevels( unsigned int );

  /**
   * Set the number of fitting levels in each parametric dimension separately.
   * Starting with the mesh size implied by setting the number of control
   * points, the mesh size is doubled at each fitting level.  Default = 1 in all
   * parametric dimensions.
   */
  void SetNumberOfLevels( const ArrayType & );

  /**
   * Get the number of fitting levels for all parametric dimensions. Starting
   * with the mesh size implied by setting the number of control points, the
   * mesh size is doubled at each fitting level.  Default = 1 in all parametric
   * dimensions.
   */
  itkGetConstReferenceMacro( NumberOfLevels, ArrayType );

  /**
   * This array of 0/1 values defines whether a particular dimension of the
   * parametric space is to be considered periodic or not. For example, if you
   * are using interpolating along a 1D closed curve, the array type will have
   * size 1, and you should set the first element of this array to the value
   * "1". In the case that you were interpolating in a planar surface with
   * cylindrical topology, the array type will have two components, and you
   * should set to "1" the component that goes around the cylinder, and set to
   * "0" the component that goes from the top of the cylinder to the bottom.
   * This will indicate the periodity of that parameter to the filter.
   * Internally, in order to make periodic the domain of the parameter, the
   * filter will reuse some of the points at the beginning of the domain as if
   * they were also located at the end of the domain. The number of points to
   * be reused will depend on the spline order. As a user, you don't need to
   * replicate the points, the filter will do this for you.
   */
  itkSetMacro( CloseDimension, ArrayType );

  /**
   * This array of 0/1 values defines whether a particular dimension of the
   * parametric space is to be considered periodic or not. For example, if you
   * are using interpolating along a 1D closed curve, the array type will have
   * size 1, and you should set the first element of this array to the value
   * "1". In the case that you were interpolating in a planar surface with
   * cylindrical topology, the array type will have two components, and you
   * should set to "1" the component that goes around the cylinder, and set to
   * "0" the component that goes from the top of the cylinder to the bottom.
   * This will indicate the periodity of that parameter to the filter.
   * Internally, in order to make periodic the domain of the parameter, the
   * filter will reuse some of the points at the beginning of the domain as if
   * they were also located at the end of the domain. The number of points to
   * be reused will depend on the spline order. As a user, you don't need to
   * replicate the points, the filter will do this for you.
   */
  itkGetConstReferenceMacro( CloseDimension, ArrayType );

  /**
   * A weighted fitting is possible where each input point is assigned a
   * relative weighting.
   */
  void SetPointWeights( WeightsContainerType *weights );

  /**
   * The result of the fitting process is an n-D grid of control points which
   * describe the continuous B-spline object.  This boolean value determines
   * whether or not this sampled B-spline object is constructed.
   */
  itkSetMacro( GenerateOutputImage, bool );

  /**
   * The result of the fitting process is an n-D grid of control points which
   * describe the continuous B-spline object.  This boolean value determines
   * whether or not this sampled B-spline object is constructed.
   */
  itkGetConstReferenceMacro( GenerateOutputImage, bool );

  /**
   * The result of the fitting process is an n-D grid of control points which
   * describe the continuous B-spline object.  This boolean value determines
   * whether or not this sampled B-spline object is constructed.
   */
  itkBooleanMacro( GenerateOutputImage );

  /**
   * Set the control point lattice produced by a previous fitting process.
   */
  itkSetMacro( PhiLattice, PointDataImagePointer );

  /**
   * Get the control point lattice produced by the fitting process.
   */
  itkGetConstMacro( PhiLattice, PointDataImagePointer );

protected:
  BSplineScatteredDataPointSetToImageFilter();
  virtual ~BSplineScatteredDataPointSetToImageFilter();

  void PrintSelf(std::ostream & os, Indent indent) const;

  void ThreadedGenerateData( const RegionType &, int );

  void BeforeThreadedGenerateData();

  void AfterThreadedGenerateData();

  int SplitRequestedRegion( int, int, RegionType & );

  void GenerateData();

private:

  //purposely not implemented
  BSplineScatteredDataPointSetToImageFilter( const Self & );
  void operator=( const Self & );

  /**
   * Function used to propagate the fitting solution at one fitting level
   * to the next level with the mesh resolution doubled.
   */
  void RefineControlPointLattice();

  /**
   * Determine the residuals after fitting to one level.
   */
  void UpdatePointSet();

  /**
   * This function is not used as it requires an evaluation of all
   * (SplineOrder+1)^ImageDimensions B-spline weights for each evaluation.
   */
  void GenerateOutputImage();

  /**
   * Function used to generate the sampled B-spline object quickly.
   */
  void ThreadedGenerateDataForFitting( const RegionType &, int  );

  /**
   * Function used to generate the sampled B-spline object quickly.
   */
  void ThreadedGenerateDataForReconstruction( const RegionType &, int  );

  /**
   * Sub-function used by GenerateOutputImageFast() to generate the sampled
   * B-spline object quickly.
   */
  void CollapsePhiLattice( PointDataImageType *, PointDataImageType *,
    const RealType, const unsigned int );

  /**
   * Set the grid parametric domain parameters such as the origin, size,
   * spacing, and direction.
   */
  void SetPhiLatticeParametricDomainParameters();

  bool                                         m_DoMultilevel;
  bool                                         m_GenerateOutputImage;
  bool                                         m_UsePointWeights;
  unsigned int                                 m_MaximumNumberOfLevels;
  unsigned int                                 m_CurrentLevel;
  ArrayType                                    m_NumberOfControlPoints;
  ArrayType                                    m_CurrentNumberOfControlPoints;
  ArrayType                                    m_CloseDimension;
  ArrayType                                    m_SplineOrder;
  ArrayType                                    m_NumberOfLevels;

  typename WeightsContainerType::Pointer       m_PointWeights;

  typename PointDataImageType::Pointer         m_PhiLattice;
  typename PointDataImageType::Pointer         m_PsiLattice;

  vnl_matrix<RealType>     m_RefinedLatticeCoefficients[ImageDimension];

  typename PointDataContainerType::Pointer     m_InputPointData;
  typename PointDataContainerType::Pointer     m_OutputPointData;

  typename KernelType::Pointer                 m_Kernel[ImageDimension];

  typename KernelOrder0Type::Pointer           m_KernelOrder0;
  typename KernelOrder1Type::Pointer           m_KernelOrder1;
  typename KernelOrder2Type::Pointer           m_KernelOrder2;
  typename KernelOrder3Type::Pointer           m_KernelOrder3;

  std::vector<RealImagePointer>                m_OmegaLatticePerThread;
  std::vector<PointDataImagePointer>           m_DeltaLatticePerThread;

  RealType                                     m_BSplineEpsilon;
  bool                                         m_IsFittingComplete;

  inline typename RealImageType::IndexType
  NumberToIndex( unsigned int number, typename RealImageType::SizeType size )
    {
    typename RealImageType::IndexType k;
    k[0] = 1;

    for ( unsigned int i = 1; i < ImageDimension; i++ )
      {
      k[i] = size[ImageDimension - i - 1] * k[i - 1];
      }
    typename RealImageType::IndexType index;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      index[ImageDimension - i - 1] =
        static_cast< unsigned int >( number / k[ImageDimension - i - 1] );
      number %= k[ImageDimension - i - 1];
      }
    return index;
    }
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBSplineScatteredDataPointSetToImageFilter.txx"
#endif

#endif
