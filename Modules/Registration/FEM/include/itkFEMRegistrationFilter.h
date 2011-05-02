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
#ifndef __itkFEMRegistrationFilter_h
#define __itkFEMRegistrationFilter_h

#include "itkFEMLinearSystemWrapperItpack.h"
#include "itkFEMLinearSystemWrapperDenseVNL.h"
#include "itkFEMGenerateMesh.h"
#include "itkFEMSolverCrankNicolson.h"
#include "itkFEMMaterialLinearElasticity.h"
#include "itkFEMImageMetricLoad.h"
#include "itkFEMFiniteDifferenceFunctionLoad.h"

#include "itkVector.h"
#include "itkVectorCastImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkImageToImageMetric.h"
#include "itkVectorExpandImageFilter.h"

#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkFEMLoadLandmark.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector_fixed.h"

#include <iostream>
#include <string>

namespace itk
{
namespace fem
{
/** \class FEMRegistrationFilter
 *    \brief FEM Image registration filter.
 * The image registration problem is modeled here with the finite
 * element method. Image registration is, in general, an ill-posed
 * problem.  Thus, we use an optimization scheme where the
 * optimization criterion is given by a regularized variational
 * energy. The variational energy arises from modeling the image as a
 * physical body on which external forces act.  The body is allowed to
 * deform so as to minimize the applied force.  The resistance of the
 * physical body to deformation, determined by the physics associated
 * with the body, serves to regularize the solution. The forces
 * applied to the body are, generally, highly non-linear and so the
 * body is allowed to deform slowly and incrementally.  The direction
 * it deforms follows the gradient of the potential energy (the force)
 * we define.  The potential energies we may choose from are given by
 * the itk image-to-image metrics. The choices and the associated
 * direction of descent are :
 *       Mean Squares (minimize),
 *       Normalized Cross-Correlation (maximize)
 *       Mutual Information (maximize).
 *    Note that we have to set the direction (SetDescentDirection)
 *    when we choose a metric.
 * The forces driving the problem may also be given by user-supplied
 * landmarks. The corners of the image, in this example, are always
 * pinned.  This example is designed for 2D or 3D images.  A
 * rectilinear mesh is generated automatically given the correct
 * element type (Quadrilateral or Hexahedral). Our specific Solver for
 * this example uses trapezoidal time stepping.  This is a method for
 * solving a second-order PDE in time.  The solution is penalized by
 * the zeroth (mass matrix) and first derivatives (stiffness matrix)
 * of the shape functions.  There is an option to perform a line
 * search on the energy after each iteration.  Optimal parameter
 * settings require experimentation.
 *   The following approach tends to work well :
 *       Choose the relative size of density  to elasticity (e.g. Rho
 *       / E ~= 1.) such that the image deforms locally and
 *       slowly. This also affects the stability of the
 *       solution. Choose the time step to control the size of the
 *       deformation at each step. Choose enough iterations to allow
 *       the solution to converge (this may be automated).
 *
 *    Reading images is up to the user.  Either set the images using
 *    SetMoving/FixedImage or see the ReadImages function.
 *
 *   \note This code works for only 2 or 3 dimensions b/c we do not
 *   have > 3D elements.
 *
 *   \note TODO :  Keep the full field around (if using
 *   re-gridding). Introduce compensation for kinematic non-linearity
 *   in time (if using Eulerian frame).
 * \ingroup ITK-FEMRegistration
 */

template< class TMovingImage, class TFixedImage >
class ITK_EXPORT FEMRegistrationFilter:public ImageToImageFilter< TMovingImage, TFixedImage >
{
public:
  typedef FEMRegistrationFilter                           Self;
  typedef ImageToImageFilter< TMovingImage, TFixedImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(FEMRegistrationFilter, ImageToImageFilter);

  typedef TMovingImage                       MovingImageType;
  typedef TFixedImage                        FixedImageType;
  typedef typename FixedImageType::PixelType PixelType;
  typedef typename FixedImageType::SizeType  ImageSizeType;

  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int, FixedImageType::ImageDimension);

  typedef Image< float, itkGetStaticConstMacro(ImageDimension) > FloatImageType;
  typedef LinearSystemWrapperItpack                              LinearSystemSolverType;
  typedef SolverCrankNicolson                                    SolverType;
  enum Sign { positive = 1, negative = -1 };
  typedef double          Float;
  typedef Load::ArrayType LoadArray;

  typedef std::vector< typename LoadLandmark::Pointer >                         LandmarkArrayType;
  typedef itk::Vector< float, itkGetStaticConstMacro(ImageDimension) >          VectorType;
  typedef itk::Image< VectorType, itkGetStaticConstMacro(ImageDimension) >      FieldType;
  typedef itk::WarpImageFilter< MovingImageType, FixedImageType, FieldType >    WarperType;
  typedef MaterialLinearElasticity                                              MaterialType;
  typedef itk::ImageRegionIteratorWithIndex< FixedImageType >                   ImageIterator;
  typedef itk::ImageRegionIteratorWithIndex< FloatImageType >                   FloatImageIterator;
  typedef itk::ImageRegionIteratorWithIndex< FieldType >                        FieldIterator;
  typedef itk::VectorIndexSelectionCastImageFilter< FieldType, FloatImageType > IndexSelectCasterType;

  /** Typedef support for the interpolation function */
  typedef double                                                          CoordRepType;
  typedef VectorInterpolateImageFunction< FieldType, CoordRepType >       InterpolatorType;
  typedef typename InterpolatorType::Pointer                              InterpolatorPointer;
  typedef VectorLinearInterpolateImageFunction< FieldType, CoordRepType > DefaultInterpolatorType;

  /** Set the interpolator function. */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the interpolator function. */
  itkGetObjectMacro(Interpolator, InterpolatorType);

  typedef itk::VectorExpandImageFilter< FieldType, FieldType > ExpanderType;
  typedef typename ExpanderType::ExpandFactorsType             ExpandFactorsType;

  typedef itk::RecursiveMultiResolutionPyramidImageFilter< FixedImageType, FixedImageType > FixedPyramidType;

/** Instantiate the load class with the correct image type. */
//#define USEIMAGEMETRIC
#ifdef  USEIMAGEMETRIC
  typedef ImageToImageMetric< ImageType, FixedImageType > MetricBaseType;
  typedef  ImageMetricLoad< ImageType, ImageType >        ImageMetricLoadType;
#else
  typedef  FiniteDifferenceFunctionLoad< MovingImageType, FixedImageType >
  ImageMetricLoadType;
  typedef PDEDeformableRegistrationFunction< FixedImageType, MovingImageType, FieldType >
  MetricBaseType;
#endif
  typedef typename MetricBaseType::Pointer MetricBaseTypePointer;
  /* Main functions */

  /** Read the configuration file to set up the example parameters */
  bool      ReadConfigFile(const char *);

  /** Call this to register two images. */
  void      RunRegistration(void);

  /** Call this to write out images - a counter is attached to the
   *  file name so we can output a numbered sequence tracking the deformation.
   */
  void      WriteWarpedImage(const char *fn);

  /** The solution loop */
  void      IterativeSolve(SolverType & S);

  /** The solution loop for a simple multi-resolution strategy. */
  void      MultiResSolve();

  /** Applies the warp to the input image. */
  void      WarpImage(const MovingImageType *R);

  /** Writes the displacement field to a file. */
  int       WriteDisplacementField(unsigned int index);

  /** Writes the displacement field to a file as a single volume with multiple
    components. */
  int       WriteDisplacementFieldMultiComponent();

  /** One can set the reference file names to read images from files.
   \deprecated  This method currently doesn't have any effect. */
  void      SetMovingFile(const char *r)
  {
    m_MovingFileName = r;
  }

  std::string GetMovingFile()
  {
    return m_MovingFileName;
  }

  /** \deprecated  This method doesn't have any effect */
  void      SetFixedFile(const char *t) { m_FixedFileName = t; }

  std::string GetFixedFile() { return m_FixedFileName; }

  /** One can set the images directly to input images in an application */

  /** Define the reference (moving) image. */
  void SetMovingImage(MovingImageType *R);

  /** Define the target (fixed) image. */
  void SetFixedImage(FixedImageType *T);

  MovingImageType * GetMovingImage(){ return m_MovingImage; }
  MovingImageType * GetOriginalMovingImage(){ return m_OriginalMovingImage; }

  FixedImageType * GetFixedImage(){ return m_FixedImage; }

  /** Get the reference image warped to the target image.
      Must first apply the warp using WarpImage() */
  FixedImageType * GetWarpedImage(){ return m_WarpedImage; }

  /** Compute the jacobian of the current deformation field. */
  void ComputeJacobian(float sign = 1.0, FieldType *field = NULL, float smooth = 0.0);

  /** Get the image that gives the jacobian of the deformation field. */
  FloatImageType * GetJacobianImage(){ return m_FloatImage; }

  /** Outputs the FE deformation field interpolated over the entire
   * image domain. */
  FieldType * GetDeformationField(){ return m_Field; }
  /** Sets the FE deformation field. */
  void SetDeformationField(FieldType *F)
  {
    m_FieldSize = F->GetLargestPossibleRegion().GetSize();
    m_Field = F;
  }

  /** These functions control the use of landmark constraints.  Currently,
      landmarks must be read in from a file. */
  void      SetLandmarkFile(const char *l)
  {
    m_LandmarkFileName = l;
  }

  /** This determines if the landmark file will be read */
  void      UseLandmarks(bool b)
  {
    m_UseLandmarks = b;
  }

  /** We check the jacobian of the current deformation field.
      If it is < threshold, we begin diffeomorphism enforcement:
        1)  Warp the moving image.
        2)  Set the vector field to zero.
        3)  Set the warped moving image as the new moving image,
            resizing if necessary.
    */
  void      EnforceDiffeomorphism(float thresh, SolverType & S,  bool onlywriteimages);

  /** The warped reference image will be written to this file name with
      the extension "11.img" appended to it.  One can also output the
      image after every iteration, yielding result11.img, result12.img, etc.
      by uncommenting the code at the end of IterativeSolve. */
  void      SetResultsFile(const char *r)
  {
    m_ResultsFileName = r;
  }

  void      SetResultsFileName(const char *f)
  {
    m_ResultsFileName = f;
  }

  std::string GetResultsFileName()
  {
    return m_ResultsFileName;
  }

  /** Sets the filename for the vector field component images. */
  void      SetDisplacementsFile(const char *r)
  {
    m_DisplacementsFileName = r;
  }

  /** The FEM filter can generate its own mesh for 2 or 3 dimensions, if none is provided.
      The mesh is generated for quadrilaterals in 2D and hexahedra in 3D.  This function
      sets the number of elements generated along each dimension at the resolution
      designated by "which".
      E.g. to generate 10 pixels per element in each dimension in the 1st resolution, use SetMeshResolution(10,0);.
    */
  void      SetMeshPixelsPerElementAtEachResolution(unsigned int i, unsigned int which = 0)
  {
    m_MeshPixelsPerElementAtEachResolution[which] = i;
  }

  /** This determines the number of integration points to use at each resolution.
      These integration points are used to generate the force.  The actual number
      used will be i^d, where d is the number of parameters in the elements local domain. */
  void      SetNumberOfIntegrationPoints(unsigned int i, unsigned int which = 0)
  {
    m_NumberOfIntegrationPoints[which] = i;
  }

  /** The metric region allows one to compute the derivative (force) of the similarity metric
    * using a region of size [i,i] in 2D [i,i,i] in 3D.
    * \param i number of elements
    * \param which determines the region at a given resolution of the solution process.
    */
  void      SetWidthOfMetricRegion(unsigned int i, unsigned int which = 0)
  {
    m_MetricWidth[which] = i;
  }

  unsigned int      GetWidthOfMetricRegion(unsigned int which = 0)
  {
    return m_MetricWidth[which];
  }

  /** Setting the maximum iterations stops the solution after i iterations regardless of energy.
    * \param i number of elements
    * \param which determines the resolution of the solution process the call is applied to.
    */
  void      SetMaximumIterations(unsigned int i, unsigned int which)
  {
    m_Maxiters[which] = i;
  }

  /** Setting the time step - usually 1.0.  We prefer to use rho to control step sizes.
    */
  void      SetTimeStep(Float i)
  {
    m_Dt = i;
  }

  /** Set alpha for the trapezoidal rule (usually 1.0 in our experiments). */
  void      SetAlpha(Float a)
  {
    m_Alpha = a;
  }

  /** Sets the energy below which we decide the solution has converged.
    */
  void      SetEnergyReductionFactor(Float i)
  {
    m_EnergyReductionFactor = i;
  }

  /** Sets the stiffness Matrix weight. */
  void      SetElasticity(Float i, unsigned int which = 0)
  {
    m_E[which] = i;
  }

  /** Gets the stiffness Matrix weight. */
  Float     GetElasticity(unsigned int which = 0)
  {
    return m_E[which];
  }

  /** Mass matrix weight */
  void      SetRho(Float r, unsigned int which = 0)
  {
    m_Rho[which] = r;
  }

  /** Image similarity energy weight */
  void      SetGamma(Float r, unsigned int which = 0)
  {
    m_Gamma[which] = r;
  }

  /** Tries to minimize energy */
  void      SetDescentDirectionMinimize()
  {
    m_DescentDirection = positive;
  }

  /** Tries to maximize energy */
  void      SetDescentDirectionMaximize()
  {
    m_DescentDirection = negative;
  }

  /** Finds the minimum energy between the current and next solution
   * by linear search. */
  void      DoLineSearch(unsigned int b)
  {
    m_DoLineSearchOnImageEnergy = b;
  }

  /** Sets the use of multi-resolution strategy.  The control file always uses
    multi-res. */
  void      DoMultiRes(bool b)
  {
    m_DoMultiRes = b;
  }

  /** Sets the use of multi-resolution strategy.  The control file always uses
    multi-res. */
  void      EmployRegridding(unsigned int b)
  {
    m_EmployRegridding = b;
  }

  /** This sets the line search's max iterations. */
  void      SetLineSearchMaximumIterations(unsigned int f)
  {
    m_LineSearchMaximumIterations = f;
  }

  /** Sets the boolean for writing the displacement field to a file. */
  void      SetWriteDisplacements(bool b)
  {
    m_WriteDisplacementField = b;
  }

  /** Sets the boolean for writing the displacement field to a file. */
  bool      GetWriteDisplacements()
  {
    return m_WriteDisplacementField;
  }

  /** Sets the file name for the FEM multi-resolution registration.
      One can also set the parameters in code. */
  void      SetConfigFileName(const char *f){ m_ConfigFileName = f; }

  std::string GetConfigFileName() { return m_ConfigFileName; }

  ImageSizeType GetImageSize(){ return m_FullImageSize; }

  /** Set/Get the Metric.  */
  MetricBaseTypePointer    GetMetric() { return m_Metric; }
  void      SetMetric(MetricBaseTypePointer MP) { m_Metric = MP; }

  /** Choose the metric by parameter : 0= mean squares, 1=cross correlation,
      2=pattern intensity, 3 = mutual information. */
  void      ChooseMetric(float whichmetric);

  /** This function allows one to set the element and its material externally.
    */
  void      SetElement(Element::Pointer e) { m_Element = e; }

  /** This sets the pointer to the material. */
  void      SetMaterial(MaterialType::Pointer m) { m_Material = m; }

  void      PrintVectorField(unsigned int modnum = 1000);

  void      SetNumLevels(unsigned int i) { m_NumLevels = i; }
  void      SetMaxLevel(unsigned int i) { m_MaxLevel = i; }

  void      SetTemp(Float i) { m_Temp = i; }

  /** de/constructor */
  FEMRegistrationFilter();
  ~FEMRegistrationFilter();

// HELPER FUNCTIONS
protected:

  /**
   * \class FEMOF
   * A non-templated class to access FEMObjectFactory
   * Easy access to the FEMObjectFactory. We create a new class
   * whose name is shorter and it's not templated...
   * \ingroup ITK-FEMRegistration
   */
  class FEMOF:public FEMObjectFactory< FEMLightObject >
  {
protected:
    FEMOF();
    ~FEMOF();
  };

  /** This function generates a regular mesh of ElementsPerSide^D size */
  void      CreateMesh(double ElementsPerSide, Solver & S, ImageSizeType sz);

  /** The non-image loads are entered into the solver. */
  void      ApplyLoads(SolverType & S, ImageSizeType Isz, double *spacing = NULL);

  /** The image loads are entered into the solver. */
  void      ApplyImageLoads(SolverType & S, MovingImageType *i1, FixedImageType *i2);

  /**  Builds the itpack linear system wrapper with appropriate parameters.
       Currently undefined */
  void      CreateLinearSystemSolver();

  /** Evaluates the image similarity energy by calling the image metric */
  Float     EvaluateEnergy();

  /** Interpolates the vector field over the domain.
    * Our convention is to always keep the vector field
    * at the scale of the original images.
    */
  void      InterpolateVectorField(SolverType & S);

  /** Calculates the metric over the domain given the vector field.
    */
  FloatImageType *      GetMetricImage(FieldType *F);

  /** Re-size the vector field (smaller to larger). */
  typedef  typename FieldType::Pointer FieldPointer;
  FieldPointer ExpandVectorField(ExpandFactorsType *expandFactors, FieldType *f);

  /** This is used for changing between mesh resolutions. */
  void      SampleVectorFieldAtNodes(SolverType & S);

  Float EvaluateResidual(SolverType & mySolver, Float t);

  /* Finds a triplet that brackets the energy minimum.  From Numerical
    Recipes.*/
  void FindBracketingTriplet(SolverType & mySolver, Float *a, Float *b, Float *c);

  /** Finds the optimum value between the last two solutions
    * and sets the current solution to that value.  Uses Evaluate Residual;
    */
  Float GoldenSection(SolverType & mySolver, Float tol = 0.01, unsigned int MaxIters = 25);

  /** Set the solver's current load. */
//  itkSetMacro( Load, ImageMetricLoadType* );
  itkGetConstMacro(Load, ImageMetricLoadType *);

  void PrintSelf(std::ostream & os, Indent indent) const;

private:

  void InitializeField();

  FEMRegistrationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);        //purposely not implemented

  std::string m_ConfigFileName;
  std::string m_ResultsFileName;
  std::string m_MovingFileName;        // This variable is currently not being
                                       // used.
  std::string m_FixedFileName;         // This variable is currently not being
                                       // used.
  std::string m_LandmarkFileName;
  std::string m_DisplacementsFileName;
  std::string m_MeshFileName;

  unsigned int m_DoLineSearchOnImageEnergy;
  unsigned int m_LineSearchMaximumIterations;

  vnl_vector< unsigned int > m_NumberOfIntegrationPoints;  // resolution of
                                                           // integration
  vnl_vector< unsigned int > m_MetricWidth;
  vnl_vector< unsigned int > m_Maxiters;   // max iterations
  unsigned int               m_TotalIterations;
  unsigned int               m_NumLevels;   // Number of Resolution Levels
  unsigned int               m_MaxLevel;    // Maximum Level (NumLevels is
                                            // original resolution).
  unsigned int m_MeshLevels;                // Number of Mesh Resolutions (
                                            // should be >= 1)
  unsigned int m_MeshStep;                  // Ratio Between Mesh Resolutions (
                                            // currently set to 2, should be >=
                                            // 1)
  unsigned int m_FileCount;                 // keeps track of number of files
                                            // written
  unsigned int m_CurrentLevel;

  typename FixedImageType::SizeType m_CurrentLevelImageSize;

  unsigned int m_WhichMetric;

  /** Stores the number of  pixels per element  of the mesh for each
      resolution of the multi-resolution pyramid */
  vnl_vector< unsigned int > m_MeshPixelsPerElementAtEachResolution;

  Float               m_Dt;          // time step
  vnl_vector< Float > m_E;           // elasticity
  vnl_vector< Float > m_Rho;         // mass matrix weight
  vnl_vector< Float > m_Gamma;       // image similarity weight
  Float               m_Energy;      // current value of energy
  Float               m_MinE;        // minimum recorded energy
  Float               m_MinJacobian; // minimum recorded energy
  Float               m_Alpha;       // difference parameter
  /** Factor we want to reduce the energy by - determines convergence. */
  Float m_EnergyReductionFactor;
  Float m_Temp;

  bool         m_WriteDisplacementField;
  bool         m_DoMultiRes;
  bool         m_UseLandmarks;
  bool         m_ReadMeshFile;
  bool         m_UseMassMatrix;
  unsigned int m_EmployRegridding;
  Sign         m_DescentDirection;

  ImageSizeType m_FullImageSize;   // image size
  ImageSizeType m_ImageOrigin;     // image size
  /** Gives the ratio of original image size to current image size - for dealing
    with multi-res.*/
  ImageSizeType m_ImageScaling;
  ImageSizeType m_CurrentImageScaling;

  typename FieldType::RegionType m_FieldRegion;

  typename FieldType::SizeType m_FieldSize;

  typename FieldType::Pointer m_Field;
  //
  // only use TotalField if re-gridding is employed.
  typename FieldType::Pointer m_TotalField;

  ImageMetricLoadType *m_Load;             // Defines the load to use

  // define the warper
  typename WarperType::Pointer m_Warper;

  // declare a new image to hold the warped  reference
  typename FixedImageType::Pointer m_WarpedImage;
  typename FloatImageType::Pointer m_FloatImage;

  typename FixedImageType::RegionType m_Wregion;

  typename FixedImageType::IndexType m_Windex;

  // declare images for target and reference
  typename MovingImageType::Pointer m_MovingImage;
  typename MovingImageType::Pointer m_OriginalMovingImage;

  typename FixedImageType::Pointer m_FixedImage;

  // element and metric pointers
  typename Element::Pointer m_Element;

  typename MaterialType::Pointer m_Material;

  MetricBaseTypePointer m_Metric;

  LandmarkArrayType   m_LandmarkArray;
  InterpolatorPointer m_Interpolator;
};
}
}  // end namespace itk::fem

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFEMRegistrationFilter.txx"
#endif

#endif
