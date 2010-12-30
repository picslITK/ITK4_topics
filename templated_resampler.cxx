#include <iostream>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <time.h>
#include "itkVector.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTranslationTransform.h"
#include "itkTemplatedResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNoInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

int main(int argc, char *argv[])
{
  typedef double RealType;
  typedef itk::Image<float,3> ImageType;
  typedef itk::TranslationTransform<double,3>  TransformType;
  typedef  itk::Vector<double,3>     VectorType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::LinearInterpolateImageFunction<ImageType> FuncTypeL;
  typedef itk::NoInterpolateImageFunction<ImageType> FuncTypeNo;
  typedef itk::BSplineInterpolateImageFunction<ImageType, RealType>
    FuncType;

  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[1] );
  reader2->Update();

  FuncTypeNo::Pointer noop=FuncTypeNo::New();

  FuncTypeL::Pointer line1=FuncTypeL::New();
  //  line1->SetInputImage( reader1->GetOutput());
  FuncTypeL::Pointer line2=FuncTypeL::New();
  //line2->SetInputImage( reader2->GetOutput());
  
  FuncType::Pointer bSplineInterpolator1
    = FuncType::New();
  bSplineInterpolator1->SetSplineOrder( 3 );
  FuncType::Pointer bSplineInterpolator2
    = FuncType::New();
  bSplineInterpolator2->SetSplineOrder( 3 );


  TransformType::Pointer  trans_tran = TransformType::New();
  VectorType trans;
  trans.Fill(3.5);
  trans_tran->Translate(trans);
  bool write=false;
  clock_t t3=clock();
  typedef itk::ResampleImageFilter<ImageType, ImageType, double> ResamplerType;
  ResamplerType::Pointer resampleFilter = ResamplerType::New();
  resampleFilter->SetNumberOfThreads( 1 );
  resampleFilter->SetInput( reader1->GetOutput() );
  resampleFilter->SetOutputParametersFromImage( reader1->GetOutput() );
  resampleFilter->SetTransform( trans_tran );
  //resampleFilter->SetInterpolator( bSplineInterpolator1 );
  //resampleFilter->SetInterpolator( noop );
  resampleFilter->SetInterpolator( line1 );
  resampleFilter->ReleaseDataFlagOn();
  if (write)
  {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( "test1.nii.gz" );
    writer->SetInput( resampleFilter->GetOutput() );
    writer->Update();
  }
  else resampleFilter->Update();
  clock_t t4=clock();

  std::cout << " T-1 " << std::endl;
  clock_t t1=clock();
  typedef itk::TemplatedResampleImageFilter<ImageType, ImageType, FuncTypeL, TransformType> ResamplerType2;
  ResamplerType2::Pointer resampleFilter2 = ResamplerType2::New();
  resampleFilter2->SetNumberOfThreads( 1 );
  resampleFilter2->SetInput( reader2->GetOutput() );
  resampleFilter2->SetOutputParametersFromImage( reader2->GetOutput() );
  resampleFilter2->SetTransform( trans_tran );
  resampleFilter2->SetInterpolator( line2 );
  //resampleFilter2->SetInterpolator( bSplineInterpolator2 );
  //resampleFilter2->SetInterpolator( noop );
  resampleFilter2->ReleaseDataFlagOn();
  if (write)
  {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( "test2.nii.gz" );
    writer->SetInput( resampleFilter2->GetOutput() );
    writer->Update();
  }
  else resampleFilter2->Update();
  clock_t t2=clock();

  double testtime1=(t2-t1)/1.e5; 
  double testtime2=(t4-t3)/1.e5; 
  std::cout<< " templated-op-time " << testtime1 << " varifunc-op-time " << testtime2 << " %speed-dif " << 100.*(1-testtime1/testtime2) << std::endl;
}

