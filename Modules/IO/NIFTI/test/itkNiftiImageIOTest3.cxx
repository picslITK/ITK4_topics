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

#include "itksys/SystemTools.hxx"
#include "itkNiftiImageIO.h"
#include "itkNiftiImageIOTest.h"

template <class ScalarType, unsigned TVecLength, unsigned TDimension>
int
TestImageOfVectors(const std::string &fname)
{
  const int dimsize = 2;
  /** Deformation field pixel type. */
  typedef typename itk::Vector<ScalarType,TVecLength> FieldPixelType;

  /** Deformation field type. */
  typedef typename itk::Image<FieldPixelType,TDimension> VectorImageType;

  //
  // swizzle up a random vector image.
  typename VectorImageType::RegionType imageRegion;
  typename VectorImageType::SizeType size;
  typename VectorImageType::IndexType index;
  typename VectorImageType::SpacingType spacing;
  typename VectorImageType::PointType origin;
  typename VectorImageType::DirectionType myDirection;
  myDirection.Fill(0.0);
  // original test case was destined for failure.  NIfTI always writes out 3D
  // orientation.  The only sensible matrices you could pass in would be of the form
  // A B C 0
  // D E F 0
  // E F G 0
  // 0 0 0 1
  // anything in the 4th dimension that didn't follow that form would just come up scrambled.
  //NOTE: Nifti only reports upto 3D images correctly for direction cosigns.  It is implicitly assumed
  //      that the direction for dimensions 4 or greater come diagonal elements including a 1 in the
  //      direction matrix.
  switch(TDimension)
    {
    case 1:
      myDirection[0][0] = -1.0;
      break;
    case 2:
      myDirection[0][1] = 1.0;
      myDirection[1][0] = -1.0;
      break;
    case 3:
      myDirection[0][2] = 1.0;
      myDirection[1][0] = -1.0;
      myDirection[2][1] = 1.0;
      break;
    case 4:
      myDirection[0][2] = 1.0;
      myDirection[1][0] = -1.0;
      myDirection[2][1] = 1.0;
      myDirection[3][3] = 1.0;
      break;
    }

  std::cout << " === Testing VectorLength: " << TVecLength << " Image Dimension " << static_cast<int>(TDimension) << std::endl;
  std::cout << "======================== Initialized Direction" << std::endl;
  std::cout << myDirection << std::endl;

  for(unsigned i = 0; i < TDimension; i++)
    {
    size[i] = dimsize;
    index[i] = 0;
    spacing[i] = 1.0;
    origin[i] = 0;
    }

  imageRegion.SetSize(size);
  imageRegion.SetIndex(index);
  typename VectorImageType::Pointer vi =
    itk::IOTestHelper::AllocateImageFromRegionAndSpacing<VectorImageType>(imageRegion, spacing);
  vi->SetOrigin(origin);
  vi->SetDirection(myDirection);

  typedef itk::ImageRegionIterator<VectorImageType>      IteratorType;
  typedef itk::ImageRegionConstIterator<VectorImageType> ConstIteratorType;

  int dims[7];
  int _index[7];
  for(unsigned i = 0; i < TDimension; i++)
    {
    dims[i] = size[i];
    }
  for(unsigned i = TDimension; i < 7; i++)
    {
    dims[i] = 1;
    }

  int incr_value=0;
  //  for(fillIt.GoToBegin(); !fillIt.IsAtEnd(); ++fillIt)
  for(int l = 0; l < dims[6]; l++)
    {
    _index[6] = l;
    for(int m = 0; m < dims[5]; m++)
      {
      _index[5] = m;
      for(int n = 0; n < dims[4]; n++)
        {
        _index[4] = n;
        for(int p = 0; p < dims[3]; p++)
          {
          _index[3] = p;
          for(int i = 0; i < dims[2]; i++)
            {
            _index[2] = i;
            for(int j = 0; j < dims[1]; j++)
              {
              _index[1] = j;
              for(int k = 0; k < dims[0]; k++)
                {
                _index[0] = k;
                FieldPixelType pixel;
                float lowrange(100.00),highrange(200.00);
                for(unsigned int q = 0; q < TVecLength; q++)
                  {
                  //pixel[q] = randgen.drand32(lowrange,highrange);
                  pixel[q] = incr_value++;
                  lowrange += 100.0;
                  highrange += 100.0;
                  }
                for(unsigned int q = 0; q < TDimension; q++)
                  {
                  index[q] = _index[q];
                  }
                vi->SetPixel(index,pixel);
                }
              }
            }
          }
        }
      }
    }
  try
    {
    itk::IOTestHelper::WriteImage<VectorImageType,itk::NiftiImageIO>(vi,fname);
    }
  catch(itk::ExceptionObject &ex)
    {
    std::string message;
    message = "Problem found while writing image ";
    message += fname; message += "\n";
    message += ex.GetLocation(); message += "\n";
    message += ex.GetDescription(); std::cout << message << std::endl;
    itk::IOTestHelper::Remove(fname.c_str());
    return EXIT_FAILURE;
    }
  //
  // read it back in.
  typename VectorImageType::Pointer readback;
  try
    {
    readback = itk::IOTestHelper::ReadImage<VectorImageType>(fname);
    }
  catch(itk::ExceptionObject &ex)
    {
    std::string message;
    message = "Problem found while reading image ";
    message += fname; message += "\n";
    message += ex.GetLocation(); message += "\n";
    message += ex.GetDescription(); std::cout << message << std::endl;
    itk::IOTestHelper::Remove(fname.c_str());
    return EXIT_FAILURE;
    }
  bool same = true;
  if(readback->GetOrigin() != vi->GetOrigin() )
    {
    std::cout << "Origin is different: " << readback->GetOrigin() << " != " << vi->GetOrigin()  << std::endl;
    same = false;
    }
  if(readback->GetSpacing() != vi->GetSpacing() )
    {
    std::cout << "Spacing is different: " << readback->GetSpacing() << " != " << vi->GetSpacing()  << std::endl;
    same = false;
    }
  for(unsigned int r=0;r<TDimension;r++)
    {
    for(unsigned int c=0;c<TDimension;c++)
      {
      if(vcl_abs(readback->GetDirection()[r][c] - vi->GetDirection()[r][c]) > 1e-7 )
        {
        std::cout << "Direction is different:\n " << readback->GetDirection() << "\n != \n" << vi->GetDirection()  << std::endl;
        same = false;
        break;
        }
      }
    }
  std::cout << "Original vector Image  ?=   vector Image read from disk " << std::endl;
  for(int l = 0; l < dims[6]; l++)
    {
    _index[6] = l;
    for(int m = 0; m < dims[5]; m++)
      {
      _index[5] = m;
      for(int n = 0; n < dims[4]; n++)
        {
        _index[4] = n;
        for(int p = 0; p < dims[3]; p++)
          {
          _index[3] = p;
          for(int i = 0; i < dims[2]; i++)
            {
            _index[2] = i;
            for(int j = 0; j < dims[1]; j++)
              {
              _index[1] = j;
              for(int k = 0; k < dims[0]; k++)
                {
                _index[0] = k;
                FieldPixelType p1,p2;
                for(unsigned int q = 0; q < TDimension; q++)
                  {
                  index[q] = _index[q];
                  }
                p1 = vi->GetPixel(index);
                p2 = readback->GetPixel(index);
                if(p1 != p2)
                  {
                  same = false;
                  std::cout << p1 << " != " << p2 <<  "    ERROR! " << std::endl;
                  }
                else
                  {
                  std::cout << p1 << " == " << p2 << std::endl;
                  }
                }
              }
            }
          }
        }
      }
    }
  if(same)
    {
    itk::IOTestHelper::Remove(fname.c_str());
    }
  else
    {
    std::cout << "Failing image can be found at: " << fname << std::endl;
    }
  return same ? 0 : EXIT_FAILURE;
}

/** Test writing and reading a Vector Image
 */
int itkNiftiImageIOTest3(int ac, char* av[])
{
  //
  // first argument is passing in the writable directory to do all testing
  if(ac > 1)
    {
    char *testdir = *++av;
    itksys::SystemTools::ChangeDirectory(testdir);
    }
  else
    {
    return EXIT_FAILURE;
    }
  int success(0);

  success |= TestImageOfVectors<float,3,1>(std::string("testVectorImage_float_3_1.nii.gz"));
  success |= TestImageOfVectors<float,3,2>(std::string("testVectorImage_float_3_2.nii.gz"));
  success |= TestImageOfVectors<float,3,3>(std::string("testVectorImage_float_3_3.nii.gz"));
  success |= TestImageOfVectors<float,4,3>(std::string("testVectorImage_float_4_3.nii.gz"));
  success |= TestImageOfVectors<float,4,4>(std::string("testVectorImage_float_4_4.nii.gz"));
  success |= TestImageOfVectors<double,3,3>(std::string("testVectorImage_double_3_3.nii.gz"));

  return success;
}
