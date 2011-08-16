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

#include <stdio.h>
#include <ctype.h>
#include "itkMetaArrayWriter.h"
#include "itkMetaArrayReader.h"
#include "itksys/SystemTools.hxx"

int testMetaArray(int argc, char * argv[])
  {
  if (argc > 1)
    {
    itksys::SystemTools::ChangeDirectory(argv[1]);
    }

  std::cout << "Array" << std::endl;
  itk::Array<short> arr;
  arr.SetSize(5);
  arr[0] = 1;
  arr[1] = 2;
  arr[2] = 3;
  arr[3] = 2;
  arr[4] = 1;

  // Write them
  itk::MetaArrayWriter::Pointer arrayWriter = itk::MetaArrayWriter::New();
  arrayWriter->SetFileName("test.mva");
  arrayWriter->SetInput(MET_SHORT, &arr);
  arrayWriter->Update();

  std::cout << "Fixed array" << std::endl;
  itk::FixedArray<short, 5> farr;
  farr[0] = 1;
  farr[1] = 2;
  farr[2] = 3;
  farr[3] = 2;
  farr[4] = 1;
  arrayWriter->SetFileName("test_far.mva");
  arrayWriter->SetInput(MET_SHORT, &farr);
  arrayWriter->Update();

  std::cout << "Vector" << std::endl;
  itk::Vector<float, 5> vec;
  vec[0] = 1;
  vec[1] = 2;
  vec[2] = 3;
  vec[3] = 2;
  vec[4] = 1;
  arrayWriter->SetFileName("test_vec.mva");
  arrayWriter->SetInput(MET_FLOAT, &vec);
  arrayWriter->Update();

  std::cout << "CovariantVector" << std::endl;
  itk::CovariantVector<float, 5> cvec;
  cvec[0] = 1;
  cvec[1] = 2;
  cvec[2] = 3;
  cvec[3] = 2;
  cvec[4] = 1;
  arrayWriter->SetFileName("test_cvec.mvh");
  arrayWriter->SetInput(MET_FLOAT, &cvec);
  arrayWriter->Update();

  std::cout << "VariableLengthVector" << std::endl;
  itk::VariableLengthVector<float> vvec;
  vvec.Reserve(5);
  vvec[0] = 1;
  vvec[1] = 2;
  vvec[2] = 3;
  vvec[3] = 2;
  vvec[4] = 1;
  arrayWriter->SetFileName("test_vvec.mvh");
  arrayWriter->SetInput(MET_FLOAT, &vvec);
  arrayWriter->Update();

  // Read them
  std::cout << "Read VariableLengthVector short" << std::endl;
  itk::VariableLengthVector<short> rvecs;
  itk::MetaArrayReader::Pointer arrayReader = itk::MetaArrayReader::New();
  arrayReader->SetFileName("test.mva");
  arrayReader->Update();
  arrayReader->GetOutput(MET_SHORT, & rvecs);
  std::cout << "  vec short = " << rvecs << std::endl;

  std::cout << "Read VariableLengthVector float" << std::endl;
  itk::VariableLengthVector<float> rvecf;
  arrayReader->GetOutput(MET_FLOAT, & rvecf);
  std::cout << "  rvec float = " << rvecf << std::endl;

  std::cout << "Read fixed array" << std::endl;
  itk::FixedArray<float, 5> farray;
  arrayReader->GetOutput(MET_FLOAT, & farray);
  std::cout << "  fixed array float = " << farray << std::endl;

  std::cout << "Read vector" << std::endl;
  itk::Vector<float, 5> rvecf5;
  arrayReader->GetOutput(MET_FLOAT, & rvecf5);
  std::cout << "  vector float = " << rvecf5 << std::endl;

  std::cout << "Read CovariantVector" << std::endl;
  itk::CovariantVector<float, 5> rcovec;
  arrayReader->GetOutput(MET_FLOAT, & rcovec);
  std::cout << "  covariant vector float = " << rcovec << std::endl;

  return EXIT_SUCCESS;
  }
