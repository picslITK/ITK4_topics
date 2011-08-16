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

/**
 * This is a test file for the itkCylinderSpatialObject class.
 */
#include "itkCylinderSpatialObject.h"

int itkCylinderSpatialObjectTest(int, char* [])
{
  typedef itk::CylinderSpatialObject   CylinderType;

  CylinderType::Pointer myCylinder = CylinderType::New();
  double radius = 3.0;
  double height = 12.0;

  std::cout << "Testing Print after construction " << std::endl;
  myCylinder->Print( std::cout );

  std::cout << "Testing radius : ";

  myCylinder->SetRadius(radius);
  double radius2 = myCylinder->GetRadius();
  if(radius2 != radius)
    {
    std::cout << "[FAILURE]" << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "[PASSED]" << std::endl;

  myCylinder->SetHeight(height);

  // Point consistency
  std::cout << "Is Inside: ";
  itk::Point<double,3> in;
  in[0]=1;in[1]=2;in[2]=0;
  itk::Point<double,3> out;
  out[0]=0;out[1]=0;out[2]=3.1;

  if(!myCylinder->IsInside(in))
    {
    std::cout<<"[FAILED]"<<std::endl;
    return EXIT_FAILURE;
    }

  if(myCylinder->IsInside(out))
    {
    std::cout<<"OUT [FAILED]"<<std::endl;
    return EXIT_FAILURE;
    }

  std::cout<<"[PASSED]"<<std::endl;

  std::cout << "ComputeBoundingBox: ";
  myCylinder->ComputeBoundingBox();
  CylinderType::BoundingBoxType * boundingBox = myCylinder->GetBoundingBox();


  if( (boundingBox->GetBounds()[0] != -3 )
    || (boundingBox->GetBounds()[1] != 3 )
    || (boundingBox->GetBounds()[2] != -6 )
    || (boundingBox->GetBounds()[3] != 6 )
    || (boundingBox->GetBounds()[4] != -3 )
    || (boundingBox->GetBounds()[5] != 3 )
    )
   {
   std::cout<<"[FAILED]"<<std::endl;
   return EXIT_FAILURE;
   }

  std::cout << "Testing Print after all modifications " << std::endl;
  myCylinder->Print( std::cout );

  std::cout<<"[PASSED]"<<std::endl;
  return EXIT_SUCCESS;

}
