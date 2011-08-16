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

#include <iostream>
#include "itkPolyLineParametricPath.h"

int itkPolyLineParametricPathTest(int, char* [])
{
  typedef  itk::PolyLineParametricPath<2>  PathType;
  typedef  PathType::InputType             InputType;
  typedef  PathType::IndexType             IndexType;
  typedef  PathType::OffsetType            OffsetType;
  typedef  PathType::VertexType            VertexType;

  bool passed = true;

  InputType   input;
  OffsetType  offset;
  VertexType  v;

  PathType::Pointer path = PathType::New();


  v.Fill(1);
  path->AddVertex(v);

  v[0]=2;
  v[1]=3;
  path->AddVertex(v);

  v[0]=3;
  v[1]=4;
  path->AddVertex(v);

  v[0]=0;
  v[1]=5;
  path->AddVertex(v);

  v.Fill(1);
  path->AddVertex(v);

  std::cout << "Evaluating at 0, 0.5, and 4.0:  " << path->Evaluate(0) << ", "
       << path->Evaluate(0.5) << ", " << path->Evaluate(4.0) << std::endl;

  std::cout << "Evaluating to an index at 0, 0.5, and 1.0:  "
       << path->EvaluateToIndex(0) << ", " << path->EvaluateToIndex(0.5)
       << ", " << path->EvaluateToIndex(1.0) << std::endl;
  if( path->EvaluateToIndex(4.0) != path->EvaluateToIndex(0) )
    {
    std::cout << "PolyLineParametricPathTest:  EvaluateToIndex() Failed" << std::endl;
    passed = false;
    }

  std::cout << "Evaluating the derivative at 0, 0.5, and 1.0:  "
       << path->EvaluateDerivative(0) << ", " << path->EvaluateDerivative(0.5)
       << ", " << path->EvaluateDerivative(1.0) << std::endl;
  if( int(0.5+1000*(path->EvaluateDerivative(0.5))[0]) != 1000 ||
      int(0.5+1000*(path->EvaluateDerivative(0.5))[1]) != 2000 )
    {
    std::cout << "PolyLineParametricPathTest:  EvaluateDerivative() Failed" << std::endl;
    passed = false;
    }

  input = 0;
  offset = path->IncrementInput( input );
  std::cout << "Incrementing the input from 0 to " << input << ":  " << offset << std::endl;

  input = 0.5;
  offset = path->IncrementInput( input );
  std::cout << "Incrementing the input from 0.5 to " << input << ":  " << offset << std::endl;

  if (passed)
    {
    std::cout << "PolyLineParametricPath tests passed" << std::endl;
    return EXIT_SUCCESS;
    }
  else
    {
    std::cout << "PolyLineParametricPath tests failed" << std::endl;
    return EXIT_FAILURE;
    }
}
