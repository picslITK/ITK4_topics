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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>
#include <vector>

#include "itkMaximumDecisionRule.h"


int itkMaximumDecisionRuleTest(int,char *[] )
{

  std::cout << "==================================" << std::endl;
  std::cout << "Testing MaximumDecionRule " << std::endl << std::endl;

  typedef itk::MaximumDecisionRule  DecisionRuleType;
  DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

  std::vector< double > discriminantScores;
  discriminantScores.resize( 3 );

  discriminantScores[0] = 0.0;
  discriminantScores[1] = 1.0;
  discriminantScores[2] = 2.0;

  if ( decisionRule->Evaluate( discriminantScores ) != 2 )
    {
    std::cout << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  DecisionRuleType::VectorType discriminantScores2;
  discriminantScores2.resize( 3 );

  discriminantScores2[0] = 0.0;
  discriminantScores2[1] = 1.0;
  discriminantScores2[2] = 2.0;

  if ( decisionRule->Evaluate( discriminantScores2 ) != 2 )
    {
    std::cout << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  DecisionRuleType::ArrayType discriminantScores3(3);

  discriminantScores3[0] = 0.0;
  discriminantScores3[1] = 1.0;
  discriminantScores3[2] = 2.0;

  if ( decisionRule->Evaluate( discriminantScores3 ) != 2 )
    {
    std::cout << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "[SUCCEEDED]" << std::endl;
  return EXIT_SUCCESS;
}
