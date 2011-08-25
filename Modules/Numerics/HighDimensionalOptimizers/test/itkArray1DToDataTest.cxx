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

#include "itkArray1DToData.h"

using namespace itk;

namespace{

class HolderClass;
typedef Array1DToData<HolderClass>                  Array1DToDataType;
typedef Array1DToDataType::IndexRangeType           IndexRangeType;

/*
 * Helper class that holds callback and stores results from threads.
 */
class HolderClass
{
public:
  typedef HolderClass Self;

  /* Callback used during threaded operation.
   * An instance of this class is referenced through
   * holder, which is passed in via the threader's user data. */
  static void ThreadedCallback(
                                const IndexRangeType& rangeForThread,
                                ThreadIdType threadId,
                                Self *holder )
  {
    if( threadId >= holder->m_NumberOfThreads )
      {
      holder->m_GotMoreThreadsThanExpected = true;
      return;
      }
    holder->m_RangeInCallback[threadId] = rangeForThread;
  }

  void Init(ThreadIdType numberOfThreads)
  {
    m_GotMoreThreadsThanExpected = false;
    m_RangeInCallback.resize(numberOfThreads);
    IndexRangeType emptyRange;
    emptyRange.Fill(-1);
    for(ThreadIdType i=0; i<m_RangeInCallback.size(); i++)
      {
      m_RangeInCallback[i] = emptyRange;
      }
  }

  bool                         m_GotMoreThreadsThanExpected;
  std::vector<IndexRangeType>  m_RangeInCallback;
  ThreadIdType                 m_NumberOfThreads;

};

/*
 * Run the actually test
 */
int RunTest( Array1DToDataType::Pointer& threader, ThreadIdType numberOfThreads,
              IndexRangeType& fullRange, HolderClass& holder)
{
  std::cout << "Testing with " << numberOfThreads
            << " threads and range " << fullRange << "..." << std::endl;

  /* Try to set the requested number of threads */
  threader->SetNumberOfThreads( numberOfThreads );
  if( threader->GetNumberOfThreads() != numberOfThreads )
    {
    std::cerr << "Failed setting requested number of threads: "
              << numberOfThreads << std::endl
              << "threader->GetNumberOfThreads(): "
              << threader->GetNumberOfThreads() << std::endl;
    return EXIT_FAILURE;
    }
  holder.m_NumberOfThreads = numberOfThreads;

  /* Reset the holder */
  holder.Init( numberOfThreads );

  /* Tell the threader the whole range over which to split and thread */
  threader->SetOverallIndexRange( fullRange );

  /* Run the threader */
  threader->GenerateData();

  /* Did we somehow get more threads than requested? */
  if( holder.m_GotMoreThreadsThanExpected )
    {
    std::cerr << "Threader received more threads than requested." << std::endl;
    return EXIT_FAILURE;
    }

  /* Did we use as many threads as requested? */
  std::cout << "requested numberOfThreads: " << numberOfThreads << std::endl
            << "actual: threader->GetNumberOfThreadsUsed(): "
            << threader->GetNumberOfThreadsUsed() << std::endl;

  /* Check the results */
  IndexRangeType::IndexValueType previousEndIndex = -1;
  for( ThreadIdType i=0; i < threader->GetNumberOfThreadsUsed(); i++ )
    {
    IndexRangeType subRange = holder.m_RangeInCallback[i];
    /* Check that the sub range was assigned something at all */
    if( subRange[0] == -1 ||
        subRange[1] == -1 )
        {
        std::cerr << "Error: subRange " << i << " is -1: "
                  << subRange[i];
        return EXIT_FAILURE;
        }
    /* Check that we got the begin of the range */
    if( i == 0 && subRange[0] != fullRange[0] )
        {
        std::cerr << "Error: subRange[0][0] should be " << fullRange[0]
                  << ", but it's " << subRange[0] << ".";
        return EXIT_FAILURE;
        }
    /* Check that we got the end of the range */
    if( i == numberOfThreads-1 && subRange[1] != fullRange[1] )
        {
        std::cerr << "Error: subRange[N-1][1] should be " << fullRange[1]
                  << ", but it's " << subRange[1] << ".";
        return EXIT_FAILURE;
        }
    /* Check that the sub-range endings and beginnings are continuous */
    if( i > 0 )
      {
      if( previousEndIndex + 1 != subRange[0] )
        {
        std::cerr << "Error: subRange " << i << " is not continuous with "
                  << "previous subRange." << std::endl
                  << "previousEndIndex: " << previousEndIndex << std::endl
                  << "subRange[0]: " << subRange[0] << std::endl;
        return EXIT_FAILURE;
        }
      }
    previousEndIndex = subRange[1];
    }

  return EXIT_SUCCESS;
}

}//namespace

/*
 * Main test entry function
 */
int itkArray1DToDataTest(int , char* [])
{
  Array1DToDataType::Pointer threader = Array1DToDataType::New();
  HolderClass   holder;

  int result = EXIT_SUCCESS;

  /* Check # of threads */
  std::cout << "GetGlobalMaximumNumberOfThreads: "
            << threader->GetMultiThreader()->GetGlobalMaximumNumberOfThreads()
            << std::endl;
  std::cout << "GetGlobalDefaultNumberOfThreads: "
            << threader->GetMultiThreader()->GetGlobalDefaultNumberOfThreads()
            << std::endl;
  std::cout << "threader->NumberOfThreads(): " << threader->GetNumberOfThreads()
            << std::endl;

  /* Set the callback for the threader to use */
  threader->SetThreadedGenerateData( HolderClass::ThreadedCallback );
  /* Set the data/method holder that threader will pass to callback */
  threader->SetHolder( &holder );

  IndexRangeType fullRange;

  /* Test with single thread */
  fullRange[0] = 0;
  fullRange[1] = 102; //set total range to prime to test uneven division
  ThreadIdType numberOfThreads = 1;
  if( RunTest( threader, numberOfThreads, fullRange, holder ) != EXIT_SUCCESS )
    {
    result = EXIT_FAILURE;
    }

  /* Test with range that doesn't start at 0 */
  fullRange[0] = 2;
  fullRange[1] = 104; //set total range to prime to test uneven division
  numberOfThreads = 1;
  if( RunTest( threader, numberOfThreads, fullRange, holder ) != EXIT_SUCCESS )
    {
    result = EXIT_FAILURE;
    }


  /* Test with multiple threads */
  if( threader->GetMultiThreader()->GetGlobalMaximumNumberOfThreads() > 1 )
    {
    /* Test with default number of threads. */
    fullRange[0] = 6;
    fullRange[1] = 108; //set total range to prime to test uneven division
    numberOfThreads =
      threader->GetMultiThreader()->GetGlobalDefaultNumberOfThreads();
    if( RunTest( threader, numberOfThreads, fullRange, holder )
          != EXIT_SUCCESS )
      {
      result = EXIT_FAILURE;
      }

    /* Test with max number of threads and check that we only used as
     * many as is reasonable. */
    ThreadIdType maxNumberOfThreads =
      threader->GetMultiThreader()->GetGlobalMaximumNumberOfThreads();
    fullRange[0] = 6;
    fullRange[1] = fullRange[0]+maxNumberOfThreads-2;
    if( RunTest( threader, maxNumberOfThreads, fullRange, holder )
          != EXIT_SUCCESS )
      {
      result = EXIT_FAILURE;
      }
    if( threader->GetNumberOfThreadsUsed() != maxNumberOfThreads-1 )
      {
      std::cerr << "Error: Expected to use only " << maxNumberOfThreads-1
                << "threads, but used " << threader->GetNumberOfThreadsUsed()
                << "." << std::endl;
      }
    }
  else
    {
    std::cout << "No multi-threading available. " << std::endl;
    }

  /* Test that malformed range is caught */
  fullRange[0] = 6;
  fullRange[1] = 2;
  numberOfThreads = 1;
  std::cout << "Test with BAD range. Expect to catch exception ... "
            << std::endl;
  bool caught = false;
  try
    {
    RunTest( threader, numberOfThreads, fullRange, holder );
    }
  catch ( ExceptionObject & err )
    {
    caught = true;
    std::cout << err << std::endl << std::endl
              << "...caught expected exception." << std::endl;
    }
  if( ! caught )
    {
    std::cerr << "Exception not caught for malformed range." << std::endl;
    result = EXIT_FAILURE;
    }

  return result;
}
