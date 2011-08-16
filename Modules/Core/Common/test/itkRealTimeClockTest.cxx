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
#include "itkRealTimeClock.h"
#include "vcl_cmath.h"

int itkRealTimeClockTest( int, char * [] )
{
  try
    {

    // Create an ITK RealTimeClock
    itk::RealTimeClock::Pointer clock = itk::RealTimeClock::New();

    std::cout << "Testing itk::RealTimeClock" << std::endl;
    std::cout << "Frequency: " << clock->GetFrequency() << std::endl;

    std::cout.precision(30);

    unsigned int i;

    itk::RealTimeStamp timestamps[5];

    std::cout << "Printing timestamps got one by one" << std::endl;

    for( i = 0; i < 5; ++i )
      {
      std::cout << clock->GetRealTimeStamp() << std::endl;
      }

    for( i = 0; i < 5; ++i )
      {
      timestamps[i] = clock->GetRealTimeStamp();
      }

    std::cout << "Printing timestamps buffered" << std::endl;
    for( i = 0; i < 5; ++i )
      {
      std::cout << timestamps[i] << std::endl;
      }

    // Print out several time stamps
    itk::RealTimeStamp realStamp1 = clock->GetRealTimeStamp();
    itk::RealTimeStamp realStamp2 = clock->GetRealTimeStamp();
    std::cout << "Current Time " << realStamp2 << std::endl;

    typedef itk::RealTimeStamp::TimeRepresentationType    TimeRepresentationType;

    TimeRepresentationType tolerance = 1e6;

    for( i = 0; i < 10; ++i )
      {
      realStamp1 = realStamp2;
      realStamp2 = clock->GetRealTimeStamp();
      itk::RealTimeInterval difference = realStamp2 - realStamp1;
      itk::RealTimeStamp::TimeRepresentationType seconds1 = realStamp1.GetTimeInSeconds();
      itk::RealTimeStamp::TimeRepresentationType seconds2 = realStamp2.GetTimeInSeconds();
      itk::RealTimeStamp::TimeRepresentationType secondsD = difference.GetTimeInSeconds();
      itk::RealTimeStamp::TimeRepresentationType secondsE = seconds2 - seconds1;
      std::cout << realStamp2 << " - " << realStamp1 << " = ";
      std::cout << secondsD << " = " << secondsE << std::endl;

      if( vcl_abs( secondsD - secondsE ) / secondsE > tolerance )
        {
        std::cerr << "Precision error in time difference" << std::endl;
        std::cerr << "Expected " << secondsE << " seconds " << std::endl;
        std::cerr << "But got  " << secondsD << " seconds " << std::endl;
        return EXIT_FAILURE;
        }
      }

    }
  catch(...)
    {
    std::cerr << "Exception catched !!" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "[PASSED]" << std::endl;
  return EXIT_SUCCESS;
}


