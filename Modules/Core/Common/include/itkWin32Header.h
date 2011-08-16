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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef __itkWin32Header_h
#define __itkWin32Header_h

#include "itkConfigure.h"

/** Disable some common warnings in MS VC++ */
#if defined( _MSC_VER )

// 'conversion' conversion from 'type1' to 'type2', possible loss of data
#pragma warning ( disable : 4244 )

// 'identifier' : truncation from 'type1' to 'type2'
#pragma warning ( disable : 4305 )

// 'conversion' : truncation of constant value
#pragma warning ( disable : 4309 )

// decorated name length exceeded, name was truncated
#pragma warning ( disable : 4503 )

// 'type' : forcing value to bool 'true' or 'false' (performance warning)
#pragma warning ( disable : 4800 )

// 'identifier' : class 'type' needs to have dll-interface to be used by
// clients of class 'type2'
#pragma warning ( disable : 4251 )

// non dll-interface class 'type' used as base for dll-interface class 'type2'
#pragma warning ( disable : 4275 )

// C++ exception specification ignored except to indicate a
// function is not __declspec(nothrow)
#pragma warning ( disable : 4290 )

// 'type' : inconsistent dll linkage.  dllexport assumed.
#pragma warning ( disable : 4273 )

// conditional expression is constant
#pragma warning ( disable : 4127 )

// unreferenced local function has been removed
#pragma warning ( disable : 4505 )

#endif

// As of MSVS++ 7.1 and greater, typename is supported in templates
#define ITK_TYPENAME typename

// When a class definition has ITK_EXPORT, the class will be
// checked automatically, by Utilities/Dart/PrintSelfCheck.tcl
#define ITK_EXPORT

#if ( defined( _WIN32 ) || defined( WIN32 ) ) && !defined( ITKSTATIC )
#ifdef ITKCommon_EXPORTS
#define ITKCommon_EXPORT __declspec(dllexport)
#else
#define ITKCommon_EXPORT __declspec(dllimport)
#endif  /* ITKCommon_EXPORTS */

#else
/* unix needs nothing */
#define ITKCommon_EXPORT
#endif

#endif
