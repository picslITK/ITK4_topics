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
#ifndef __itkNumericTraits_h
#define __itkNumericTraits_h

#include "itkMacro.h"

#undef min
#undef max

#define itkNUMERIC_TRAITS_MIN_MAX_MACRO()          \
  static ValueType min()                           \
    {                                              \
    return vcl_numeric_limits< ValueType >::min(); \
    }                                              \
  static ValueType max()                           \
    {                                              \
    return vcl_numeric_limits< ValueType >::max(); \
    }                                              \
  static ValueType min(ValueType)                  \
    {                                              \
    return vcl_numeric_limits< ValueType >::min(); \
    }                                              \
  static ValueType max(ValueType)                  \
    {                                              \
    return vcl_numeric_limits< ValueType >::max(); \
    }                                              \


#include "vcl_limits.h" // for vcl_numeric_limits
#include <complex>

namespace itk
{

// forward decare to avoid circular dependencies
template< typename TValueType, unsigned int VLength>  class FixedArray;

/** \class NumericTraits
 * \brief Define additional traits for native types such as int or float.
 *
 * NumericTraits is used to extend the traits associated with native types
 * such as float, char, int, and so on. These traits are extensions of the
 * standard numeric_limits defined by the C++ compilers. Some of the added
 * traits include minimum and maximum value; accumulation type; etc.
 *
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 * \wikiexample{SimpleOperations/NumericTraits,Get some basic information about a type}
 */
template< class T >
class NumericTraits:public vcl_numeric_limits< T >
{
public:
  /** The type of this limits trait object. */
  typedef vcl_numeric_limits< T > TraitsType;

  /** Return the type of this native type. */
  typedef T ValueType;

  /** Return the type that can be printed. */
  typedef T PrintType;

  /** Return value of vcl_abs(). */
  typedef T AbsType;

  /** Accumulation of addition and multiplication. */
  typedef double AccumulateType;

  /** Measurement vector type */
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  /** Typedef for operations that use floating point instead of real precision
   *  to save memory */
  typedef float FloatType;

  /** Additive identity. */
  static const T Zero;

  /** Multiplicative identity. */
  static const T One;

  /** Smallest (most nonpositive) value */
  static T NonpositiveMin() { return TraitsType::min(); }

  /** Is a given value positive? */
  static bool IsPositive(T val) { return val > Zero; }

  /** Is a given value nonpositive? */
  static bool IsNonpositive(T val) { return val <= Zero; }

  /** Is a given value negative? */
  static bool IsNegative(T val) { return val < Zero; }

  /** Is a given value nonnegative? */
  static bool IsNonnegative(T val) { return val >= Zero; }

  /** Return zero value. This function should be used to support
   *  RGBPixel type and standard types (not vectors) */
  static T ZeroValue() { return Zero; }

  /** Return one value. This function should be used to support
   *  RGBPixel type and standard types (not vectors) */
  static T OneValue() { return One; }

  /* Provide a default implementation of the max() method with
   * argument. This API is needed for VariableLengthVector because
   * its length is only known at run-time. Specializations of the
   * VariableLengthVector will provide a different implementation
   * where a vector of the correct size is built. */
  static T max(const T & val) { return TraitsType::max(); }
  static T min(const T & val) { return TraitsType::min(); }
};

/** \cond HIDE_SPECIALIZATION_DOCUMENTATION */

/** \class NumericTraits<bool>
 * \brief Define traits for type bool.
 *
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */

template< >
class NumericTraits< bool > :public vcl_numeric_limits< bool >
{
public:
  typedef bool                     ValueType;
  typedef bool                     PrintType;
  typedef unsigned char            AbsType;
  typedef unsigned char            AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const bool ITKCommon_EXPORT Zero;
  static const bool ITKCommon_EXPORT One;

  static bool min() { return false; }
  static bool max() { return true; }
  static bool min(bool) { return min(); }
  static bool max(bool) { return max(); }
  static bool NonpositiveMin() { return false; }
  static bool IsPositive(bool val) { return val; }
  static bool IsNonpositive(bool val) { return !val; }
  static bool IsNegative(bool val) { return val ? false : false; }
  static bool IsNonnegative(bool val) { return val ? true : true; }
  static bool ZeroValue() { return Zero; }
  static bool OneValue() { return One; }
};

/** \class NumericTraits<char>
 * \brief Define traits for type char.
 * NOTE: char is not guaranteed to be signed. On SGI's, the default is unsigned
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< char > :public vcl_numeric_limits< char >
{
public:
  typedef char                     ValueType;
  typedef int                      PrintType;
  typedef unsigned char            AbsType;
  typedef short                    AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const char ITKCommon_EXPORT Zero;
  static const char ITKCommon_EXPORT One;

#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable: 4310) // cast truncates constant value
#endif
  static char min() { return char(255) < 0 ? -128 : 0; }
  static char max() { return char(255) < 0 ? 127 : 255; }
#ifdef _MSC_VER
#pragma warning (pop)
#endif

  static char min(char) { return min(); }
  static char max(char) { return max(); }
  static char NonpositiveMin() { return min(); }
  static bool IsPositive(char val) { return val > Zero; }
  static bool IsNonpositive(char val) { return val <= Zero; }
  static bool IsNegative(char val) { return val < Zero; }
  static bool IsNonnegative(char val) { return val >= Zero; }
  static char ZeroValue() { return Zero; }
  static char OneValue() { return One; }
};

/** \class NumericTraits<char>
 * \brief Define traits for type char.
 * NOTE: char is not guaranteed to be signed. On SGI's, the default is unsigned
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< signed char > :public vcl_numeric_limits< signed char >
{
public:
  typedef signed char              ValueType;
  typedef int                      PrintType;
  typedef unsigned char            AbsType;
  typedef short                    AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const signed char ITKCommon_EXPORT Zero;
  static const signed char ITKCommon_EXPORT One;

  static signed char min() { return -128; }
  static signed char max() { return 127; }
  static signed char min(signed char) { return min(); }
  static signed char max(signed char) { return max(); }
  static signed char NonpositiveMin() { return min(); }
  static bool IsPositive(signed char val) { return val > Zero; }
  static bool IsNonpositive(signed char val) { return val <= Zero; }
  static bool IsNegative(signed char val) { return val < Zero; }
  static bool IsNonnegative(signed char val) { return val >= Zero; }
  static signed char  ZeroValue() { return Zero; }
  static signed char OneValue() { return One; }
};

/** \class NumericTraits<unsigned char>
 * \brief Define traits for type unsigned char.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< unsigned char > :public vcl_numeric_limits< unsigned char >
{
public:
  typedef unsigned char            ValueType;
  typedef int                      PrintType;
  typedef unsigned char            AbsType;
  typedef unsigned short           AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const unsigned char ITKCommon_EXPORT Zero;
  static const unsigned char ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();

  static unsigned char NonpositiveMin() { return vcl_numeric_limits< ValueType >::min(); }
  static bool IsPositive(unsigned char val) { return val != Zero; }
  static bool IsNonpositive(unsigned char val) { return val == Zero; }
  static bool IsNegative(unsigned char val) { return val ? false : false; }
  static bool IsNonnegative(unsigned char val) { return val ? true : true; }
  static unsigned char  ZeroValue() { return Zero; }
  static unsigned char OneValue() { return One; }
};

/** \class NumericTraits<short>
 * \brief Define traits for type short.
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< short > :public vcl_numeric_limits< short >
{
public:
  typedef short                    ValueType;
  typedef short                    PrintType;
  typedef unsigned short           AbsType;
  typedef int                      AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const short ITKCommon_EXPORT Zero;
  static const short ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static short NonpositiveMin() { return vcl_numeric_limits< ValueType >::min(); }
  static bool IsPositive(short val) { return val > Zero; }
  static bool IsNonpositive(short val) { return val <= Zero; }
  static bool IsNegative(short val) { return val < Zero; }
  static bool IsNonnegative(short val) { return val >= Zero; }
  static short  ZeroValue() { return Zero; }
  static short OneValue() { return One; }
};

/** \class NumericTraits<unsigned short>
 * \brief Define traits for type unsigned short.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< unsigned short > :public vcl_numeric_limits< unsigned short >
{
public:
  typedef unsigned short           ValueType;
  typedef unsigned short           PrintType;
  typedef unsigned short           AbsType;
  typedef unsigned int             AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const unsigned short ITKCommon_EXPORT Zero;
  static const unsigned short ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static unsigned short NonpositiveMin() { return vcl_numeric_limits< ValueType >::min(); }
  static bool IsPositive(unsigned short val) { return val != Zero; }
  static bool IsNonpositive(unsigned short val) { return val == Zero; }
  static bool IsNegative(unsigned short val) { return val ? false : false; }
  static bool IsNonnegative(unsigned short val) { return val ? true : true; }
  static unsigned short  ZeroValue() { return Zero; }
  static unsigned short OneValue() { return One; }
};

/** \class NumericTraits<int>
 * \brief Define traits for type int.
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< int > :public vcl_numeric_limits< int >
{
public:
  typedef int                      ValueType;
  typedef int                      PrintType;
  typedef unsigned int             AbsType;
  typedef long                     AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const int ITKCommon_EXPORT Zero;
  static const int ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static int NonpositiveMin() { return vcl_numeric_limits< ValueType >::min(); }
  static bool IsPositive(int val) { return val > Zero; }
  static bool IsNonpositive(int val) { return val <= Zero; }
  static bool IsNegative(int val) { return val < Zero; }
  static bool IsNonnegative(int val) { return val >= Zero; }
  static int  ZeroValue() { return Zero; }
  static int OneValue() { return One; }
};

/** \class NumericTraits<unsigned int>
 * \brief Define traits for type unsigned int.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< unsigned int > :public vcl_numeric_limits< unsigned int >
{
public:
  typedef unsigned int             ValueType;
  typedef unsigned int             PrintType;
  typedef unsigned int             AbsType;
  typedef unsigned int             AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const unsigned int ITKCommon_EXPORT Zero;
  static const unsigned int ITKCommon_EXPORT One;

  static unsigned int min(void) { return 0; }
  static unsigned int max(void) { return static_cast< unsigned int >( -1 ); }
  static unsigned int min(unsigned int) { return vcl_numeric_limits< ValueType >::min(); }
  static unsigned int max(unsigned int) { return vcl_numeric_limits< ValueType >::max(); }
  static unsigned int NonpositiveMin() { return 0; }
  static bool IsPositive(unsigned int val) { return val != Zero; }
  static bool IsNonpositive(unsigned int val) { return val == Zero; }
  static bool IsNegative(unsigned int val) { return val ? false : false; }
  static bool IsNonnegative(unsigned int val) { return val ? true : true; }
  static unsigned int  ZeroValue() { return Zero; }
  static unsigned int OneValue() { return One; }
};

/** \class NumericTraits<long>
 * \brief Define traits for type long.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< long > :public vcl_numeric_limits< long >
{
public:
  typedef long                     ValueType;
  typedef long                     PrintType;
  typedef unsigned long            AbsType;
  typedef long                     AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const long ITKCommon_EXPORT Zero;
  static const long ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static long NonpositiveMin() { return vcl_numeric_limits< ValueType >::min(); }
  static bool IsPositive(long val) { return val > Zero; }
  static bool IsNonpositive(long val) { return val <= Zero; }
  static bool IsNegative(long val) { return val < Zero; }
  static bool IsNonnegative(long val) { return val >= Zero; }
  static long  ZeroValue() { return Zero; }
  static long OneValue() { return One; }
};

/** \class NumericTraits<unsigned long>
 * \brief Define traits for type unsigned long.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< unsigned long > :public vcl_numeric_limits< unsigned long >
{
public:
  typedef unsigned long            ValueType;
  typedef unsigned long            PrintType;
  typedef unsigned long            AbsType;
  typedef unsigned long            AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const unsigned long ITKCommon_EXPORT Zero;
  static const unsigned long ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static unsigned long NonpositiveMin() { return vcl_numeric_limits< ValueType >::min(); }
  static bool IsPositive(unsigned long val) { return val != Zero; }
  static bool IsNonpositive(unsigned long val) { return val == Zero; }
  static bool IsNegative(unsigned long) { return false; }
  static bool IsNonnegative(unsigned long) { return true; }
  static unsigned long  ZeroValue() { return Zero; }
  static unsigned long  OneValue() { return One; }
};

/** \class NumericTraits<float>
 * \brief Define traits for type float.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< float > :public vcl_numeric_limits< float >
{
public:
  typedef float                    ValueType;
  typedef float                    PrintType;
  typedef float                    AbsType;
  typedef double                   AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const float ITKCommon_EXPORT Zero;
  static const float ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static float NonpositiveMin() { return -vcl_numeric_limits< ValueType >::max(); }
  static bool IsPositive(float val) { return val > Zero; }
  static bool IsNonpositive(float val) { return val <= Zero; }
  static bool IsNegative(float val) { return val < Zero; }
  static bool IsNonnegative(float val) { return val >= Zero; }
  static float  ZeroValue() { return Zero; }
  static float  OneValue() { return One; }
};

/** \class NumericTraits<double>
 * \brief Define traits for type double.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< double > :public vcl_numeric_limits< double >
{
public:
  typedef double                   ValueType;
  typedef double                   PrintType;
  typedef double                   AbsType;
  typedef double                   AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const double ITKCommon_EXPORT Zero;
  static const double ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static double NonpositiveMin() { return -vcl_numeric_limits< ValueType >::max(); }
  static bool IsPositive(double val) { return val > Zero; }
  static bool IsNonpositive(double val) { return val <= Zero; }
  static bool IsNegative(double val) { return val < Zero; }
  static bool IsNonnegative(double val) { return val >= Zero; }
  static double  ZeroValue() { return Zero; }
  static double  OneValue() { return One; }
};

/** \class NumericTraits<long double>
 * \brief Define traits for type long double.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< long double > :public vcl_numeric_limits< long double >
{
public:
  typedef long double ValueType;
#if defined( __SUNPRO_CC ) && defined( _ILP32 )
  // sun studio in 32 bit mode is unable to print long double values: it
  // segfaults.
  // conversion to double will give usable results if the value is in the double
  // range - better than nothing.
  typedef double                   PrintType;
#else
  typedef long double              PrintType;
#endif
  typedef long double              AbsType;
  typedef long double              AccumulateType;
  typedef long double              RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const long double ITKCommon_EXPORT Zero;
  static const long double ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static long double NonpositiveMin() { return -vcl_numeric_limits< ValueType >::max(); }
  static bool IsPositive(long double val) { return val > Zero; }
  static bool IsNonpositive(long double val) { return val <= Zero; }
  static bool IsNegative(long double val) { return val < Zero; }
  static bool IsNonnegative(long double val) { return val >= Zero; }
  static long double ZeroValue() { return Zero; }
  static long double OneValue() { return One; }
};

/** \class NumericTraits< std::complex<float> >
 * \brief Define traits for type std::complex<float>.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< std::complex< float > >
{
public:
  typedef std::complex< float >  TheType;
  typedef float                  ValueType;
  typedef TheType                PrintType;
  typedef double                 AbsType;
  typedef TheType                AccumulateType;
  typedef std::complex< double > RealType;
  typedef double                 ScalarRealType;
  typedef std::complex< float >  FloatType;
  typedef FixedArray<TheType, 1> MeasurementVectorType;

  static const TheType ITKCommon_EXPORT Zero;
  static const TheType ITKCommon_EXPORT One;

  static TheType min(TheType) { return vcl_numeric_limits< ValueType >::min(); }
  static TheType max(TheType) { return vcl_numeric_limits< ValueType >::max(); }
  static TheType NonpositiveMin()
  {
    return TheType(-NumericTraits< float >::NonpositiveMin(), 0.0f);
  }

  static bool IsPositive(TheType val) { return val.real() > 0.0; }
  static bool IsNonpositive(TheType val) { return val.real() <= 0.0; }
  static bool IsNegative(TheType val) { return val.real() < 0.0; }
  static bool IsNonnegative(TheType val) { return val.real() >= 0.0; }
  static TheType ZeroValue() { return Zero; }
  static TheType OneValue() { return One; }
};

/** \class NumericTraits< std::complex<double> >
 * \brief Define traits for type std::complex<double>.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< std::complex< double > >
{
public:
  typedef std::complex< double > TheType;
  typedef double                 ValueType;
  typedef TheType                PrintType;
  typedef double                 AbsType;
  typedef TheType                AccumulateType;
  typedef std::complex< double > RealType;
  typedef double                 ScalarRealType;
  typedef std::complex< float >  FloatType;
  typedef FixedArray<TheType, 1> MeasurementVectorType;

  static const TheType ITKCommon_EXPORT Zero;
  static const TheType ITKCommon_EXPORT One;

  static TheType min(TheType) { return vcl_numeric_limits< ValueType >::min(); }
  static TheType max(TheType) { return vcl_numeric_limits< ValueType >::max(); }
  static TheType NonpositiveMin()
  {
    return TheType(-NumericTraits< double >::NonpositiveMin(), 0.0);
  }

  static bool IsPositive(TheType val) { return val.real() > 0.0; }
  static bool IsNonpositive(TheType val) { return val.real() <= 0.0; }
  static bool IsNegative(TheType val) { return val.real() < 0.0; }
  static bool IsNonnegative(TheType val) { return val.real() >= 0.0; }
  static TheType ZeroValue() { return Zero; }
  static TheType OneValue() { return One; }
};

/** \class NumericTraits<long long>
 * \brief Define traits for type long long.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< long long > :
  public vcl_numeric_limits< long long >
{
public:
  typedef long long                ValueType;
  typedef long long                PrintType;
  typedef long long                AbsType;
  typedef long long                AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const ValueType ITKCommon_EXPORT Zero;
  static const ValueType ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static ValueType NonpositiveMin() { return vcl_numeric_limits< ValueType >::min(); }
  static bool IsPositive(ValueType val) { return val > Zero; }
  static bool IsNonpositive(ValueType val) { return val <= Zero; }
  static bool IsNegative(ValueType val) { return val < Zero; }
  static bool IsNonnegative(ValueType val) { return val >= Zero; }
  static ValueType  ZeroValue() { return Zero; }
  static ValueType  OneValue() { return One; }
};

/** \class NumericTraits<unsigned long long>
 * \brief Define traits for type unsigned long long.
 * \ingroup DataRepresentation
 * \ingroup ITK-Common
 */
template< >
class NumericTraits< unsigned long long > :
  public vcl_numeric_limits< unsigned long long >
{
public:
  typedef unsigned long long       ValueType;
  typedef unsigned long long       PrintType;
  typedef unsigned long long       AbsType;
  typedef unsigned long long       AccumulateType;
  typedef double                   RealType;
  typedef RealType                 ScalarRealType;
  typedef float                    FloatType;
  typedef FixedArray<ValueType, 1> MeasurementVectorType;

  static const ValueType ITKCommon_EXPORT Zero;
  static const ValueType ITKCommon_EXPORT One;

  itkNUMERIC_TRAITS_MIN_MAX_MACRO();
  static ValueType NonpositiveMin() { return vcl_numeric_limits< ValueType >::min(); }
  static bool IsPositive(ValueType val) { return val != Zero; }
  static bool IsNonpositive(ValueType val) { return val == Zero; }
  static bool IsNegative(ValueType) { return false; }
  static bool IsNonnegative(ValueType) { return true; }
  static ValueType ZeroValue() { return Zero; }
  static ValueType OneValue() { return One; }
};

/** \endcond */

} // end namespace itk

#include "itkFixedArray.h"

#endif // __itkNumericTraits_h
