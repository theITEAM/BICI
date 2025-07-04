//  Copyright John Maddock 2006-7.
//  Copyright Paul A. Bristow 2007.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_GAMMA_HPP
#define BOOST_MATH_SF_GAMMA_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/config.hpp>
#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4127 4701)
//  // For lexical_cast, until fixed in 1.35?
//  // conditional expression is constant &
//  // Potentially uninitialized local variable 'name' used
#endif
#include <boost/lexical_cast.hpp>
#ifdef BOOST_MSVC
# pragma warning(pop)
#endif
#include <boost/math/tools/series.hpp>
#include <boost/math/tools/fraction.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp>
#include <boost/math/special_functions/lanczos.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/detail/igamma_large.hpp>
#include <boost/math/special_functions/detail/unchecked_factorial.hpp>
#include <boost/math/special_functions/detail/lgamma_small.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/assert.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/greater.hpp>

#include <boost/config/no_tr1/cmath.hpp>
#include <algorithm>

#ifdef BOOST_MATH_INSTRUMENT
#include <iostream>
#include <iomanip>
#include <typeinfo>
#endif

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4702) // unreachable code (return after domain_error throw).
# pragma warning(disable: 4127) // conditional expression is constant.
# pragma warning(disable: 4100) // unreferenced formal parameter.
// Several variables made comments,
// but some difficulty as whether referenced on not may depend on macro values.
// So to be safe, 4100 warnings suppressed.
// TODO - revisit this?
#endif

namespace boost{ namespace math{

namespace detail{

template <class T>
inline bool is_odd(T v, const boost::true_type&)
{
   int i = static_cast<int>(v);
   return i&1;
}
template <class T>
inline bool is_odd(T v, const boost::false_type&)
{
   // Oh dear can't cast T to int!
   BOOST_MATH_STD_USING
   T modulus = v - 2 * floor(v/2);
   return static_cast<bool>(modulus != 0);
}
template <class T>
inline bool is_odd(T v)
{
   return is_odd(v, ::boost::is_convertible<T, int>());
}

template <class T>
T sinpx(T z)
{
   // Ad hoc function calculates x * sin(pi * x),
   // taking extra care near when x is near a whole number.
   BOOST_MATH_STD_USING
   int sign = 1;
   if(z < 0)
   {
      z = -z;
   }
   else
   {
      sign = -sign;
   }
   T fl = floor(z);
   T dist;
   if(is_odd(fl))
   {
      fl += 1;
      dist = fl - z;
      sign = -sign;
   }
   else
   {
      dist = z - fl;
   }
   BOOST_ASSERT(fl >= 0);
   if(dist > 0.5)
      dist = 1 - dist;
   T result = sin(dist*boost::math::constants::pi<T>());
   return sign*z*result;
} // template <class T> T sinpx(T z)
//
// tgamma(z), with Lanczos support:
//
template <class T, class Policy, class L>
T gamma_imp(T z, const Policy& pol, const L& l)
{
   BOOST_MATH_STD_USING

   T result = 1;

#ifdef BOOST_MATH_INSTRUMENT
   static bool b = false;
   if(!b)
   {
      std::cout << "tgamma_imp called with " << typeid(z).name() << " " << typeid(l).name() << std::endl;
      b = true;
   }
#endif
   static const char* function = "boost::math::tgamma<%1%>(%1%)";

   if(z <= 0)
   {
      if(floor(z) == z)
         return policies::raise_pole_error<T>(function, "Evaluation of tgamma at a negative integer %1%.", z, pol);
      if(z <= -20)
      {
         result = gamma_imp(T(-z), pol, l) * sinpx(z);
         if((fabs(result) < 1) && (tools::max_value<T>() * fabs(result) < boost::math::constants::pi<T>()))
            return policies::raise_overflow_error<T>(function, "Result of tgamma is too large to represent.", pol);
         result = -boost::math::constants::pi<T>() / result;
         if(result == 0)
            return policies::raise_underflow_error<T>(function, "Result of tgamma is too small to represent.", pol);
         if((boost::math::fpclassify)(result) == (int)FP_SUBNORMAL)
            return policies::raise_denorm_error<T>(function, "Result of tgamma is denormalized.", result, pol);
         return result;
      }

      // shift z to > 1:
      while(z < 0)
      {
         result /= z;
         z += 1;
      }
   }
   if((floor(z) == z) && (z < max_factorial<T>::value))
   {
      result *= unchecked_factorial<T>(itrunc(z, pol) - 1);
   }
   else
   {
      result *= L::lanczos_sum(z);
      if(z * log(z) > tools::log_max_value<T>())
      {
         // we're going to overflow unless this is done with care:
         T zgh = (z + static_cast<T>(L::g()) - boost::math::constants::half<T>());
         if(log(zgh) * z / 2 > tools::log_max_value<T>())
            return policies::raise_overflow_error<T>(function, "Result of tgamma is too large to represent.", pol);
         T hp = pow(zgh, (z / 2) - T(0.25));
         result *= hp / exp(zgh);
         if(tools::max_value<T>() / hp < result)
            return policies::raise_overflow_error<T>(function, "Result of tgamma is too large to represent.", pol);
         result *= hp;
      }
      else
      {
         T zgh = (z + static_cast<T>(L::g()) - boost::math::constants::half<T>());
         result *= pow(zgh, z - boost::math::constants::half<T>()) / exp(zgh);
      }
   }
   return result;
}
//
// lgamma(z) with Lanczos support:
//
template <class T, class Policy, class L>
T lgamma_imp(T z, const Policy& pol, const L& l, int* sign = 0)
{
#ifdef BOOST_MATH_INSTRUMENT
   static bool b = false;
   if(!b)
   {
      std::cout << "lgamma_imp called with " << typeid(z).name() << " " << typeid(l).name() << std::endl;
      b = true;
   }
#endif

   BOOST_MATH_STD_USING

   static const char* function = "boost::math::lgamma<%1%>(%1%)";

   T result = 0;
   int sresult = 1;
   if(z <= 0)
   {
      // reflection formula:
      if(floor(z) == z)
         return policies::raise_pole_error<T>(function, "Evaluation of lgamma at a negative integer %1%.", z, pol);

      T t = sinpx(z);
      z = -z;
      if(t < 0)
      {
         t = -t;
      }
      else
      {
         sresult = -sresult;
      }
      result = log(boost::math::constants::pi<T>()) - lgamma_imp(z, pol, l) - log(t);
   }
   else if(z < 15)
   {
      typedef typename policies::precision<T, Policy>::type precision_type;
      typedef typename mpl::if_<
         mpl::and_<
            mpl::less_equal<precision_type, mpl::int_<64> >, 
            mpl::greater<precision_type, mpl::int_<0> > 
         >,
         mpl::int_<64>,
         typename mpl::if_<
            mpl::and_<
               mpl::less_equal<precision_type, mpl::int_<113> >,
               mpl::greater<precision_type, mpl::int_<0> > 
            >,
            mpl::int_<113>, mpl::int_<0> >::type
          >::type tag_type;
      result = lgamma_small_imp<T>(z, z - 1, z - 2, tag_type(), pol, l);
   }
   else if((z >= 3) && (z < 100))
   {
      // taking the log of tgamma reduces the error, no danger of overflow here:
      result = log(gamma_imp(z, pol, l));
   }
   else
   {
      // regular evaluation:
      T zgh = static_cast<T>(z + L::g() - boost::math::constants::half<T>());
      result = log(zgh) - 1;
      result *= z - 0.5f;
      result += log(L::lanczos_sum_expG_scaled(z));
   }

   if(sign)
      *sign = sresult;
   return result;
}

//
// Incomplete gamma functions follow:
//
template <class T>
struct upper_incomplete_gamma_fract
{
private:
   T z, a;
   int k;
public:
   typedef std::pair<T,T> result_type;

   upper_incomplete_gamma_fract(T a1, T z1)
      : z(z1-a1+1), a(a1), k(0)
   {
   }

   result_type operator()()
   {
      ++k;
      z += 2;
      return result_type(k * (a - k), z);
   }
};

template <class T>
inline T upper_gamma_fraction(T a, T z, T eps)
{
   // Multiply result by z^a * e^-z to get the full
   // upper incomplete integral.  Divide by tgamma(z)
   // to normalise.
   upper_incomplete_gamma_fract<T> f(a, z);
   return 1 / (z - a + 1 + boost::math::tools::continued_fraction_a(f, eps));
}

template <class T>
struct lower_incomplete_gamma_series
{
private:
   T a, z, result;
public:
   typedef T result_type;
   lower_incomplete_gamma_series(T a1, T z1) : a(a1), z(z1), result(1){}

   T operator()()
   {
      T r = result;
      a += 1;
      result *= z/a;
      return r;
   }
};

template <class T, class Policy>
inline T lower_gamma_series(T a, T z, const Policy& pol, T init_value = 0)
{
   // Multiply result by ((z^a) * (e^-z) / a) to get the full
   // lower incomplete integral. Then divide by tgamma(a)
   // to get the normalised value.
   lower_incomplete_gamma_series<T> s(a, z);
   boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
   T factor = policies::get_epsilon<T, Policy>();
   T result = boost::math::tools::sum_series(s, factor, max_iter, init_value);
   policies::check_series_iterations("boost::math::detail::lower_gamma_series<%1%>(%1%)", max_iter, pol);
   return result;
}

//
// Fully generic tgamma and lgamma use the incomplete partial
// sums added together:
//
template <class T, class Policy>
T gamma_imp(T z, const Policy& pol, const lanczos::undefined_lanczos& l)
{
   static const char* function = "boost::math::tgamma<%1%>(%1%)";
   BOOST_MATH_STD_USING
   if((z <= 0) && (floor(z) == z))
      return policies::raise_pole_error<T>(function, "Evaluation of tgamma at a negative integer %1%.", z, pol);
   if(z <= -20)
   {
      T result = gamma_imp(-z, pol, l) * sinpx(z);
      if((fabs(result) < 1) && (tools::max_value<T>() * fabs(result) < boost::math::constants::pi<T>()))
         return policies::raise_overflow_error<T>(function, "Result of tgamma is too large to represent.", pol);
      result = -boost::math::constants::pi<T>() / result;
      if(result == 0)
         return policies::raise_underflow_error<T>(function, "Result of tgamma is too small to represent.", pol);
      if((boost::math::fpclassify)(result) == (int)FP_SUBNORMAL)
         return policies::raise_denorm_error<T>(function, "Result of tgamma is denormalized.", result, pol);
      return result;
   }
   //
   // The upper gamma fraction is *very* slow for z < 6, actually it's very
   // slow to converge everywhere but recursing until z > 6 gets rid of the
   // worst of it's behaviour.
   //
   T prefix = 1;
   while(z < 6)
   {
      prefix /= z;
      z += 1;
   }
   BOOST_MATH_INSTRUMENT_CODE(prefix);
   if((floor(z) == z) && (z < max_factorial<T>::value))
   {
      prefix *= unchecked_factorial<T>(itrunc(z, pol) - 1);
   }
   else
   {
      prefix = prefix * pow(z / boost::math::constants::e<T>(), z);
      BOOST_MATH_INSTRUMENT_CODE(prefix);
      T sum = detail::lower_gamma_series(z, z, pol) / z;
      BOOST_MATH_INSTRUMENT_CODE(sum);
      sum += detail::upper_gamma_fraction(z, z, ::boost::math::policies::get_epsilon<T, Policy>());
      BOOST_MATH_INSTRUMENT_CODE(sum);
      if(fabs(tools::max_value<T>() / prefix) < fabs(sum))
         return policies::raise_overflow_error<T>(function, "Result of tgamma is too large to represent.", pol);
      BOOST_MATH_INSTRUMENT_CODE((sum * prefix));
      return sum * prefix;
   }
   return prefix;
}

template <class T, class Policy>
T lgamma_imp(T z, const Policy& pol, const lanczos::undefined_lanczos& l, int*sign)
{
   BOOST_MATH_STD_USING

   static const char* function = "boost::math::lgamma<%1%>(%1%)";
   T result = 0;
   int sresult = 1;
   if(z <= 0)
   {
      if(floor(z) == z)
         return policies::raise_pole_error<T>(function, "Evaluation of tgamma at a negative integer %1%.", z, pol);
      T t = detail::sinpx(z);
      z = -z;
      if(t < 0)
      {
         t = -t;
      }
      else
      {
         sresult = -sresult;
      }
      result = log(boost::math::constants::pi<T>()) - lgamma_imp(z, pol, l, 0) - log(t);
   }
   else if((z != 1) && (z != 2))
   {
      T limit = (std::max)(z+1, T(10));
      T prefix = z * log(limit) - limit;
      T sum = detail::lower_gamma_series(z, limit, pol) / z;
      sum += detail::upper_gamma_fraction(z, limit, ::boost::math::policies::get_epsilon<T, Policy>());
      result = log(sum) + prefix;
   }
   if(sign)
      *sign = sresult;
   return result;
}
//
// This helper calculates tgamma(dz+1)-1 without cancellation errors,
// used by the upper incomplete gamma with z < 1:
//
template <class T, class Policy, class L>
T tgammap1m1_imp(T dz, Policy const& pol, const L& l)
{
   BOOST_MATH_STD_USING

   typedef typename policies::precision<T,Policy>::type precision_type;

   typedef typename mpl::if_<
      mpl::or_<
         mpl::less_equal<precision_type, mpl::int_<0> >,
         mpl::greater<precision_type, mpl::int_<113> >
      >,
      typename mpl::if_<
         is_same<L, lanczos::lanczos24m113>,
         mpl::int_<113>,
         mpl::int_<0>
      >::type,
      typename mpl::if_<
         mpl::less_equal<precision_type, mpl::int_<64> >,
         mpl::int_<64>, mpl::int_<113> >::type
       >::type tag_type;

   T result;
   if(dz < 0)
   {
      if(dz < -0.5)
      {
         // Best method is simply to subtract 1 from tgamma:
         result = boost::math::tgamma(1+dz, pol) - 1;
         BOOST_MATH_INSTRUMENT_CODE(result);
      }
      else
      {
         // Use expm1 on lgamma:
         result = boost::math::expm1(-boost::math::log1p(dz, pol) 
            + lgamma_small_imp<T>(dz+2, dz + 1, dz, tag_type(), pol, l));
         BOOST_MATH_INSTRUMENT_CODE(result);
      }
   }
   else
   {
      if(dz < 2)
      {
         // Use expm1 on lgamma:
         result = boost::math::expm1(lgamma_small_imp<T>(dz+1, dz, dz-1, tag_type(), pol, l), pol);
         BOOST_MATH_INSTRUMENT_CODE(result);
      }
      else
      {
         // Best method is simply to subtract 1 from tgamma:
         result = boost::math::tgamma(1+dz, pol) - 1;
         BOOST_MATH_INSTRUMENT_CODE(result);
      }
   }

   return result;
}

template <class T, class Policy>
inline T tgammap1m1_imp(T dz, Policy const& pol,
                 const ::boost::math::lanczos::undefined_lanczos& l)
{
   BOOST_MATH_STD_USING // ADL of std names
   //
   // There should be a better solution than this, but the
   // algebra isn't easy for the general case....
   // Start by subracting 1 from tgamma:
   //
   T result = gamma_imp(1 + dz, pol, l) - 1;
   BOOST_MATH_INSTRUMENT_CODE(result);
   //
   // Test the level of cancellation error observed: we loose one bit
   // for each power of 2 the result is less than 1.  If we would get
   // more bits from our most precise lgamma rational approximation, 
   // then use that instead:
   //
   BOOST_MATH_INSTRUMENT_CODE((dz > -0.5));
   BOOST_MATH_INSTRUMENT_CODE((dz < 2));
   BOOST_MATH_INSTRUMENT_CODE((ldexp(1.0, boost::math::policies::digits<T, Policy>()) * fabs(result) < 1e34));
   if((dz > -0.5) && (dz < 2) && (ldexp(1.0, boost::math::policies::digits<T, Policy>()) * fabs(result) < 1e34))
   {
      result = tgammap1m1_imp(dz, pol, boost::math::lanczos::lanczos24m113());
      BOOST_MATH_INSTRUMENT_CODE(result);
   }
   return result;
}

//
// Series representation for upper fraction when z is small:
//
template <class T>
struct small_gamma2_series
{
   typedef T result_type;

   small_gamma2_series(T a_, T x_) : result(-x_), x(-x_), apn(a_+1), n(1){}

   T operator()()
   {
      T r = result / (apn);
      result *= x;
      result /= ++n;
      apn += 1;
      return r;
   }

private:
   T result, x, apn;
   int n;
};
//
// calculate power term prefix (z^a)(e^-z) used in the non-normalised
// incomplete gammas:
//
template <class T, class Policy>
T full_igamma_prefix(T a, T z, const Policy& pol)
{
   BOOST_MATH_STD_USING

   T prefix;
   T alz = a * log(z);

   if(z >= 1)
   {
      if((alz < tools::log_max_value<T>()) && (-z > tools::log_min_value<T>()))
      {
         prefix = pow(z, a) * exp(-z);
      }
      else if(a >= 1)
      {
         prefix = pow(z / exp(z/a), a);
      }
      else
      {
         prefix = exp(alz - z);
      }
   }
   else
   {
      if(alz > tools::log_min_value<T>())
      {
         prefix = pow(z, a) * exp(-z);
      }
      else if(z/a < tools::log_max_value<T>())
      {
         prefix = pow(z / exp(z/a), a);
      }
      else
      {
         prefix = exp(alz - z);
      }
   }
   //
   // This error handling isn't very good: it happens after the fact
   // rather than before it...
   //
   if((boost::math::fpclassify)(prefix) == (int)FP_INFINITE)
      policies::raise_overflow_error<T>("boost::math::detail::full_igamma_prefix<%1%>(%1%, %1%)", "Result of incomplete gamma function is too large to represent.", pol);

   return prefix;
}
//
// Compute (z^a)(e^-z)/tgamma(a)
// most if the error occurs in this function:
//
template <class T, class Policy, class L>
T regularised_gamma_prefix(T a, T z, const Policy& pol, const L& l)
{
   BOOST_MATH_STD_USING
   T agh = a + static_cast<T>(L::g()) - T(0.5);
   T prefix;
   T d = ((z - a) - static_cast<T>(L::g()) + T(0.5)) / agh;

   if(a < 1)
   {
      //
      // We have to treat a < 1 as a special case because our Lanczos
      // approximations are optimised against the factorials with a > 1,
      // and for high precision types especially (128-bit reals for example)
      // very small values of a can give rather eroneous results for gamma
      // unless we do this:
      //
      // TODO: is this still required?  Lanczos approx should be better now?
      //
      if(z <= tools::log_min_value<T>())
      {
         // Oh dear, have to use logs, should be free of cancellation errors though:
         return exp(a * log(z) - z - lgamma_imp(a, pol, l));
      }
      else
      {
         // direct calculation, no danger of overflow as gamma(a) < 1/a
         // for small a.
         return pow(z, a) * exp(-z) / gamma_imp(a, pol, l);
      }
   }
   else if((fabs(d*d*a) <= 100) && (a > 150))
   {
      // special case for large a and a ~ z.
      prefix = a * boost::math::log1pmx(d, pol) + z * static_cast<T>(0.5 - L::g()) / agh;
      prefix = exp(prefix);
   }
   else
   {
      //
      // general case.
      // direct computation is most accurate, but use various fallbacks
      // for different parts of the problem domain:
      //
      T alz = a * log(z / agh);
      T amz = a - z;
      if(((std::min)(alz, amz) <= tools::log_min_value<T>()) || ((std::max)(alz, amz) >= tools::log_max_value<T>()))
      {
         T amza = amz / a;
         if(((std::min)(alz, amz)/2 > tools::log_min_value<T>()) && ((std::max)(alz, amz)/2 < tools::log_max_value<T>()))
         {
            // compute square root of the result and then square it:
            T sq = pow(z / agh, a / 2) * exp(amz / 2);
            prefix = sq * sq;
         }
         else if(((std::min)(alz, amz)/4 > tools::log_min_value<T>()) && ((std::max)(alz, amz)/4 < tools::log_max_value<T>()) && (z > a))
         {
            // compute the 4th root of the result then square it twice:
            T sq = pow(z / agh, a / 4) * exp(amz / 4);
            prefix = sq * sq;
            prefix *= prefix;
         }
         else if((amza > tools::log_min_value<T>()) && (amza < tools::log_max_value<T>()))
         {
            prefix = pow((z * exp(amza)) / agh, a);
         }
         else
         {
            prefix = exp(alz + amz);
         }
      }
      else
      {
         prefix = pow(z / agh, a) * exp(amz);
      }
   }
   prefix *= sqrt(agh / boost::math::constants::e<T>()) / L::lanczos_sum_expG_scaled(a);
   return prefix;
}
//
// And again, without Lanczos support:
//
template <class T, class Policy>
T regularised_gamma_prefix(T a, T z, const Policy& pol, const lanczos::undefined_lanczos&)
{
   BOOST_MATH_STD_USING

   T limit = (std::max)(T(10), a);
   T sum = detail::lower_gamma_series(a, limit, pol) / a;
   sum += detail::upper_gamma_fraction(a, limit, ::boost::math::policies::get_epsilon<T, Policy>());

   if(a < 10)
   {
      // special case for small a:
      T prefix = pow(z / 10, a);
      prefix *= exp(10-z);
      if(0 == prefix)
      {
         prefix = pow((z * exp((10-z)/a)) / 10, a);
      }
      prefix /= sum;
      return prefix;
   }

   T zoa = z / a;
   T amz = a - z;
   T alzoa = a * log(zoa);
   T prefix;
   if(((std::min)(alzoa, amz) <= tools::log_min_value<T>()) || ((std::max)(alzoa, amz) >= tools::log_max_value<T>()))
   {
      T amza = amz / a;
      if((amza <= tools::log_min_value<T>()) || (amza >= tools::log_max_value<T>()))
      {
         prefix = exp(alzoa + amz);
      }
      else
      {
         prefix = pow(zoa * exp(amza), a);
      }
   }
   else
   {
      prefix = pow(zoa, a) * exp(amz);
   }
   prefix /= sum;
   return prefix;
}
//
// Upper gamma fraction for very small a:
//
template <class T, class Policy>
inline T tgamma_small_upper_part(T a, T x, const Policy& pol, T* pgam = 0, bool invert = false, T* pderivative = 0)
{
   BOOST_MATH_STD_USING  // ADL of std functions.
   //
   // Compute the full upper fraction (Q) when a is very small:
   //
   T result;
   result = boost::math::tgamma1pm1(a, pol);
   if(pgam)
      *pgam = (result + 1) / a;
   T p = boost::math::powm1(x, a, pol);
   result -= p;
   result /= a;
   detail::small_gamma2_series<T> s(a, x);
   boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>() - 10;
   p += 1;
   if(pderivative)
      *pderivative = p / (*pgam * exp(x));
   T init_value = invert ? *pgam : 0;
   result = -p * tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter, (init_value - result) / p);
   policies::check_series_iterations("boost::math::tgamma_small_upper_part<%1%>(%1%, %1%)", max_iter, pol);
   if(invert)
      result = -result;
   return result;
}
//
// Upper gamma fraction for integer a:
//
template <class T, class Policy>
inline T finite_gamma_q(T a, T x, Policy const& pol, T* pderivative = 0)
{
   //
   // Calculates normalised Q when a is an integer:
   //
   BOOST_MATH_STD_USING
   T e = exp(-x);
   T sum = e;
   if(sum != 0)
   {
      T term = sum;
      for(unsigned n = 1; n < a; ++n)
      {
         term /= n;
         term *= x;
         sum += term;
      }
   }
   if(pderivative)
   {
      *pderivative = e * pow(x, a) / boost::math::unchecked_factorial<T>(itrunc(T(a - 1), pol));
   }
   return sum;
}
//
// Upper gamma fraction for half integer a:
//
template <class T, class Policy>
T finite_half_gamma_q(T a, T x, T* p_derivative, const Policy& pol)
{
   //
   // Calculates normalised Q when a is a half-integer:
   //
   BOOST_MATH_STD_USING
   T e = boost::math::erfc(sqrt(x), pol);
   if((e != 0) && (a > 1))
   {
      T term = exp(-x) / sqrt(constants::pi<T>() * x);
      term *= x;
      static const T half = T(1) / 2;
      term /= half;
      T sum = term;
      for(unsigned n = 2; n < a; ++n)
      {
         term /= n - half;
         term *= x;
         sum += term;
      }
      e += sum;
      if(p_derivative)
      {
         *p_derivative = 0;
      }
   }
   else if(p_derivative)
   {
      // We'll be dividing by x later, so calculate derivative * x:
      *p_derivative = sqrt(x) * exp(-x) / constants::root_pi<T>();
   }
   return e;
}
//
// Main incomplete gamma entry point, handles all four incomplete gamma's:
//
template <class T, class Policy>
T gamma_incomplete_imp(T a, T x, bool normalised, bool invert, 
                       const Policy& pol, T* p_derivative)
{
   static const char* function = "boost::math::gamma_p<%1%>(%1%, %1%)";
   if(a <= 0)
      policies::raise_domain_error<T>(function, "Argument a to the incomplete gamma function must be greater than zero (got a=%1%).", a, pol);
   if(x < 0)
      policies::raise_domain_error<T>(function, "Argument x to the incomplete gamma function must be >= 0 (got x=%1%).", x, pol);

   BOOST_MATH_STD_USING

   typedef typename lanczos::lanczos<T, Policy>::type lanczos_type;

   T result;

   BOOST_ASSERT((p_derivative == 0) || (normalised == true));

   bool is_int, is_half_int;
   bool is_small_a = (a < 30) && (a <= x + 1);
   if(is_small_a)
   {
      T fa = floor(a);
      is_int = (fa == a);
      is_half_int = is_int ? false : (fabs(fa - a) == 0.5f);
   }
   else
   {
      is_int = is_half_int = false;
   }

   int eval_method;
   
   if(is_int && (x > 0.6))
   {
      // calculate Q via finite sum:
      invert = !invert;
      eval_method = 0;
   }
   else if(is_half_int && (x > 0.2))
   {
      // calculate Q via finite sum for half integer a:
      invert = !invert;
      eval_method = 1;
   }
   else if(x < 0.5)
   {
      //
      // Changeover criterion chosen to give a changeover at Q ~ 0.33
      //
      if(-0.4 / log(x) < a)
      {
         eval_method = 2;
      }
      else
      {
         eval_method = 3;
      }
   }
   else if(x < 1.1)
   {
      //
      // Changover here occurs when P ~ 0.75 or Q ~ 0.25:
      //
      if(x * 0.75f < a)
      {
         eval_method = 2;
      }
      else
      {
         eval_method = 3;
      }
   }
   else
   {
      //
      // Begin by testing whether we're in the "bad" zone
      // where the result will be near 0.5 and the usual
      // series and continued fractions are slow to converge:
      //
      bool use_temme = false;
      if(normalised && std::numeric_limits<T>::is_specialized && (a > 20))
      {
         T sigma = fabs((x-a)/a);
         if((a > 200) && (policies::digits<T, Policy>() <= 113))
         {
            //
            // This limit is chosen so that we use Temme's expansion
            // only if the result would be larger than about 10^-6.
            // Below that the regular series and continued fractions
            // converge OK, and if we use Temme's method we get increasing
            // errors from the dominant erfc term as it's (inexact) argument
            // increases in magnitude.
            //
            if(20 / a > sigma * sigma)
               use_temme = true;
         }
         else if(policies::digits<T, Policy>() <= 64)
         {
            // Note in this zone we can't use Temme's expansion for 
            // types longer than an 80-bit real:
            // it would require too many terms in the polynomials.
            if(sigma < 0.4)
               use_temme = true;
         }
      }
      if(use_temme)
      {
         eval_method = 5;
      }
      else
      {
         //
         // Regular case where the result will not be too close to 0.5.
         //
         // Changeover here occurs at P ~ Q ~ 0.5
         // Note that series computation of P is about x2 faster than continued fraction
         // calculation of Q, so try and use the CF only when really necessary, especially
         // for small x.
         //
         if(x - (1 / (3 * x)) < a)
         {
            eval_method = 2;
         }
         else
         {
            eval_method = 4;
            invert = !invert;
         }
      }
   }

   switch(eval_method)
   {
   case 0:
      {
         result = finite_gamma_q(a, x, pol, p_derivative);
         if(normalised == false)
            result *= boost::math::tgamma(a, pol);
         break;
      }
   case 1:
      {
         result = finite_half_gamma_q(a, x, p_derivative, pol);
         if(normalised == false)
            result *= boost::math::tgamma(a, pol);
         if(p_derivative && (*p_derivative == 0))
            *p_derivative = regularised_gamma_prefix(a, x, pol, lanczos_type());
         break;
      }
   case 2:
      {
         // Compute P:
         result = normalised ? regularised_gamma_prefix(a, x, pol, lanczos_type()) : full_igamma_prefix(a, x, pol);
         if(p_derivative)
            *p_derivative = result;
         if(result != 0)
         {
            T init_value = 0;
            if(invert)
            {
               init_value = -a * (normalised ? 1 : boost::math::tgamma(a, pol)) / result;
            }
            result *= detail::lower_gamma_series(a, x, pol, init_value) / a;
            if(invert)
            {
               invert = false;
               result = -result;
            }
         }
         break;
      }
   case 3:
      {
         // Compute Q:
         invert = !invert;
         T g;
         result = tgamma_small_upper_part(a, x, pol, &g, invert, p_derivative);
         invert = false;
         if(normalised)
            result /= g;
         break;
      }
   case 4:
      {
         // Compute Q:
         result = normalised ? regularised_gamma_prefix(a, x, pol, lanczos_type()) : full_igamma_prefix(a, x, pol);
         if(p_derivative)
            *p_derivative = result;
         if(result != 0)
            result *= upper_gamma_fraction(a, x, policies::get_epsilon<T, Policy>());
         break;
      }
   case 5:
      {
         //
         // Use compile time dispatch to the appropriate
         // Temme asymptotic expansion.  This may be dead code
         // if T does not have numeric limits support, or has
         // too many digits for the most precise version of
         // these expansions, in that case we'll be calling
         // an empty function.
         //
         typedef typename policies::precision<T, Policy>::type precision_type;

         typedef typename mpl::if_<
            mpl::or_<mpl::equal_to<precision_type, mpl::int_<0> >,
            mpl::greater<precision_type, mpl::int_<113> > >,
            mpl::int_<0>,
            typename mpl::if_<
               mpl::less_equal<precision_type, mpl::int_<53> >,
               mpl::int_<53>,
               typename mpl::if_<
                  mpl::less_equal<precision_type, mpl::int_<64> >,
                  mpl::int_<64>,
                  mpl::int_<113>
               >::type
            >::type
         >::type tag_type;

         result = igamma_temme_large(a, x, pol, static_cast<tag_type const*>(0));
         if(x >= a)
            invert = !invert;
         if(p_derivative)
            *p_derivative = regularised_gamma_prefix(a, x, pol, lanczos_type());
         break;
      }
         }

   if(normalised && (result > 1))
      result = 1;
   if(invert)
   {
      T gam = normalised ? 1 : boost::math::tgamma(a, pol);
      result = gam - result;
   }
   if(p_derivative)
   {
      //
      // Need to convert prefix term to derivative:
      //
      if((x < 1) && (tools::max_value<T>() * x < *p_derivative))
      {
         // overflow, just return an arbitrarily large value:
         *p_derivative = tools::max_value<T>() / 2;
      }

      *p_derivative /= x;
   }

   return result;
}

//
// Ratios of two gamma functions:
//
template <class T, class Policy, class L>
T tgamma_delta_ratio_imp_lanczos(T z, T delta, const Policy& pol, const L&)
{
   BOOST_MATH_STD_USING
   T zgh = z + L::g() - constants::half<T>();
   T result;
   if(fabs(delta) < 10)
   {
      result = exp((constants::half<T>() - z) * boost::math::log1p(delta / zgh, pol));
   }
   else
   {
      result = pow(zgh / (zgh + delta), z - constants::half<T>());
   }
   result *= pow(constants::e<T>() / (zgh + delta), delta);
   result *= L::lanczos_sum(z) / L::lanczos_sum(z + delta);
   return result;
}
//
// And again without Lanczos support this time:
//
template <class T, class Policy>
T tgamma_delta_ratio_imp_lanczos(T z, T delta, const Policy& pol, const lanczos::undefined_lanczos&)
{
   BOOST_MATH_STD_USING
   //
   // The upper gamma fraction is *very* slow for z < 6, actually it's very
   // slow to converge everywhere but recursing until z > 6 gets rid of the
   // worst of it's behaviour.
   //
   T prefix = 1;
   T zd = z + delta;
   while((zd < 6) && (z < 6))
   {
      prefix /= z;
      prefix *= zd;
      z += 1;
      zd += 1;
   }
   if(delta < 10)
   {
      prefix *= exp(-z * boost::math::log1p(delta / z, pol));
   }
   else
   {
      prefix *= pow(z / zd, z);
   }
   prefix *= pow(constants::e<T>() / zd, delta);
   T sum = detail::lower_gamma_series(z, z, pol) / z;
   sum += detail::upper_gamma_fraction(z, z, ::boost::math::policies::get_epsilon<T, Policy>());
   T sumd = detail::lower_gamma_series(zd, zd, pol) / zd;
   sumd += detail::upper_gamma_fraction(zd, zd, ::boost::math::policies::get_epsilon<T, Policy>());
   sum /= sumd;
   if(fabs(tools::max_value<T>() / prefix) < fabs(sum))
      return policies::raise_overflow_error<T>("boost::math::tgamma_delta_ratio<%1%>(%1%, %1%)", "Result of tgamma is too large to represent.", pol);
   return sum * prefix;
}

template <class T, class Policy>
T tgamma_delta_ratio_imp(T z, T delta, const Policy& pol)
{
   BOOST_MATH_STD_USING

   if(z <= 0)
      policies::raise_domain_error<T>("boost::math::tgamma_delta_ratio<%1%>(%1%, %1%)", "Gamma function ratios only implemented for positive arguments (got a=%1%).", z, pol);
   if(z+delta <= 0)
      policies::raise_domain_error<T>("boost::math::tgamma_delta_ratio<%1%>(%1%, %1%)", "Gamma function ratios only implemented for positive arguments (got b=%1%).", z+delta, pol);

   if(floor(delta) == delta)
   {
      if(floor(z) == z)
      {
         //
         // Both z and delta are integers, see if we can just use table lookup
         // of the factorials to get the result:
         //
         if((z <= max_factorial<T>::value) && (z + delta <= max_factorial<T>::value))
         {
            return unchecked_factorial<T>((unsigned)itrunc(z, pol) - 1) / unchecked_factorial<T>((unsigned)itrunc(T(z + delta), pol) - 1);
         }
      }
      if(fabs(delta) < 20)
      {
         //
         // delta is a small integer, we can use a finite product:
         //
         if(delta == 0)
            return 1;
         if(delta < 0)
         {
            z -= 1;
            T result = z;
            while(0 != (delta += 1))
            {
               z -= 1;
               result *= z;
            }
            return result;
         }
         else
         {
            T result = 1 / z;
            while(0 != (delta -= 1))
            {
               z += 1;
               result /= z;
            }
            return result;
         }
      }
   }
   typedef typename lanczos::lanczos<T, Policy>::type lanczos_type;
   return tgamma_delta_ratio_imp_lanczos(z, delta, pol, lanczos_type());
}

template <class T, class Policy>
T gamma_p_derivative_imp(T a, T x, const Policy& pol)
{
   //
   // Usual error checks first:
   //
   if(a <= 0)
      policies::raise_domain_error<T>("boost::math::gamma_p_derivative<%1%>(%1%, %1%)", "Argument a to the incomplete gamma function must be greater than zero (got a=%1%).", a, pol);
   if(x < 0)
      policies::raise_domain_error<T>("boost::math::gamma_p_derivative<%1%>(%1%, %1%)", "Argument x to the incomplete gamma function must be >= 0 (got x=%1%).", x, pol);
   //
   // Now special cases:
   //
   if(x == 0)
   {
      return (a > 1) ? 0 :
         (a == 1) ? 1 : policies::raise_overflow_error<T>("boost::math::gamma_p_derivative<%1%>(%1%, %1%)", 0, pol);
   }
   //
   // Normal case:
   //
   typedef typename lanczos::lanczos<T, Policy>::type lanczos_type;
   T f1 = detail::regularised_gamma_prefix(a, x, pol, lanczos_type());
   if((x < 1) && (tools::max_value<T>() * x < f1))
   {
      // overflow:
      return policies::raise_overflow_error<T>("boost::math::gamma_p_derivative<%1%>(%1%, %1%)", 0, pol);
   }

   f1 /= x;

   return f1;
}

template <class T, class Policy>
inline typename tools::promote_args<T>::type 
   tgamma(T z, const Policy& /* pol */, const mpl::true_)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename lanczos::lanczos<value_type, Policy>::type evaluation_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, forwarding_policy>(detail::gamma_imp(static_cast<value_type>(z), forwarding_policy(), evaluation_type()), "boost::math::tgamma<%1%>(%1%)");
}

template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type
   tgamma(T1 a, T2 z, const Policy&, const mpl::false_)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename lanczos::lanczos<value_type, Policy>::type evaluation_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, forwarding_policy>(
      detail::gamma_incomplete_imp(static_cast<value_type>(a),
      static_cast<value_type>(z), false, true,
      forwarding_policy(), static_cast<value_type*>(0)), "boost::math::tgamma<%1%>(%1%, %1%)");
}

template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type
   tgamma(T1 a, T2 z, const mpl::false_ tag)
{
   return tgamma(a, z, policies::policy<>(), tag);
}

} // namespace detail

template <class T>
inline typename tools::promote_args<T>::type 
   tgamma(T z)
{
   return tgamma(z, policies::policy<>());
}

template <class T, class Policy>
inline typename tools::promote_args<T>::type 
   lgamma(T z, int* sign, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename lanczos::lanczos<value_type, Policy>::type evaluation_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, forwarding_policy>(detail::lgamma_imp(static_cast<value_type>(z), forwarding_policy(), evaluation_type(), sign), "boost::math::lgamma<%1%>(%1%)");
}

template <class T>
inline typename tools::promote_args<T>::type 
   lgamma(T z, int* sign)
{
   return lgamma(z, sign, policies::policy<>());
}

template <class T, class Policy>
inline typename tools::promote_args<T>::type 
   lgamma(T x, const Policy& pol)
{
   return ::boost::math::lgamma(x, 0, pol);
}

template <class T>
inline typename tools::promote_args<T>::type 
   lgamma(T x)
{
   return ::boost::math::lgamma(x, 0, policies::policy<>());
}

template <class T, class Policy>
inline typename tools::promote_args<T>::type 
   tgamma1pm1(T z, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename lanczos::lanczos<value_type, Policy>::type evaluation_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<typename remove_cv<result_type>::type, forwarding_policy>(detail::tgammap1m1_imp(static_cast<value_type>(z), forwarding_policy(), evaluation_type()), "boost::math::tgamma1pm1<%!%>(%1%)");
}

template <class T>
inline typename tools::promote_args<T>::type 
   tgamma1pm1(T z)
{
   return tgamma1pm1(z, policies::policy<>());
}

//
// Full upper incomplete gamma:
//
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type
   tgamma(T1 a, T2 z)
{
   //
   // Type T2 could be a policy object, or a value, select the 
   // right overload based on T2:
   //
   typedef typename policies::is_policy<T2>::type maybe_policy;
   return detail::tgamma(a, z, maybe_policy());
}
template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type
   tgamma(T1 a, T2 z, const Policy& pol)
{
   return detail::tgamma(a, z, pol, mpl::false_());
}
//
// Full lower incomplete gamma:
//
template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type
   tgamma_lower(T1 a, T2 z, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename lanczos::lanczos<value_type, Policy>::type evaluation_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, forwarding_policy>(
      detail::gamma_incomplete_imp(static_cast<value_type>(a),
      static_cast<value_type>(z), false, false,
      forwarding_policy(), static_cast<value_type*>(0)), "tgamma_lower<%1%>(%1%, %1%)");
}
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type
   tgamma_lower(T1 a, T2 z)
{
   return tgamma_lower(a, z, policies::policy<>());
}
//
// Regularised upper incomplete gamma:
//
template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type
   gamma_q(T1 a, T2 z, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename lanczos::lanczos<value_type, Policy>::type evaluation_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, forwarding_policy>(
      detail::gamma_incomplete_imp(static_cast<value_type>(a),
      static_cast<value_type>(z), true, true,
      forwarding_policy(), static_cast<value_type*>(0)), "gamma_q<%1%>(%1%, %1%)");
}
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type
   gamma_q(T1 a, T2 z)
{
   return gamma_q(a, z, policies::policy<>());
}
//
// Regularised lower incomplete gamma:
//
template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type
   gamma_p(T1 a, T2 z, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename lanczos::lanczos<value_type, Policy>::type evaluation_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, forwarding_policy>(
      detail::gamma_incomplete_imp(static_cast<value_type>(a),
      static_cast<value_type>(z), true, false,
      forwarding_policy(), static_cast<value_type*>(0)), "gamma_p<%1%>(%1%, %1%)");
}
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type
   gamma_p(T1 a, T2 z)
{
   return gamma_p(a, z, policies::policy<>());
}

// ratios of gamma functions:
template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type 
   tgamma_delta_ratio(T1 z, T2 delta, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, forwarding_policy>(detail::tgamma_delta_ratio_imp(static_cast<value_type>(z), static_cast<value_type>(delta), forwarding_policy()), "boost::math::tgamma_delta_ratio<%1%>(%1%, %1%)");
}
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type 
   tgamma_delta_ratio(T1 z, T2 delta)
{
   return tgamma_delta_ratio(z, delta, policies::policy<>());
}
template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type 
   tgamma_ratio(T1 a, T2 b, const Policy&)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, forwarding_policy>(detail::tgamma_delta_ratio_imp(static_cast<value_type>(a), static_cast<value_type>(static_cast<value_type>(b) - static_cast<value_type>(a)), forwarding_policy()), "boost::math::tgamma_delta_ratio<%1%>(%1%, %1%)");
}
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type 
   tgamma_ratio(T1 a, T2 b)
{
   return tgamma_ratio(a, b, policies::policy<>());
}

template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type 
   gamma_p_derivative(T1 a, T2 x, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, forwarding_policy>(detail::gamma_p_derivative_imp(static_cast<value_type>(a), static_cast<value_type>(x), forwarding_policy()), "boost::math::gamma_p_derivative<%1%>(%1%, %1%)");
}
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type 
   gamma_p_derivative(T1 a, T2 x)
{
   return gamma_p_derivative(a, x, policies::policy<>());
}

} // namespace math
} // namespace boost

#ifdef BOOST_MSVC
# pragma warning(pop)
#endif

#include <boost/math/special_functions/detail/igamma_inverse.hpp>
#include <boost/math/special_functions/detail/gamma_inva.hpp>
#include <boost/math/special_functions/erf.hpp>

#endif // BOOST_MATH_SF_GAMMA_HPP
