#ifndef MATHPARSER_HPP
#define MATHPARSER_HPP

/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#if defined PHOENIX_LIMIT && PHOENIX_LIMIT < 6
#undef PHOENIX_LIMIT
#endif

#ifndef PHOENIX_LIMIT
#define PHOENIX_LIMIT 6
#endif

#define YAC_SINGLE_COMPILATION_UNIT

#include "external/yac/yac_gnuplot_grammar.hpp"
#include "external/yac/yac_skip_grammar.hpp"
#include "external/yac/yac_virtual_machine.hpp"

#include <boost/spirit/iterator/file_iterator.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>

#include <iostream>
#include <sstream>
#include <string>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::string;

namespace spirit = boost::spirit;

namespace {

int    bitwise_and(int a, int b) { return a & b; }
int    bitwise_or(int a, int b) { return a | b; }
int    bitwise_xor(int a, int b) { return a ^ b; }
bool   logical_and(bool a, bool b) { return a && b; }
double logical_not(double a) { return !a; }
double equal(double a, double b) { return a == b; }
double not_equal(double a, double b) { return a != b; }
double less(double a, double b) { return a < b; }
double greater(double a, double b) { return a > b; }
double less_equal(double a, double b) { return a <= b; }
double greater_equal(double a, double b) { return a >= b; }
int    shift_left(int a, int b) { return a << b; }
int    shift_right(int a, int b) { return a >> b; }
double add(double a, double b) { return a + b; }
double subtract(double a, double b) { return a - b; }
double multiply(double a, double b) { return a * b; }
double divide(double a, double b) { return a / b; }
int    mod(int a, int b) { return a % b; }
double negate(double a) { return -a; }
double sgn(double a) { return a >= 0.0 ? 1.0 : -1.0; }

typedef int(* ifunc2_ptr_t)(int, int);
typedef bool(* bfunc2_ptr_t)(bool, bool);
typedef double(* dfunc2_ptr_t)(double, double);
typedef double(* dfunc1_ptr_t)(double);

} // namespace anon

namespace yac {

boost::shared_ptr<function> make_func(string const & name, dfunc1_ptr_t ptr)
{
    return boost::shared_ptr<function>(new function1<double>(name, ptr));
}


boost::shared_ptr<function> make_func(string const & name, dfunc2_ptr_t ptr)
{
    return boost::shared_ptr<function>(new function2<double>(name, ptr));
}

boost::shared_ptr<function> make_func(string const & name, bfunc2_ptr_t ptr)
{
    return boost::shared_ptr<function>(new function2<bool>(name, ptr));
}

boost::shared_ptr<function> make_func(string const & name, ifunc2_ptr_t ptr)
{
    return boost::shared_ptr<function>(new function2<int>(name, ptr));
}

class controller {
private:
	  virtual_machine vm;
	  std::vector<RealType> x_vector, y_vector;
public:
  RealType time;
  controller():
    time(0.)
    {
      std::vector<RealType> X;
      std::vector<RealType> Y;
      X.push_back(0.);
      Y.push_back(0.);
      x_vector = X;
      y_vector = Y;

      char const * const keywords[] = {
	"x", "y",
	"return",
	"pi",
	"abs", "acos", "acosh", "asin", "asinh", "atan",
	"atanh", "besj0", "besj1", "besy0", "besy1", "ceil",
	"cos", "cosh",
#if defined(__GNUC__)
	"erf", "erfc",
#endif
	"exp", "floor", "log", "log10", "sgn", "sin", "sinh",
	"sqrt", "tan", "tanh", "atan2",
      };
      std::size_t const keywords_size =
	sizeof(keywords) / sizeof(keywords[0]);

      for (std::size_t i = 0; i != keywords_size; ++i)
	name_grammar::reserved_keywords.add(keywords[i]);

      vm.global_vars.add
	("pi", 3.1415926536);

      vm.funcs.add
	("bitwise_and##2",   make_func("bitwise_and#",   &::bitwise_and))
	("bitwise_or##2",    make_func("bitwise_or#",    &::bitwise_or))
	("bitwise_xor##2",   make_func("bitwise_xor#",   &::bitwise_xor))
	("logical_and##2",   make_func("logical_and#",   &::logical_and))
	("logical_not##1",   make_func("logical_not#",   &::logical_not))
	("equal##2",         make_func("equal#",         &::equal))
	("not_equal##2",     make_func("not_equal#",     &::not_equal))
	("less##2",          make_func("less#",          &::less))
	("greater##2",       make_func("greater#",       &::greater))
	("less_equal##2",    make_func("less_equal#",    &::less_equal))
	("greater_equal##2", make_func("greater_equal#", &::greater_equal))
	("shift_left##2",    make_func("shift_left#",    &::shift_left))
	("shift_right##2",   make_func("shift_right#",   &::shift_right))
	("add##2",           make_func("add#",           &::add))
	("subtract##2",      make_func("subtract#",      &::subtract))
	("multiply##2",      make_func("multiply#",      &::multiply))
	("divide##2",        make_func("divide#",        &::divide))
	("mod##2",           make_func("mod#",           &::mod))
	("negate##1",        make_func("negate#",        &::negate))
	("abs#1",   make_func("abs",   &std::abs))
	("acos#1",  make_func("acos",  &std::acos))
	("asin#1",  make_func("asin",  &std::asin))
	("atan#1",  make_func("atan",  &std::atan))
	("acosh#1", make_func("acosh", &boost::math::acosh))
	("asinh#1", make_func("asinh", &boost::math::asinh))
	("atanh#1", make_func("atanh", &boost::math::atanh))
	("besj0#1", make_func("besj0", &j0))
	("besj1#1", make_func("besj1", &j1))
	("besy0#1", make_func("besy0", &y0))
	("besy1#1", make_func("besy1", &y1))
	("ceil#1",  make_func("ceil",  &std::ceil))
	("cos#1",   make_func("cos",   &std::cos))
	("cosh#1",  make_func("cosh",  &std::cosh))
#if defined(__GNUC__)
	("erf#1",   make_func("erf",   &erf))
	("erfc#1",  make_func("erfc",  &erfc))
#endif
	("exp#1",   make_func("exp",   &std::exp))
	("floor#1", make_func("floor", &std::floor))
	("log#1",   make_func("log",   &std::log))
	("log10#1", make_func("log10", &std::log10))
	("sgn#1",   make_func("sgn",   &sgn))
	("sin#1",   make_func("sin",   &std::sin))
	("sinh#1",  make_func("sinh",  &std::sinh))
	("sqrt#1",  make_func("sqrt",  &std::sqrt))
	("tan#1",   make_func("tan",   &std::tan))
	("tanh#1",  make_func("tanh",  &std::tanh))
	("atan2#2", make_func("atan2", &std::atan2));
    }
  
  controller(std::vector<RealType> &X, std::vector<RealType> &Y, RealType &time_)
    : x_vector(X), y_vector(Y), time(time_)
  {
        char const * const keywords[] = {
	  "x", "y", "t",
	  "return",
	  "pi",
            "abs", "acos", "acosh", "asin", "asinh", "atan",
            "atanh", "besj0", "besj1", "besy0", "besy1", "ceil",
            "cos", "cosh",
#if defined(__GNUC__)
            "erf", "erfc",
#endif
            "exp", "floor", "log", "log10", "sgn", "sin", "sinh",
            "sqrt", "tan", "tanh", "atan2",
        };
        std::size_t const keywords_size =
            sizeof(keywords) / sizeof(keywords[0]);

        for (std::size_t i = 0; i != keywords_size; ++i)
                name_grammar::reserved_keywords.add(keywords[i]);

        vm.global_vars.add
	  ("pi", 3.1415926536);

        vm.funcs.add
            ("bitwise_and##2",   make_func("bitwise_and#",   &::bitwise_and))
            ("bitwise_or##2",    make_func("bitwise_or#",    &::bitwise_or))
            ("bitwise_xor##2",   make_func("bitwise_xor#",   &::bitwise_xor))
            ("logical_and##2",   make_func("logical_and#",   &::logical_and))
            ("logical_not##1",   make_func("logical_not#",   &::logical_not))
            ("equal##2",         make_func("equal#",         &::equal))
            ("not_equal##2",     make_func("not_equal#",     &::not_equal))
            ("less##2",          make_func("less#",          &::less))
            ("greater##2",       make_func("greater#",       &::greater))
            ("less_equal##2",    make_func("less_equal#",    &::less_equal))
            ("greater_equal##2", make_func("greater_equal#", &::greater_equal))
            ("shift_left##2",    make_func("shift_left#",    &::shift_left))
            ("shift_right##2",   make_func("shift_right#",   &::shift_right))
            ("add##2",           make_func("add#",           &::add))
            ("subtract##2",      make_func("subtract#",      &::subtract))
            ("multiply##2",      make_func("multiply#",      &::multiply))
            ("divide##2",        make_func("divide#",        &::divide))
            ("mod##2",           make_func("mod#",           &::mod))
            ("negate##1",        make_func("negate#",        &::negate))
            ("abs#1",   make_func("abs",   &std::abs))
            ("acos#1",  make_func("acos",  &std::acos))
            ("asin#1",  make_func("asin",  &std::asin))
            ("atan#1",  make_func("atan",  &std::atan))
            ("acosh#1", make_func("acosh", &boost::math::acosh))
            ("asinh#1", make_func("asinh", &boost::math::asinh))
            ("atanh#1", make_func("atanh", &boost::math::atanh))
            ("besj0#1", make_func("besj0", &j0))
            ("besj1#1", make_func("besj1", &j1))
            ("besy0#1", make_func("besy0", &y0))
            ("besy1#1", make_func("besy1", &y1))
            ("ceil#1",  make_func("ceil",  &std::ceil))
            ("cos#1",   make_func("cos",   &std::cos))
            ("cosh#1",  make_func("cosh",  &std::cosh))
#if defined(__GNUC__)
            ("erf#1",   make_func("erf",   &erf))
            ("erfc#1",  make_func("erfc",  &erfc))
#endif
            ("exp#1",   make_func("exp",   &std::exp))
            ("floor#1", make_func("floor", &std::floor))
            ("log#1",   make_func("log",   &std::log))
            ("log10#1", make_func("log10", &std::log10))
            ("sgn#1",   make_func("sgn",   &sgn))
            ("sin#1",   make_func("sin",   &std::sin))
            ("sinh#1",  make_func("sinh",  &std::sinh))
            ("sqrt#1",  make_func("sqrt",  &std::sqrt))
            ("tan#1",   make_func("tan",   &std::tan))
            ("tanh#1",  make_func("tanh",  &std::tanh))
            ("atan2#2", make_func("atan2", &std::atan2));
    }

  void updateReservedVariables( std::vector<RealType> &X, 
				std::vector<RealType> &Y )
  {
    x_vector = X;
    y_vector = Y;
  }

  void updateTime( RealType& t ) {
    time = t;
  }

    template <typename ItT>
    std::vector<RealType> parse( ItT first, ItT last )
  {
        using phoenix::arg1;
        using phoenix::var;

        typedef spirit::parse_info<ItT> parse_info_t;

        gnuplot_grammar calculator(vm.funcs, vm.global_vars);

        parse_info_t info = spirit::parse(first, last,
                                          calculator[var(vm.stk) = arg1],
                                          skip_grammar());

        if (info.full) {
	  std::vector<RealType> result = vectorizedEval(vm.stk,
						      x_vector,
						      y_vector, time );
	  return result;
        } else {
            cerr << "-------------------------\n"
                 << "Math parsing failed\n"
                 << "-------------------------\n";
	    return std::vector<RealType> ( x_vector.size(), 0.0 );
        }
    }
};

} // namespace yac

#endif // MATHPARSER_HPP
