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

#include "yac_gnuplot_grammar.hpp"
#include "yac_skip_grammar.hpp"
#include "yac_virtual_machine.hpp"

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
public:
  controller(std::vector<double> X)
    : x_vector(X)
    {
        char const * const keywords[] = {
	  "x",
            "return",
            "pi", "e",
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
            ("pi", 3.1415926536)
            ("e",  2.7182818285);

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

    template <typename ItT>
    std::vector<double> parse(ItT first, ItT last)
    {
        using phoenix::arg1;
        using phoenix::var;

        typedef spirit::parse_info<ItT> parse_info_t;

        gnuplot_grammar calculator(vm.funcs, vm.global_vars);

        parse_info_t info = spirit::parse(first, last,
                                          calculator[var(vm.stk) = arg1],
                                          skip_grammar());

        if (info.full) {
	  //evaluate(vm.stk, x_vector.at(0));
	  std::vector<double> result = vector_eval(vm.stk, x_vector);
	  return result;
        } else {
            cerr << "-------------------------\n"
                 << "Parsing failed\n"
                 << "-------------------------\n";
	    return std::vector<double> (x_vector.size(), 0.0);
        }
    }
private:
    virtual_machine vm;
  std::vector<double> x_vector;
};

} // namespace yac

////////////////////////////////////////////////////////////////////////////
//
//  What follows below is glue to handle data input in different ways,
//  from file, as piped input from stdin or interactive input.
//
////////////////////////////////////////////////////////////////////////////
namespace {

char const * prompt()
{
    return "yac> ";
}


int usage()
{
    cerr << "Usage: yac [Input mode]\n"
         << "\tInput mode\n"
         << "\t-i, --interactive    Data input interactively\n"
         << "\t-f, --file FILE      Data input from FILE\n"
         << "\t-p, --piped          Data piped through stdin\n";
    return 0;
}


enum input_mode {
    file_mode,
    piped_mode,
    interactive_mode,
    unrecognized_mode
};


input_mode get_input_mode(string const & mode_arg)
{
    if (mode_arg == "-f" || mode_arg == "--file")
        return file_mode;
    if (mode_arg == "-p" || mode_arg == "--piped")
        return piped_mode;
    if (mode_arg == "-i" || mode_arg == "--interactive")
        return interactive_mode;
    return unrecognized_mode;
}

} // namespace anon


int main(int argc, char * argv[])
{
    if (argc < 2)
        return usage();

    input_mode const mode = get_input_mode(argv[1]);
    if (mode == unrecognized_mode)
        return usage();

    if (mode == file_mode) {
        if (argc != 3)
            return usage();
    } else if (argc != 2)
        return usage();

    unsigned int n = 20;
    std::vector<double> X(n+1, 0.0);
    for (unsigned int i=0; i<=n; ++i) {
      X.at(i) = (double)i/n;
    }
    yac::controller calc(X);

    switch (mode) {
    case file_mode:
    {
        char const * const file = argv[2];
        spirit::file_iterator<> file_it(file);

        if (!file_it) {
                cerr << "Unable to open file '" << file << "'!"
                     << endl;
                return 1;
        }

	std::vector<double> result = calc.parse(file_it, file_it.make_end());
	for (unsigned int i=0; i<result.size();++i) {
	  std::cerr << result.at(i) << std::endl;
	}
        break;
    }
    case piped_mode:
    {
        string data;
        std::ostringstream ss;
        while (getline(cin, data))
            ss << data << '\n';
        data = ss.str();

        calc.parse(data.begin(), data.end());
        break;
    }
    case interactive_mode:
    {
        typedef spirit::parse_info<string::const_iterator> parse_info_t;
        using spirit::space_p;

        cout << "/////////////////////////////////////////////////////////\n\n"
             << "\t\tYet Another Calculator...\n\n"
             << "/////////////////////////////////////////////////////////\n\n"
             << "Type an expression...or [q or Q] to quit\n\n";

        string data;
        string buffer;
        cout << prompt();
        while (getline(cin, buffer)) {
            // getline strips the final '\n', so add it back.
            buffer += '\n';

            // Quit?
            // Find the first non-whitespace character.
            parse_info_t info
                = spirit::parse(buffer.begin(), buffer.end(), *space_p);

            if (info.full) {
                // Nothing but whitespace
                if (data.empty())
                    cout << prompt();
                continue;
            } else {
                string::size_type const first = info.stop - buffer.begin();
                if (buffer[first] == 'q' || buffer[first] == 'Q')
                    break;
            }

            // Continuation?
            // Find the last non-whitespace character.
            string const rbuffer = string(buffer.rbegin(), buffer.rend());
            // Cannot fail because we're testing the same data as above.
            info = spirit::parse(rbuffer.begin(), rbuffer.end(), *space_p);
            string::size_type const last = rbuffer.end() - info.stop - 1;

            data += buffer;
            if (buffer[last] == '\\')
                continue;

            // The data has been generated, so now get on and parse it properly.
            calc.parse(data.begin(), data.end());
            data.clear();
            cout << prompt();
        }
        break;
    }
    case unrecognized_mode:
        break;
    }

    return 0;
}
