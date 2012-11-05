/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#define YAC_SINGLE_COMPILATION_UNIT
#include "yac_virtual_machine.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::string;
using std::vector;

using boost::shared_ptr;

namespace {

void evaluate1();
void evaluate2();
void evaluate3();
void evaluate4();

} // namespace anon


int main()
{
    evaluate1();
    evaluate2();
    evaluate3();
    evaluate4();
    return 0;
}


namespace {

double add(double a, double b) { return a + b; }
double subt(double a, double b) { return a - b; }
double mult(double a, double b) { return a * b; }
double div(double a, double b) { return a / b; }
double less(double a, double b) { return a < b; }


typedef double(* func2_ptr_t)(double, double);
typedef double(* func1_ptr_t)(double);

shared_ptr<yac::function>
make_function(char const * name, func2_ptr_t func_ptr)
{
    using yac::function;
    using yac::function2;
    return shared_ptr<function>(new function2<double>(name, func_ptr));
}


shared_ptr<yac::function>
make_function(char const * name, func1_ptr_t func_ptr)
{
    using yac::function;
    using yac::function1;
    return shared_ptr<function>(new function1<double>(name, func_ptr));
}


shared_ptr<yac::function>
make_function(char const * name, int arity)
{
    using yac::function;
    using yac::user_function;
    return shared_ptr<function>(new user_function(name, arity));
}


void evaluate1()
{
    using namespace yac;
    
    /* Generate a stack for the following statement list:

       2+3*sin(sqrt(4)/5)

       Demonstrates how the virtual machine evaluates numerical expressions.
    */
    virtual_machine vm;

    vm.funcs.add
        ("add##2",  make_function("add#",  &add))
        ("div##2",  make_function("div#",  &div))
        ("mult##2", make_function("mult#", &mult))
        ("sin#1",   make_function("sin",   &sin))
        ("sqrt#1",  make_function("sqrt",  &sqrt));

    vm.stk.push_back(new const_value_node(2));
    vm.stk.push_back(new const_value_node(3));
    vm.stk.push_back(new const_value_node(4));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, "sqrt#1")));
    vm.stk.push_back(new const_value_node(5));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, "div##2")));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, "sin#1")));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, "mult##2")));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, "add##2")));

    double const result = evaluate(vm.stk);
    std::cout << "Result of (2+3*sin(sqrt(4)/5)) is "
              << result
              << '\n' << std::endl;
}


void evaluate2()
{
    using namespace yac;
    
    /* Generate a stack for the following statement list:

       x=3
       y=2*x
       print y

       Demonstrates how the virtual machine adds new variables to the
       symbol table and how it then assigns a value to these variables
       by evaluation of the stack.
    */
    virtual_machine vm;

    vm.funcs.add
        ("mult##2", make_function("mult#", &mult));

    // x=3
    vm.global_vars.add("x");

    vm.stk.push_back(new variable_node(find(vm.global_vars, "x")));
    vm.stk.push_back(new const_value_node(3));
    vm.stk.push_back(new assign_node);

    // y=2*x
    vm.global_vars.add("y");

    vm.stk.push_back(new variable_node(find(vm.global_vars, "y")));
    vm.stk.push_back(new const_value_node(2));
    vm.stk.push_back(new variable_node(find(vm.global_vars, "x")));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, "mult##2")));
    vm.stk.push_back(new assign_node);

    // print y
    vm.stk.push_back(new variable_node(find(vm.global_vars, "y")));
    vm.stk.push_back(new print_node);

    // output
    std::cout << "Result of\n"
              << "\tx=3\n"
              << "\ty=2*x\n"
              << "\tprint y\n"
              << "is:\n";
    evaluate(vm.stk);
    std::cout << std::endl;
}


void evaluate3()
{
    using namespace yac;
    
    /* Generate a stack for the following statement list:

       foo(a)=2*a
       print foo(2)
       foo(a,b)=a*b
       print foo(2,3)
       foo(a,b)=a+b
       print foo(2,3)

       Demonstrates how the virtual machine adds new functions to the
       symbol table. These functions are added at the appropriate
       point of evaluation of the statement list, just as the parser would do.

       Also demonstrates how it then assigns a value to these functions
       by evaluation of the stack.

       The test shows how name mangling is used to overload functions
       with the same name but different arity.

       Also shows that dynamic redefinition of the function body is
       handled correctly.
    */
    virtual_machine vm;

    vm.funcs.add
        ("add##2",  make_function("add#",  &add))
        ("mult##2", make_function("mult#", &mult));

    string const mangled_foo1_name = name_mangler("foo", 1);
    string const mangled_foo2_name = name_mangler("foo", 2);
    char const * const foo1_cstr = mangled_foo1_name.c_str();
    char const * const foo2_cstr = mangled_foo2_name.c_str();

    // foo(a)=2*a
    vm.funcs.add
        (foo1_cstr, make_function("foo", 1));

    vector<string> args;
    shared_ptr<spirit::symbols<double> > local_vars(new spirit::symbols<double>);
    stack stk;

    args.push_back("a");
    local_vars->add("a");

    stk.push_back(new const_value_node(2));
    stk.push_back(new variable_node(find(*local_vars, "a")));
    stk.push_back(new sys_function_node(*find(vm.funcs, "mult##2")));

    vm.stk.push_back(new func_def_node(
                         *find(vm.funcs, foo1_cstr),
                         user_func_expression(args, local_vars, stk)));

    // print foo(2)
    vm.stk.push_back(new const_value_node(2));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, foo1_cstr)));
    vm.stk.push_back(new print_node);

    // foo(a,b)=a*b
    vm.funcs.add
        (foo2_cstr, make_function("foo", 2));

    args.clear();
    local_vars.reset(new spirit::symbols<double>);
    stk.clear();

    args.push_back("a");
    local_vars->add("a");
    args.push_back("b");
    local_vars->add("b");

    stk.push_back(new variable_node(find(*local_vars, "a")));
    stk.push_back(new variable_node(find(*local_vars, "b")));
    stk.push_back(new sys_function_node(*find(vm.funcs, "mult##2")));

    vm.stk.push_back(new func_def_node(
                         *find(vm.funcs, foo2_cstr),
                         user_func_expression(args, local_vars, stk)));

    // print foo(2,3)
    vm.stk.push_back(new const_value_node(2));
    vm.stk.push_back(new const_value_node(3));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, foo2_cstr)));
    vm.stk.push_back(new print_node);

    // foo(a,b)=a+b
    args.clear();
    local_vars.reset(new spirit::symbols<double>);
    stk.clear();

    args.push_back("a");
    local_vars->add("a");
    args.push_back("b");
    local_vars->add("b");

    stk.push_back(new variable_node(find(*local_vars, "a")));
    stk.push_back(new variable_node(find(*local_vars, "b")));
    stk.push_back(new sys_function_node(*find(vm.funcs, "add##2")));

    vm.stk.push_back(new func_def_node(
                         *find(vm.funcs, foo2_cstr),
                         user_func_expression(args, local_vars, stk)));

    // print foo(2,3)
    vm.stk.push_back(new const_value_node(2));
    vm.stk.push_back(new const_value_node(3));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, foo2_cstr)));
    vm.stk.push_back(new print_node);

    // output
    std::cout << "Result of\n"
              << "\tfoo(a)=2*a\n"
              << "\tprint foo(2)\n"
              << "\tfoo(a,b)=a*b\n"
              << "\tprint foo(2,3)\n"
              << "\tfoo(a,b)=a+b\n"
              << "\tprint foo(2,3)\n"
              << "is:\n";
    evaluate(vm.stk);
    std::cout << std::endl;
}


void evaluate4()
{
    using namespace yac;
    
    /* Generate a stack for the following statement list:

       factorial(x) = (x < 0.1) ? 1 : x * factorial(x - 1)
       print factorial(2)
       print factorial(5)

       Demonstrates that recursive user-defined functions are handled correctly.
    */
    virtual_machine vm;

    vm.funcs.add
        ("less##2", make_function("less#", &less))
        ("mult##2", make_function("mult#", &mult))
        ("subt##2", make_function("subt#", &subt));

    // Add factorial to the symbol table.
    string const mangled_factorial_name = name_mangler("factorial", 1);
    char const * const factorial_cstr = mangled_factorial_name.c_str();
    vm.funcs.add
        (factorial_cstr, make_function("factorial", 1));

    // The factorial definition
    vector<string> args;
    shared_ptr<spirit::symbols<double> > local_vars(new spirit::symbols<double>);
    stack stk;

    args.push_back("x");
    local_vars->add("x");

    double * const x_ptr = find(*local_vars, "x");
    BOOST_ASSERT(x_ptr);

    // x < 0.1
    stk.push_back(new variable_node(x_ptr));
    stk.push_back(new const_value_node(0.1));
    stk.push_back(new sys_function_node(*find(vm.funcs, "less##2")));

    // ? 1
    stack stk1;
    stk1.push_back(new const_value_node(1));

    // : x * factorial(x - 1)
    stack stk2;
    stk2.push_back(new variable_node(x_ptr));
    stk2.push_back(new variable_node(x_ptr));
    stk2.push_back(new const_value_node(1));
    stk2.push_back(new sys_function_node(*find(vm.funcs, "subt##2")));
    stk2.push_back(new sys_function_node(*find(vm.funcs, factorial_cstr)));
    stk2.push_back(new sys_function_node(*find(vm.funcs, "mult##2")));

    // Put it all together
    stk.push_back(new branch_node(stk1, stk2));

    vm.stk.push_back(new func_def_node(
                         *find(vm.funcs, factorial_cstr),
                         user_func_expression(args, local_vars, stk)));

    // print factorial(2)
    vm.stk.push_back(new const_value_node(2));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, factorial_cstr)));
    vm.stk.push_back(new print_node);

     // print factorial(5)
    vm.stk.push_back(new const_value_node(5));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, factorial_cstr)));
    vm.stk.push_back(new print_node);

    // output
    std::cout << "Result of\n"
              << "\tfactorial(x) = (x < 0.1) ? 1 : x * factorial(x - 1)\n"
              << "\tprint factorial(2)\n"
              << "\tprint factorial(5)\n"
              << "is:\n";
    evaluate(vm.stk);
    std::cout << std::endl;
}

} // namespace anon
