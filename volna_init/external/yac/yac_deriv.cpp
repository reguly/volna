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

} // namespace anon


int main()
{
    evaluate1();
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


//  A derivative parser using the (excellent!) YAC framework. 
     
void evaluate1()
{
    using namespace yac;
    
    /* Generate a stack for the following statement list:

       deriv(sin(x)) = cos(x)  
       deriv(cos(x)) = -sin(x) 
       print deriv(sin(2)) 
       print deriv(cos(5)) 
       
    */
    
    virtual_machine vm;

    vm.funcs.add
        ("add##2",  make_function("add#",  &add))
        ("subt##2", make_function("subt#", &subt))
        ("div##2",  make_function("div#",  &div))
        ("mult##2", make_function("mult#", &mult))
        ("sin#1",   make_function("sin",   &sin))
        ("cos#1",   make_function("cos",   &cos)); 

    // Add deriv to the symbol table.
    string const mangled_deriv_name = name_mangler("deriv", 1);
    char const * const deriv_cstr = mangled_deriv_name.c_str();
    
    string const mangled_sin_name = name_mangler("sin", 1);
    char const * const sin_cstr = mangled_sin_name.c_str();
    
    vm.funcs.add
        (deriv_cstr, make_function("deriv", 1))
        ("sin#1",   make_function("sin", &sin))
        ("cos#1",   make_function("cos", &cos)); 
        
// Args and stack 
    vector<string> args;
    shared_ptr<spirit::symbols<double> > local_vars(new spirit::symbols<double>);
    stack stk;

    args.push_back("x");
    local_vars->add("x");

//  deriv(sin(x))
//  stk.push_back(new variable_node(x_ptr));
    stk.push_back(new variable_node(find(*local_vars, "x")));
    stk.push_back(new sys_function_node(*find(vm.funcs, "sin#1")));
    stk.push_back(new sys_function_node(*find(vm.funcs, "deriv#1")));
    stk.push_back(new sys_function_node(*find(vm.funcs, "cos#1")));
    

    // Put it all together
    vm.stk.push_back(new func_def_node(
                         *find(vm.funcs, deriv_cstr),
                         user_func_expression(args, local_vars, stk)));

    // print deriv(sin(2))
    vm.stk.push_back(new const_value_node(2));
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, sin_cstr)));    
    vm.stk.push_back(new sys_function_node(*find(vm.funcs, deriv_cstr)));
    vm.stk.push_back(new print_node);

    // output
    std::cout << "Result of\n"
              << "\t deriv(sin(2)) \n"
              << "\t print deriv(sin(2)) \n"
              << "is: \n";
    evaluate(vm.stk);
    std::cout << std::endl;
}

} // namespace anon
