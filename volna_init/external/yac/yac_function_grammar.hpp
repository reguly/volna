#ifndef YAC_FUNCTION_GRAMMAR_HPP
#define YAC_FUNCTION_GRAMMAR_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include "yac_expression_grammar.hpp"
#include "yac_lazy_functions.hpp"
#include "yac_name_grammar.hpp"
#include "yac_stack_closure.hpp"
#include "yac_virtual_machine.hpp"

#include <boost/spirit/core.hpp>
#include <boost/spirit/dynamic/lazy.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/utility/lists.hpp>
#include <boost/spirit/attribute/closure.hpp>
#include <boost/spirit/phoenix.hpp>

#include <string>
#include <vector>

namespace spirit = boost::spirit;

namespace yac {

struct func_closure
    : spirit::closure<func_closure,
                      std::string,
                      std::vector<std::string>,
                      boost::shared_ptr<spirit::symbols<double> >,
                      std::string,
                      boost::shared_ptr<expression_grammar> >
{
    member1 name;
    member2 args;
    member3 local_vars;
    member4 mangled_name;
    member5 expr;
};


struct function_grammar
    : public spirit::grammar<function_grammar, stack_closure::context_t>
{
    typedef spirit::symbols<boost::shared_ptr<function> > function_table_t;
    typedef spirit::symbols<double> var_table_t;

    function_table_t & functions;
    var_table_t const & global_vars;

    function_grammar(function_table_t & f, var_table_t const & v)
        : functions(f), global_vars(v)
    {}

    template <typename ScannerT>
    struct definition {
        definition(function_grammar const & self)
        {
            using phoenix::arg1;
            using phoenix::arg2;
            using phoenix::construct_;
            using phoenix::new_;
            using phoenix::var;

            using spirit::ch_p;
            using spirit::lazy_p;
            using spirit::list_p;

            top = func_def;

            typedef boost::shared_ptr<spirit::symbols<double> >
                shared_symbols_dbl_ptr;

            func_def
                =  // Parse "foo(a,b)". If successful, place the function
                   // in the self.functions symbol table.
                   func_decl
                   [
                       func_def.mangled_name = mangled_name(func_def.name,
                                                            size(func_def.args)),

                       add_symbol(var(self.functions),
                                  func_def.mangled_name,
                                  construct_<user_function>(
                                      func_def.name,
                                      size(func_def.args))),

                       // Create the expression_grammar instance that will
                       // parse the function definition.
                       reset(func_def.expr,
                             new_<expression_grammar>(
                                 var(self.functions),
                                 var(self.global_vars),
                                 *func_def.local_vars))
                   ]
                >> '='
                >> lazy_p(*func_def.expr)
                   [
                       // Add a node to the stack which will reset the
                       // definition of the named function.

                       // The func_def_node constructor:
                       // func_def_node(shared_ptr<function>,
                       //               user_func_expression const &);
                       push_back(self.stk,
                                 new_<func_def_node>(
                                     *find_symbol(var(self.functions),
                                                  func_def.mangled_name),
                                     construct_<user_func_expression>(
                                         func_def.args,
                                         func_def.local_vars,
                                         arg1)))
                   ]
                ;

            // This rule parses "foo(a,b)", disallowing "foo(a,a)".
            // Note that the parsed data is used to fill the "func_def"
            // closure variables, not those of "func_decl".
            func_decl
                =  name
                   [
                       func_def.name = construct_<std::string>(arg1, arg2)
                   ]
                >> ch_p('(')
                   [
                       reset(func_def.local_vars, new_<var_table_t>())
                   ]
                >> !list_p((name - lazy_p(*func_def.local_vars))
                           [
                               push_back(func_def.args,
                                         construct_<std::string>(arg1, arg2)),
                               add_symbol(*func_def.local_vars,
                                          construct_<std::string>(arg1, arg2))
                           ]
                         , ',')
                >> ')';
        }

        typedef ScannerT scanner_t;
        typedef spirit::rule<scanner_t> rule_t;

        rule_t const & start() const { return top; }

    private:
        typedef
        typename spirit::rule<scanner_t, func_closure::context_t> func_rule_t;

        rule_t top, func_decl;
        func_rule_t func_def;

        name_grammar name;
    };
};

} // namespace yac

#endif // YAC_FUNCTION_GRAMMAR_HPP
