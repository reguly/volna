#ifndef YAC_VARIABLE_GRAMMAR_HPP
#define YAC_VARIABLE_GRAMMAR_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include "yac_expression_grammar.hpp"
#include "yac_name_grammar.hpp"
#include "yac_stack_closure.hpp"

#include <boost/spirit/core.hpp>
#include <boost/spirit/phoenix.hpp>
#include <boost/spirit/symbols/symbols.hpp>

namespace spirit = boost::spirit;

namespace yac {

struct var_closure
    : spirit::closure<var_closure, std::string, stack>
{
    member1 name;
    member2 stk;
};


struct variable_grammar
    : spirit::grammar<variable_grammar, stack_closure::context_t>
{
    typedef spirit::symbols<boost::shared_ptr<function> > function_table_t;
    typedef spirit::symbols<double> var_table_t;

    function_table_t const & functions;
    var_table_t & vars;

    variable_grammar(function_table_t const & f, var_table_t & v)
        : functions(f), vars(v)
    {}

    template <typename ScannerT>
    struct definition
    {
        definition(variable_grammar const & self)
            : expression(self.functions, self.vars)
        {
            using phoenix::arg1;
            using phoenix::arg2;
            using phoenix::construct_;
            using phoenix::if_;
            using phoenix::new_;
            using phoenix::var;

            top
                = step2;

            // Parse a variable assignment statement such as "x=2".
            // If successful, place the variable in the self.vars symbol table.
            // Also generates a stack representation of "x=2".
            step2
                = step1
                  [
                      if_(find_symbol(var(self.vars), step2.name) == (double*)0)
                      [
                          add_symbol(var(self.vars), step2.name)
                      ]
                  ]
                  [
                      push_back(self.stk,
                                new_<variable_node>(find_symbol(var(self.vars),
                                                                step2.name))),
                      self.stk += step2.stk,
                      push_back(self.stk, new_<assign_node>())
                  ]
                ;

            // This rule parses "x=2".
            // Note that the grammar is not recursive, so we
            // can fill the "step2" closure variables directly.
            step1
                =  name
                   [
                       step2.name = construct_<std::string>(arg1, arg2)
                   ]
                >> '='
                >> expression
                   [
                       step2.stk = arg1
                   ]
                ;
        }

        typedef ScannerT scanner_t;
        typedef typename spirit::rule<scanner_t> rule_t;

        rule_t const & start() const { return top; }

    private:
        typedef
        typename spirit::rule<scanner_t, var_closure::context_t> var_rule_t;

        rule_t top, step1;
        var_rule_t step2;

        name_grammar name;
        expression_grammar expression;
    };
};

} // namespace yac

#endif // YAC_VARIABLE_GRAMMAR_HPP
