#ifndef YAC_GNUPLOT_GRAMMAR_HPP
#define YAC_GNUPLOT_GRAMMAR_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include "yac_expression_grammar.hpp"
#include "yac_function_grammar.hpp"
#include "yac_variable_grammar.hpp"
#include "yac_virtual_machine.hpp"

#include <boost/spirit/core.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/utility/lists.hpp>
#include <boost/spirit/attribute/closure.hpp>
#include <boost/spirit/phoenix.hpp>

namespace spirit = boost::spirit;

namespace yac {

struct gnuplot_grammar
    : spirit::grammar<gnuplot_grammar, stack_closure::context_t>
{
    typedef spirit::symbols<boost::shared_ptr<function> > function_table_t;
    typedef spirit::symbols<double> var_table_t;

    function_table_t & functions;
    var_table_t & vars;

    gnuplot_grammar(function_table_t & f, var_table_t & v)
        : functions(f), vars(v)
    {}

    template <typename ScannerT>
    struct definition {
        definition(gnuplot_grammar const & self)
            : var_assign(self.functions, self.vars),
              function_definition(self.functions, self.vars),
              expression(self.functions, self.vars)
        {
            using phoenix::arg1;
            using phoenix::new_;

            using spirit::ch_p;
            using spirit::list_p;

            statement_list
                = *(statement >> +(ch_p('\n') | ';'));

            statement
                = var_assign[self.stk += arg1]
                | function_definition[self.stk += arg1]
                | "return " >> list_p(expr, ',')
                ;

            expr
                = expression
                  [
                      self.stk += arg1,
                      push_back(self.stk, new_<print_node>())
                  ]
                ;
        }

        typedef ScannerT scanner_t;
        typedef spirit::rule<scanner_t> rule_t;

        rule_t const & start() const { return statement_list; }

    private:
        typedef spirit::rule<scanner_t, stack_closure::context_t> stk_rule_t;

        rule_t expr, statement_list, statement;

        variable_grammar var_assign;
        function_grammar function_definition;
        expression_grammar expression;
    };
};

} // namespace yac

#endif // YAC_GNUPLOT_GRAMMAR_HPP
