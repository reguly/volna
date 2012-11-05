#ifndef YAC_EXPRESSION_GRAMMAR_HPP
#define YAC_EXPRESSION_GRAMMAR_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include "yac_binary_op_parser.hpp"
#include "yac_lazy_functions.hpp"
#include "yac_name_grammar.hpp"
#include "yac_stack_closure.hpp"
#include "yac_virtual_machine.hpp"

#include <boost/spirit/core.hpp>
#include <boost/spirit/phoenix.hpp>
#include <boost/spirit/utility/lists.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/dynamic/if.hpp>

#include <string>

namespace spirit = boost::spirit;

namespace yac {

boost::shared_ptr<function>
checked_find(spirit::symbols<boost::shared_ptr<function> > const & symbols,
             std::string const & name)
{
    boost::shared_ptr<function> * ptr = find(symbols, name.c_str());
    BOOST_ASSERT(ptr && ptr->get());
    return *ptr;
}

struct checked_find_impl {
    typedef boost::shared_ptr<function> function_ptr_t;
    typedef spirit::symbols<function_ptr_t> function_table_t;

    template <typename T1, typename T2>
    struct result {
        typedef function_ptr_t type;
    };

    function_ptr_t operator()(function_table_t const & symbols,
                              std::string const & name) const
    {
        function_ptr_t * ptr = find(symbols, name.c_str());
        return ptr ? *ptr : function_ptr_t();
    }
};

phoenix::function<checked_find_impl> const checked_find_ = checked_find_impl();


struct func_expr_closure
    : spirit::closure<func_expr_closure,
                      stack, std::string, int, boost::shared_ptr<function> >
{
    member1 stk;
    member2 name;
    member3 arity;
    member4 func_ptr;
};


struct logical_closure
    : spirit::closure<logical_closure, stack, bool>
{
    member1 stk;
    member2 or_op;
};


struct conditional_closure
    : spirit::closure<conditional_closure, stack, stack, stack>
{
    member1 stk;
    member2 stk1;
    member3 stk2;
};


struct expression_grammar
    : spirit::grammar<expression_grammar, stack_closure::context_t>
{
    typedef boost::shared_ptr<function> function_ptr_t;
    typedef spirit::symbols<function_ptr_t> function_table_t;
    typedef spirit::symbols<double> var_table_t;

    var_table_t const dummy_local_vars;

    function_table_t const & functions;
    var_table_t const & global_vars;
    var_table_t const & local_vars;

    expression_grammar(function_table_t const & funcs,
                       var_table_t const & gvars)
        : functions(funcs), global_vars(gvars), local_vars(dummy_local_vars)
    {}

    expression_grammar(function_table_t const & funcs,
                       var_table_t const & gvars,
                       var_table_t const & lvars)
        : functions(funcs), global_vars(gvars), local_vars(lvars)
    {}

    template <typename ScannerT>
    struct definition {
        definition(expression_grammar const & self)
            : name(false)
        {
            using phoenix::arg1;
            using phoenix::arg2;
            using phoenix::construct_;
            using phoenix::if_;
            using phoenix::new_;
            using phoenix::var;

            using spirit::ch_p;
            using spirit::epsilon_p;
            using spirit::if_p;
            using spirit::list_p;
            using spirit::nothing_p;
            using spirit::real_p;
            using spirit::str_p;

            boolean_op.add
                ("&&", false)
                ("||", true);

            add_op.add
                ("+", checked_find(self.functions, "add##2"))
                ("-", checked_find(self.functions, "subtract##2"));

            bitwise_op.add
                ("&", checked_find(self.functions, "bitwise_and##2"))
                ("|", checked_find(self.functions, "bitwise_or##2"))
                ("^", checked_find(self.functions, "bitwise_xor##2"));

            compare_op.add
                ("<",  checked_find(self.functions, "less##2"))
                (">",  checked_find(self.functions, "greater##2"))
                ("<=", checked_find(self.functions, "less_equal##2"))
                (">=", checked_find(self.functions, "greater_equal##2"));

            equality_op.add
                ("==", checked_find(self.functions, "equal##2"))
                ("!=", checked_find(self.functions, "not_equal##2"));

            mult_op.add
                ("*", checked_find(self.functions, "multiply##2"))
                ("/", checked_find(self.functions, "divide##2"))
                ("%", checked_find(self.functions, "mod##2"));

            shift_op.add
                ("<<", checked_find(self.functions, "shift_left##2"))
                (">>", checked_find(self.functions, "shift_right##2"));

            unary_op.add
                ("+", function_ptr_t())
                ("-", checked_find(self.functions, "negate##1"))
                ("!", checked_find(self.functions, "logical_not##1"));

            top = expr[self.stk = arg1];

            expr
                =  logical_expr[expr.stk = arg1]
                >> !conditional_expr_helper[expr.stk += arg1];

            conditional_expr_helper
                =  (ch_p('?')
                >>  expr
                    [
                        conditional_expr_helper.stk1 = arg1
                    ]
                >>  ':'
                >>  expr
                    [
                        conditional_expr_helper.stk2 = arg1
                    ]
                >>  !conditional_expr_helper
                    [
                        conditional_expr_helper.stk2 += arg1
                    ]
                   )
                   [
                       push_back(conditional_expr_helper.stk,
                                 new_<branch_node>(
                                     conditional_expr_helper.stk1,
                                     conditional_expr_helper.stk2))
                   ]
                ;

           logical_expr
                =  bitwise_expr
                   [
                       logical_expr.stk = arg1
                   ]
                >> *(logical_expr_helper
                     [
                         logical_expr.stk += arg1
                     ]
                    )
                ;

            logical_expr_helper
                =  boolean_op
                   [
                       logical_expr_helper.or_op = arg1
                   ]
                >> bitwise_expr
                   [
                       if_(logical_expr_helper.or_op)
                       [
                           push_back(logical_expr_helper.stk,
                                     new_<or_node>(arg1))
                       ]
                       .else_
                       [
                           logical_expr_helper.stk = arg1,
                           push_back(logical_expr_helper.stk,
                                     new_<sys_function_node>(
                                         checked_find(self.functions, "logical_and##2")))
                       ]
                   ]
                ;

            bitwise_expr
                = binary_op_p(equality_expr, bitwise_op)[bitwise_expr.stk = arg1];

            equality_expr
                = binary_op_p(compare_expr, equality_op)[equality_expr.stk = arg1];

            compare_expr
                = binary_op_p(shift_expr, compare_op)[compare_expr.stk = arg1];

            shift_expr
                = binary_op_p(add_expr, shift_op)[shift_expr.stk = arg1];

            add_expr
                = binary_op_p(mult_expr, add_op)[ add_expr.stk = arg1 ];

            mult_expr
                = binary_op_p(expr_atom, mult_op)[ mult_expr.stk = arg1 ];

            expr_atom
                = number
                  [
                      expr_atom.stk = arg1
                  ]
                | func
                  [
                      expr_atom.stk = arg1
                  ]
                | ('(' >>
                   expr
                   [
                       expr_atom.stk = arg1
                   ]
                   >> ')')
                | unary_expr
                  [
                      expr_atom.stk = arg1
                  ]
                ;

            unary_expr
                =  unary_op
                   [
                      unary_expr.func_ptr = arg1
                   ]
                >> expr_atom
                   [
                       unary_expr.stk = arg1,
                       if_(unary_expr.func_ptr != function_ptr_t())
                       [
                           push_back(unary_expr.stk,
                                     new_<sys_function_node>(unary_expr.func_ptr))
                       ]
                   ]
                ;

            number
                = real_p
                  [
                      push_back(number.stk,
                                new_<const_value_node>(arg1))
                  ]
                | self.local_vars
                  [
                      push_back(number.stk,
                                new_<variable_node>(address_of(arg1)))
                  ]
                | self.global_vars
                  [
                      push_back(number.stk,
                                new_<variable_node>(address_of(arg1)))
                  ]
	      | str_p("x")
	      [ 
	       push_back(number.stk, new_<x_node>())
		]
	      | str_p("y")
	      [ 
	       push_back(number.stk, new_<y_node>())
		]
	      | str_p("t")
	      [
	       push_back(number.stk, new_<t_node>())
	       ]
                ;

            func
                =  name
                   [
                       func.arity = 0,
                       func.name = construct_<std::string>(arg1, arg2)
                   ]
                >> ('(' >> !list_p(arg, ',') >> ')')
                   [
                       func.func_ptr
                           = checked_find_(var(self.functions),
                                           mangled_name(func.name, func.arity))
                   ]
                   // If the function is found in the symbol table, then
                   // add it to the stack. Else fail the parser.
                >> if_p(func.func_ptr != function_ptr_t())
                   [
                       epsilon_p
                       [
                           push_back(func.stk,
                                     new_<sys_function_node>(func.func_ptr))
                       ]
                   ]
                  .else_p[nothing_p]
                ;

            arg
                = expr
                  [
                      func.arity += 1,
                      func.stk += arg1
                  ]
                ;
        }

        typedef ScannerT scanner_t;
        typedef spirit::rule<scanner_t> rule_t;

        rule_t const & start() const { return top; }

    private:
        typedef spirit::rule<scanner_t, stack_closure::context_t> srule_t;
        typedef spirit::rule<scanner_t, func_expr_closure::context_t> frule_t;
        typedef spirit::rule<scanner_t, logical_closure::context_t> lrule_t;
        typedef spirit::rule<scanner_t, conditional_closure::context_t> crule_t;

        rule_t arg, top;
        srule_t add_expr, and_expr, bitwise_expr, compare_expr,
            equality_expr, expr, expr_atom, logical_expr,
            number, or_expr, or_op, mult_expr, shift_expr;
        frule_t unary_expr, func;
        crule_t conditional_expr_helper;
        lrule_t logical_expr_helper;
        name_grammar name;

        spirit::symbols<bool> boolean_op;
        function_table_t and_op, add_op, bitwise_op, compare_op, equality_op,
            shift_op, mult_op, unary_op;
    };
};

} // namespace yac

#endif // YAC_EXPRESSION_GRAMMAR_HPP
