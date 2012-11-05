#ifndef YAC_BINARY_OP_PARSER_HPP
#define YAC_BINARY_OP_PARSER_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include "yac_lazy_functions.hpp"
#include "yac_virtual_machine.hpp"

#include <boost/spirit/core.hpp>
#include <boost/spirit/phoenix.hpp>
#include <boost/spirit/utility/functor_parser.hpp>

namespace spirit = boost::spirit;

namespace yac {

template <typename TermT, typename OpT>
struct binary_op_parser {

    binary_op_parser(TermT const & term_, OpT const & op_table_)
        : term(term_), op_table(op_table_)
    {}

    typedef stack result_t;

    template <typename ScannerT>
    std::ptrdiff_t
    operator()(ScannerT const & scan, result_t & result) const
    {
        using phoenix::arg1;
        using phoenix::new_;
        using phoenix::var;

        boost::shared_ptr<function> func_store;

        typedef
            typename spirit::match_result<ScannerT, result_t>::type
            match_t;

        match_t match
            =  (term
                [
                    var(result) = arg1
                ]
            >>  *( op_table
                   [
                       var(func_store) = arg1
                   ]
                   >>
                   term
                   [
                       var(result) += arg1,
                       push_back(var(result),
                                 new_<sys_function_node>(var(func_store)))
                   ]
                 )
              ).parse(scan);

        return match.length();
    }

private:
    typename spirit::as_parser<TermT>::type::embed_t term;
    typename spirit::as_parser<OpT>::type::embed_t op_table;
};


//////////////////////////////////
struct binary_op_parser_gen
{
    template <typename TermT, typename OpT>
    spirit::functor_parser<
	binary_op_parser<
            typename spirit::as_parser<TermT>::type,
            typename spirit::as_parser<OpT>::type
        >
    >
    operator()(TermT const & term_, OpT const & op_table_) const
    {
        typedef typename spirit::as_parser<TermT>::type term_t;
        typedef typename spirit::as_parser<OpT>::type op_t;
        typedef binary_op_parser<term_t, op_t> functor_t;
        typedef spirit::functor_parser<functor_t> return_t;

        return return_t(
            functor_t(
                spirit::as_parser<TermT>::convert(term_),
                spirit::as_parser<OpT>::convert(op_table_)
            )
        );
    }
};

//////////////////////////////////
const binary_op_parser_gen binary_op_p = binary_op_parser_gen();

} // namespace yac

#endif // NOT YAC_BINARY_OP_PARSER_HPP
