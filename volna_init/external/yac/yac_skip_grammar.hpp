#ifndef YAC_SKIP_GRAMMAR_HPP
#define YAC_SKIP_GRAMMAR_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include <boost/spirit/core.hpp>

namespace spirit = boost::spirit;

namespace yac {

////////////////////////////////////////////////////////////////////////////
//
//  The skip grammar
//
////////////////////////////////////////////////////////////////////////////
struct skip_grammar : spirit::grammar<skip_grammar>
{
    template <typename ScannerT>
    struct definition {
        definition(skip_grammar const &)
        {
            using spirit::space_p;
            using spirit::anychar_p;
	    using spirit::blank_p;

            skip
                = whitespace
	      //| '\\' >> *whitespace >> '\n'
                | '#' >> *(anychar_p - '\n');

            whitespace
	      //= space_p - '\n';
	      = blank_p;
        }

        typedef typename spirit::rule<ScannerT> rule_t;

        rule_t const & start() const { return skip; }

    private:
        rule_t skip, whitespace;
    };
};

} // namespace yac

#endif // YAC_SKIP_GRAMMAR_HPP
