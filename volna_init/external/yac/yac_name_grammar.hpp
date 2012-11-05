#ifndef YAC_NAME_GRAMMAR_HPP
#define YAC_NAME_GRAMMAR_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include <boost/spirit/core.hpp>
#include <boost/spirit/symbols/symbols.hpp>

namespace spirit = boost::spirit;

namespace yac {

////////////////////////////////////////////////////////////////////////////
//
//  name_grammar.
//
////////////////////////////////////////////////////////////////////////////
class name_grammar : public spirit::grammar<name_grammar>
{
    spirit::symbols<char> dummy_reserved_keywords;

public:
    // We are only interested in the names themselves.
    // char is just the smallest space waster available.
    static spirit::symbols<char> reserved_keywords;
    spirit::symbols<char> const & keywords;

    name_grammar(bool check_reserved_keywords = true)
        : keywords(check_reserved_keywords ?
                   reserved_keywords :
                   dummy_reserved_keywords)
    {}

    template <typename ScannerT>
    struct definition {
        definition(name_grammar const & self)
        {
            using spirit::alpha_p;
            using spirit::alnum_p;
            using spirit::lexeme_d;

            name
                = lexeme_d
                  [
                      ((alpha_p | '_') >> *(alnum_p | '_'))
                  ]
                - self.keywords
                ;
        }

        typedef ScannerT scanner_t;
        typedef typename spirit::rule<scanner_t> rule_t;

        rule_t const & start() const { return name; }

    private:
        rule_t name;
    };
};

spirit::symbols<char> name_grammar::reserved_keywords;

} // namespace yac

#endif // YAC_NAME_GRAMMAR_HPP
