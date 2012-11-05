/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#include <boost/spirit/core.hpp>

#include <iostream>
#include <string>

///////////////////////////////////////////////////////////////////////////////
//
// Spirit implementation of the grammar described in
// ../doc/grammar_overview.html
//
// The code does nothing more than flag whether the input data conforms to
// the grammar.
//
///////////////////////////////////////////////////////////////////////////////

namespace spirit = boost::spirit;

using std::string;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

struct skip_grammar : spirit::grammar<skip_grammar>
{
    template <typename ScannerT>
    struct definition {
        definition(skip_grammar const &)
        {
            using spirit::space_p;
            using spirit::anychar_p;

            skip
                = whitespace
                | '\\' >> *whitespace >> '\n'
                | '#' >> *(anychar_p - '\n');

            whitespace
                = space_p - '\n';
        }

        typedef typename spirit::rule<ScannerT> rule_t;

        rule_t const & start() const { return skip; }

    private:
        rule_t skip, whitespace;
    };
};


struct expression_grammar : public spirit::grammar<expression_grammar>
{
    template <typename ScannerT>
    struct definition
    {
        definition(expression_grammar const &)
        {
            using spirit::alpha_p;
            using spirit::alnum_p;
            using spirit::ch_p;
            using spirit::real_p;
            using spirit::str_p;

            expr
                = logical_expr >> !conditional_expr_helper;
            conditional_expr_helper
                = ch_p('?') >> expr >> ':' >> expr >> !conditional_expr_helper;
            logical_expr
                = bitwise_expr >> *((str_p("&&") | "||") >> bitwise_expr);
            bitwise_expr
                = equality_expr >> *((ch_p('&') | '|' | '^') >> equality_expr);
            equality_expr
                = compare_expr >> *((str_p("==") | "!=") >> compare_expr);
            compare_expr
                = shift_expr >> *((ch_p('<') | '>' | "<=" | ">=") >> shift_expr);
            shift_expr
                = add_expr >> *((str_p("<<") | ">>") >> add_expr);
            add_expr
                = mult_expr >> *((ch_p('+') | '-') >> mult_expr);
            mult_expr
                = expr_atom >> *((ch_p('*') | '/' | '%') >> expr_atom);
            expr_atom
                = function
                | number
                | ch_p('(') >> expr >> ')'
                | (ch_p('+') | '-' | '!') >> expr_atom;
            number
                = real_p
                | local_vars
                | global_vars;
            local_vars
                = name;
            global_vars
                = name;
            function
                = name >> '(' >> !arg_list >> ')';
            arg_list
                = expr >> *(ch_p(',') >> expr);
            name
                = (alpha_p | '_') >> *(alnum_p | '_');
        }

        typedef typename spirit::rule<ScannerT> rule_t;

        rule_t const & start() const { return expr; }

    private:
        rule_t add_expr, arg_list, bitwise_expr, compare_expr,
            conditional_expr_helper, equality_expr, expr, expr_atom,
            function, global_vars, local_vars, logical_expr, mult_expr,
            name, number, shift_expr;
    };
};


struct yac_grammar : public spirit::grammar<yac_grammar>
{
    template <typename ScannerT>
    struct definition
    {
        definition(yac_grammar const &)
        {
            using spirit::alpha_p;
            using spirit::alnum_p;
            using spirit::ch_p;

            statement_list
                = *(statement >> +(ch_p('\n') | ';'));
            statement
                = var_assignment
                | function_definition
                | "print " >> expression_list;
            expression_list
                = expression >> *(',' >> expression);
            var_assignment
                = name >> '=' >> expression;
            function_definition
                = name >> '(' >> !name_list >> ')' >> '=' >> expression;
            name_list
                = name >> *(',' >> name);
            name
                = (alpha_p | '_') >> *(alnum_p | '_');
        }

        typedef typename spirit::rule<ScannerT> rule_t;

        rule_t const & start() const { return statement_list; }

    private:
        rule_t expression_list, function_definition, name, name_list,
            statement, statement_list, var_assignment;
        expression_grammar expression;
    };
};


namespace {

char const * prompt()
{
    return "yac> ";
}

} // namespace anon


int main()
{
    typedef spirit::parse_info<string::const_iterator> parse_info_t;

    yac_grammar calculator;
    skip_grammar skipper;

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
            = spirit::parse(buffer.begin(), buffer.end(), *spirit::space_p);

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
        info = spirit::parse(rbuffer.begin(), rbuffer.end(), *spirit::space_p);
        string::size_type const last = rbuffer.end() - info.stop - 1;

        data += buffer;
        if (buffer[last] == '\\')
            continue;

        // The data has been generated, so now get on and parse it properly.
        info = spirit::parse(data.begin(), data.end(), calculator, skipper);

        if (info.full) {
            cout << "-------------------------\n"
                 << "Parsing succeeded\n"
                 << "-------------------------\n";
        } else {
            cerr << "-------------------------\n"
                 << "Parsing failed\n"
                 << "-------------------------\n";
        }

        data.clear();
        cout << prompt();
    }

    return 0;
}

