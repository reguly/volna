#ifndef YAC_STACK_CLOSURE_HPP
#define YAC_STACK_CLOSURE_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include "yac_virtual_machine.hpp"

#include <boost/spirit/core.hpp>
#include <boost/spirit/attribute/closure.hpp>

namespace spirit = boost::spirit;

namespace yac {

struct stack_closure : spirit::closure<stack_closure, stack>
{
    member1 stk;
};

} // namespace yac

#endif // YAC_STACK_CLOSURE_HPP
