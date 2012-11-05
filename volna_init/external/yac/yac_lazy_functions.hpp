#ifndef YAC_LAZY_FUNCTIONS_HPP
#define YAC_LAZY_FUNCTIONS_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include "yac_virtual_machine.hpp"

#include <boost/spirit/symbols/symbols.hpp>

#include <boost/spirit/phoenix/functions.hpp>
#include <boost/shared_ptr.hpp>

#include <string>

namespace yac {

////////////////////////////////////////////////////////////////////////////
//
//  A lazy function wrapper for spirit::symbols::add.
//
////////////////////////////////////////////////////////////////////////////
struct add_symbol_impl {
    typedef boost::shared_ptr<function> function_ptr;
    typedef spirit::symbols<function_ptr> function_table_t;

    template <typename T, typename Arg1, typename Arg2 = spirit::nil_t>
    struct result
    {
        typedef void type;
    };

    template <typename SymbolT>
    void operator()(SymbolT & symbols,
                    std::string const & name) const
    {
        symbols.add(name.c_str());
    }

    void operator()(function_table_t & symbols,
                    std::string const & name,
                    user_function const & func) const
    {
        function_ptr * fp = find(symbols, name.c_str());

        if (fp)
            fp->reset(new user_function(func));
        else
            symbols.add(name.c_str(), function_ptr(new user_function(func)));
    }
};

phoenix::function<add_symbol_impl> const add_symbol = add_symbol_impl();


////////////////////////////////////////////////////////////////////////////
//
//  A lazy function wrapper for spirit::symbols::find.
//
////////////////////////////////////////////////////////////////////////////
struct find_symbol_impl {
    template <typename T, typename Arg>
    struct result
    {
        typedef typename T::symbol_data_t * type;
    };

    template <typename SymbolT>
    typename result<SymbolT, std::string>::type
    operator()(SymbolT const & symbols, std::string const & sym) const
    {
        return find(symbols, sym.c_str());
    }
};

phoenix::function<find_symbol_impl> const find_symbol = find_symbol_impl();


////////////////////////////////////////////////////////////////////////////
//
//  A lazy function wrapper for both yac::stack::push_back and for
//  the generic ContainerT::push_back.
//
////////////////////////////////////////////////////////////////////////////
struct push_back_impl {
    template <typename T, typename Arg>
    struct result
    {
        typedef void type;
    };

    void operator()(stack & stk, node * n) const
    {
        BOOST_ASSERT(n);
        stk.push_back(n);
    }

    template <typename ContainerT>
    void operator()(ContainerT & c,
                    typename ContainerT::value_type const & data) const
    {
        c.push_back(data);
    }
};

phoenix::function<push_back_impl> const push_back = push_back_impl();


////////////////////////////////////////////////////////////////////////////
//
//  A lazy function wrapper for ContainerT::size().
//
////////////////////////////////////////////////////////////////////////////
struct size_impl {

    template <typename ContainerT>
    struct result {
        typedef std::size_t type;
    };

    template <typename ContainerT>
    std::size_t operator()(ContainerT const & c) const
    {
        return c.size();
    }
};

phoenix::function<size_impl> const size = size_impl();


////////////////////////////////////////////////////////////////////////////
//
//  A lazy function wrapper for
//  std::string const name_mangler(std::string const &, std::size_t const).
//
////////////////////////////////////////////////////////////////////////////
struct mangled_name_impl {

    template <typename NameT, typename ArityT>
    struct result {
        typedef std::string const type;
    };

    std::string const operator()(std::string const & name,
                                 std::size_t const arity) const
    {
        return name_mangler(name, arity);
    }
};

phoenix::function<mangled_name_impl> const mangled_name = mangled_name_impl();


////////////////////////////////////////////////////////////////////////////
//
//  A lazy function wrapper for boost::shared_ptr<T>::reset.
//
////////////////////////////////////////////////////////////////////////////
struct reset_impl {
    template <typename T1, typename T2>
    struct result {
        typedef void type;
    };

    template <typename T>
    void operator()(boost::shared_ptr<T> & sh_ptr, T * ptr = 0) const
    {
        sh_ptr.reset(ptr);
    }
};

phoenix::function<reset_impl> const reset = reset_impl();


////////////////////////////////////////////////////////////////////////////
//
//  A lazy function wrapper for boost::reference_wrapper<double>::get_pointer.
//
////////////////////////////////////////////////////////////////////////////
struct address_of_impl {
    template <typename T>
    struct result {
        typedef double * type;
    };

    double * const operator()(boost::reference_wrapper<double> const & ref) const
    {
        return ref.get_pointer();
    }
};

phoenix::function<address_of_impl> const address_of = address_of_impl();

} // namespace yac

#endif // NOT YAC_LAZY_FUNCTIONS_HPP
