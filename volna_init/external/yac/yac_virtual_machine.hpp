#ifndef YAC_VIRTUAL_MACHINE_HPP
#define YAC_VIRTUAL_MACHINE_HPP
/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/


#include <boost/shared_ptr.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/spirit/symbols/symbols.hpp>

#include <list>
#include <string>
#include <vector>

namespace spirit = boost::spirit;

///////////////////////////////////////////////////////////////////////////////
//
//  The virtual machine
//
//      Please see ../doc/virtual_machine.html for a detailed description
//      of the virtual machine's operation.
//
///////////////////////////////////////////////////////////////////////////////

namespace yac {

//////////////////////////////////
class stack;

struct node
{
    virtual ~node() {}
    virtual boost::shared_ptr<node> clone() const = 0;
    virtual std::size_t nbranches() const { return 0; }
    virtual stack * branch(std::size_t) { return 0; }
};


//////////////////////////////////
class stack
{
    typedef boost::shared_ptr<node> node_ptr;

public:
    typedef std::list<node_ptr> container_t;
    typedef container_t::size_type size_type;
    typedef container_t::difference_type difference_type;
    typedef boost::indirect_iterator<container_t::iterator> iterator;
    typedef boost::indirect_iterator<container_t::const_iterator> const_iterator;

    stack() {}
    stack(stack const & other);

    stack & operator=(stack const & other);
    stack & operator+=(stack const & other);

    void push_back(node * new_node);
    size_type size() const;
    bool empty() const;
    void clear();

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    iterator erase(iterator where);
    iterator erase(iterator begin, iterator end);

    iterator insert(iterator where, node * new_node);

private:
    container_t data;
};

//////////////////////////////////
struct function;

struct virtual_machine
{
    stack stk;
    spirit::symbols<double> global_vars;
    spirit::symbols<boost::shared_ptr<function> > funcs;
};


  RealType evaluate(stack stk, RealType &x, RealType &y, RealType &t);

  std::vector<RealType> vectorizedEval( stack stk, 
					std::vector<RealType> &X,
					std::vector<RealType> &Y, 
					RealType &t );

std::string const name_mangler( std::string const & name,
				std::size_t const arity );


//////////////////////////////////
struct user_function;
struct user_func_expression;

struct function
{
    function();
    function(std::string const &);
    virtual ~function();

    virtual user_function * as_user_function();

    std::string const & name() const;

    typedef std::size_t size_type;
    virtual size_type arity() const = 0;

  RealType eval(std::vector<double> const & params, RealType x, RealType y,
  RealType t) const;

private:
  virtual RealType eval_priv(std::vector<double> const & params, RealType
			     x, RealType y, RealType t) const = 0;
    std::string name_;
};


///////////////////////////////////////////////////////////////////////////////
//
// Concrete implementations of function types.
//
///////////////////////////////////////////////////////////////////////////////
template <typename T>
struct function1 : public function
{
    typedef T value_type;
    typedef value_type(* func_ptr_type)(value_type);

    function1(std::string const & name, func_ptr_type func)
        : function(name),
          func_(func)
    {}

    size_type arity() const { return 1; }

private:
  RealType eval_priv(std::vector<double> const & params, RealType x,
		     RealType y, RealType t) const
    {
        return func_(static_cast<value_type>(params[0]));
    }

    func_ptr_type func_;
};


//////////////////////////////////
template <typename T>
struct function2 : public function
{
    typedef T value_type;
    typedef value_type(* func_ptr_type)(value_type, value_type);

    function2(std::string const & name, func_ptr_type func)
        : function(name),
          func_(func)
    {}

    size_type arity() const { return 2; }

private:
  RealType eval_priv(std::vector<double> const & params, RealType x,
		     RealType y, RealType t) const
    {
        return func_(static_cast<value_type>(params[0]),
                     static_cast<value_type>(params[1]));
    }

    func_ptr_type func_;
};


//////////////////////////////////
struct user_function : public function
{
    typedef spirit::symbols<double> symbol_table_t;
    typedef boost::shared_ptr<symbol_table_t> symbol_table_ptr;

    user_function();
    user_function(std::string const & name, size_type const function_arity);
    user_function(user_function const & other);

    user_function & operator=(user_func_expression const & ufe);

    user_function * as_user_function();
    spirit::symbols<double> const & local_vars() const;
    size_type arity() const;

private:
    user_function & operator=(user_function const &);

  RealType eval_priv(std::vector<double> const & params, RealType x,
		     RealType y, RealType t) const;

    size_type arity_;
    std::vector<std::string> arg_names_;
    symbol_table_t local_vars_;
    stack stk_;
};


///////////////////////////////////////////////////////////////////////////////
//
// Implementations of various node types.
//
///////////////////////////////////////////////////////////////////////////////
struct number_node : public node
{
    virtual double value() const = 0;
};


//////////////////////////////////
struct function_node : public node
{
    virtual function const & func() const = 0;
};


//////////////////////////////////
struct const_value_node : public number_node
{
    const_value_node(double val);
    boost::shared_ptr<node> clone() const;
    double value() const;
private:
    double val_;
};


//////////////////////////////////
struct variable_node : public number_node
{
    variable_node(double * ptr);
    boost::shared_ptr<node> clone() const;
    double value() const;
    double & var() const;
    double * ptr_;
};


//////////////////////////////////
struct sys_function_node : public function_node
{
    sys_function_node(boost::shared_ptr<function> const & func_ptr);
    boost::shared_ptr<node> clone() const;
    function const & func() const;
private:
    boost::shared_ptr<function> func_ptr_;
};


//////////////////////////////////
struct print_node : public node
{
    boost::shared_ptr<node> clone() const;
};


//////////////////////////////////
struct assign_node : public node
{
    boost::shared_ptr<node> clone() const;
    void assign(double & var, double val) const;
};


//////////////////////////////////
struct user_func_expression
{
    typedef boost::shared_ptr<spirit::symbols<double> > symbol_table_ptr;

    user_func_expression(std::vector<std::string> const & an,
                         symbol_table_ptr const & lv,
                         stack const & s)
        : arg_names(an), local_vars(lv), stk(s)
    {}

    std::vector<std::string> arg_names;
    symbol_table_ptr local_vars;
    stack stk;
};


//////////////////////////////////
struct func_def_node : public node
{
    func_def_node(boost::shared_ptr<function> const & func_ptr,
                  user_func_expression const & val);

    boost::shared_ptr<node> clone() const;

    void assign() const;
private:
    boost::shared_ptr<function> func_ptr_;
    user_func_expression val_;
};


//////////////////////////////////
struct or_node : public function_node
{
    or_node(stack const &);

    boost::shared_ptr<node> clone() const;

    function const & func() const;

    std::size_t nbranches() const;
    stack * branch(std::size_t);

private:

    struct or_func : public function {
        or_func(stack const &);

        size_type arity() const;
      RealType eval_priv(std::vector<double> const & params, RealType x,
			 RealType y, RealType t) const;

        stack rhs_stk_;
    };

    or_func func_;
};




//////////////////////////////////
struct branch_node : public function_node
{
    branch_node(stack const & stk1, stack const & stk2);

    boost::shared_ptr<node> clone() const;

    function const & func() const;

    std::size_t nbranches() const;
    stack * branch(std::size_t i);

private:

    struct branch_func : public function {
        branch_func(stack const & stk1, stack const & stk2);

        size_type arity() const;
      RealType eval_priv(std::vector<double> const & params, RealType x,
			 RealType y, RealType t) const;

        stack stk1_;
        stack stk2_;
    };

    branch_func func_;
};

//////////////////////////////////
struct x_node : public node
{
  x_node();
  boost::shared_ptr<node> clone() const;

};

//////////////////////////////////
struct y_node : public node
{
  y_node();
  boost::shared_ptr<node> clone() const;

};

//////////////////////////////////
struct t_node : public node
{
  t_node();
  boost::shared_ptr<node> clone() const;
  
};

} // namespace yac

#ifdef YAC_SINGLE_COMPILATION_UNIT
#include "yac_virtual_machine.cpp"
#endif

#endif // NOT YAC_VIRTUAL_MACHINE_HPP
