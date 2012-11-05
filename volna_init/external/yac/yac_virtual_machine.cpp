/*=============================================================================
    Copyright (c) 2004 Angus Leeming

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#include "yac_virtual_machine.hpp"

#include <boost/assert.hpp>
#include <boost/next_prior.hpp>

#include <iostream>
#include <map>
#include <sstream>

using std::map;
using std::ostringstream;
using std::string;
using std::vector;


namespace yac {

namespace {

void copy(stack::container_t & lhs, stack::container_t const & rhs)
{
    lhs.clear();

    stack::container_t::const_iterator it = rhs.begin();
    stack::container_t::const_iterator const end = rhs.end();
    for(; it != end; ++it)
        lhs.push_back((*it)->clone());
}

} // namespace anon


stack::stack(stack const & other)
{
    copy(data, other.data);
}


stack & stack::operator=(stack const & other)
{
    if (&other != this)
        copy(data, other.data);
    return *this;
}


stack & stack::operator+=(stack const & other)
{
    container_t other_data;
    copy(other_data, other.data);
    data.insert(data.end(), other_data.begin(), other_data.end());
    return *this;
}


void stack::push_back(node * new_node)
{
    BOOST_ASSERT(new_node);
    data.push_back(node_ptr(new_node));
}


stack::size_type stack::size() const
{
    return data.size();
}


bool stack::empty() const
{
    return data.empty();
}


void stack::clear() {
    data.clear();
}


namespace {

stack::iterator make_iterator(stack::container_t::iterator it)
{
    return boost::make_indirect_iterator(it);
}


stack::const_iterator make_iterator(stack::container_t::const_iterator it)
{
    return boost::make_indirect_iterator(it);
}

} // namespace anon


stack::iterator stack::begin()
{
    return make_iterator(data.begin());
}


stack::iterator stack::end()
{
    return make_iterator(data.end());
}


stack::const_iterator stack::begin() const
{
    return make_iterator(data.begin());
}


stack::const_iterator stack::end() const
{
    return make_iterator(data.end());
}


stack::iterator stack::erase(iterator where)
{
    return make_iterator(data.erase( where.base() ));
}


stack::iterator stack::erase(iterator begin, iterator end)
{
    return make_iterator(data.erase( begin.base(), end.base() ));
}


stack::iterator stack::insert(iterator where, node * new_node)
{
    BOOST_ASSERT(new_node);
    return make_iterator(data.insert(where.base(), node_ptr(new_node)));
}


function::function()
{}


function::function(string const & n)
    : name_(n)
{}


function::~function()
{}


user_function * function::as_user_function()
{
    return 0;
}


string const & function::name() const
{
    return name_;
}


  RealType function::eval(vector<double> const & params, RealType x,
			  RealType y, RealType t) const
{
    BOOST_ASSERT(params.size() == arity());
    return eval_priv(params, x, y, t);
}


namespace {

void update_stack(stack & lhs_stk,
                  user_function::symbol_table_t const & lhs_local_vars,
                  user_function::symbol_table_t const & rhs_local_vars,
                  vector<string> const & arg_names);

} // namespace anon

user_function::user_function()
    : arity_(0)
{}


user_function::user_function(string const & name, size_type const function_arity)
        : function(name),
          arity_(function_arity)
    {}


user_function::user_function(user_function const & other)
    : function(other),
      arity_(other.arity_),
      arg_names_(other.arg_names_),
      local_vars_(other.local_vars_),
      stk_(other.stk_)
{
    update_stack(stk_, local_vars_, other.local_vars_, arg_names_);
}


user_function & user_function::operator=(user_func_expression const & ufe)
{
    BOOST_ASSERT(ufe.arg_names.size() == arity_ && ufe.local_vars.get());

    arg_names_ = ufe.arg_names;
    local_vars_ = *ufe.local_vars;
    stk_ = ufe.stk;

    update_stack(stk_, local_vars_, *ufe.local_vars, arg_names_);

    return *this;
}


user_function * user_function::as_user_function()
{
    return this;
}


spirit::symbols<double> const & user_function::local_vars() const
{
    return local_vars_;
}


function::size_type user_function::arity() const
{
    return arity_;
}


  RealType user_function::eval_priv(vector<double> const & params,
				    RealType x, RealType y, RealType t) const
{
    BOOST_ASSERT(arg_names_.size() == arity_);

    // Create a copy of the local variables symbol table.
    symbol_table_t local_vars_copy = local_vars_;

    // Assign values to the elements of this symbol table using the
    // argument list data stored in params.
    std::size_t const size = arg_names_.size();
    for (std::size_t i = 0; i != size; ++i) {
        double * const ptr = find(local_vars_copy, arg_names_[i].c_str());
        BOOST_ASSERT(ptr);
        *ptr = params[i];
    }

    // Ensure that the variable_nodes in stk_copy point to the appropriate
    // entry in local_vars_copy (rather than those of local_vars_).
    stack stk_copy = stk_;
    update_stack(stk_copy, local_vars_copy, local_vars_, arg_names_);

    return evaluate(stk_copy, x, y, t);
}


namespace {

void update_stack(stack & stk, map<double *, double *> const & ptr_map)
{
    typedef map<double *, double *> map_t;

    stack::iterator it = stk.begin();
    stack::iterator const end = stk.end();
    for (; it != end; ++it) {
        variable_node * vn = dynamic_cast<variable_node *>(&*it);
        if (vn) {
            map_t::const_iterator mit = ptr_map.find(vn->ptr_);
            if (mit == ptr_map.end()) {
                std::cerr << "WARNING! Pointer "
                          << vn->ptr_ << ", " << *vn->ptr_
                          << " not found!\n";
                continue;
            } else {
                vn->ptr_ = mit->second;
            }
        }

        std::size_t const nbranches = it->nbranches();
        if (nbranches == 0)
            continue;

        for(std::size_t i = 0; i != nbranches; ++i) {
            stack * bstk = it->branch(i);
            if (!bstk)
                continue;
            update_stack(*bstk, ptr_map);
        }
    }
}

void update_stack(stack & lhs_stk,
                  user_function::symbol_table_t const & lhs_local_vars,
                  user_function::symbol_table_t const & rhs_local_vars,
                  vector<string> const & arg_names)
{
    if (lhs_stk.empty())
        return;

    // Generate a mapping relating the addresses of the data in
    // rhs_local_vars to those in lhs_local_vars.
    typedef map<double *, double *> map_t;
    map_t ptr_map;
    for (std::size_t i = 0; i != arg_names.size(); ++i) {
        double * ptr1 = find(rhs_local_vars, arg_names[i].c_str());
        double * ptr2 = find(lhs_local_vars, arg_names[i].c_str());
        ptr_map[ptr1] = ptr2;
    }

    // Must now loop over the stk_ and change all the pointers that
    // point into rhs_local_vars_ so that they point into local_vars_.
    update_stack(lhs_stk, ptr_map);
}

} // namespace anon


const_value_node::const_value_node(double val)
    : val_(val)
{}


boost::shared_ptr<node> const_value_node::clone() const
{
    return boost::shared_ptr<node>(new const_value_node(*this));
}


double const_value_node::value() const
{
    return val_;
}


variable_node::variable_node(double * ptr)
    : ptr_(ptr)
{
    BOOST_ASSERT(ptr_);
}


boost::shared_ptr<node> variable_node::clone() const
{
    return boost::shared_ptr<node>(new variable_node(*this));
}


double variable_node::value() const
{
    return *ptr_;
}


double & variable_node::var() const
{
    return *ptr_;
}


sys_function_node::sys_function_node(boost::shared_ptr<function> const & func_ptr)
{
    func_ptr_ = func_ptr;
}


boost::shared_ptr<node> sys_function_node::clone() const
{
    return boost::shared_ptr<node>(new sys_function_node(*this));
}


function const & sys_function_node::func() const
{
    return *func_ptr_;
}


boost::shared_ptr<node> print_node::clone() const
{
    return boost::shared_ptr<node>(new print_node(*this));
}


boost::shared_ptr<node> assign_node::clone() const
{
    return boost::shared_ptr<node>(new assign_node(*this));
}


void assign_node::assign(double & var, double val) const
{
    var = val;
}


func_def_node::func_def_node(boost::shared_ptr<function> const & func_ptr,
                             user_func_expression const & val)
    : func_ptr_(func_ptr), val_(val)
{}


boost::shared_ptr<node> func_def_node::clone() const
{
    return boost::shared_ptr<node>(new func_def_node(*this));
}


void func_def_node::assign() const
{
    BOOST_ASSERT(func_ptr_.get());
    user_function * func = func_ptr_->as_user_function();
    BOOST_ASSERT(func);
    *func = val_;
}


or_node::or_node(stack const & stk)
    : func_(stk)
{}


boost::shared_ptr<node> or_node::clone() const
{
    return boost::shared_ptr<node>(new or_node(*this));
}


function const & or_node::func() const
{
    return func_;
}


std::size_t or_node::nbranches() const
{
    return 1;
}


stack * or_node::branch(std::size_t i)
{
    return i == 0 ? &func_.rhs_stk_ : 0;
}


or_node::or_func::or_func(stack const & stk)
    : function("or#"),
      rhs_stk_(stk)
{}


function::size_type or_node::or_func::arity() const
{
    return 1;
}


  RealType or_node::or_func::eval_priv(vector<double> const & params,
				       RealType x, RealType y, RealType t) const
{
  return bool(params[0]) ? true : bool(evaluate(rhs_stk_, x, y, t));
}


branch_node::branch_node(stack const & stk1, stack const & stk2)
    : func_(stk1, stk2)
{}


boost::shared_ptr<node> branch_node::clone() const
{
    return boost::shared_ptr<node>(new branch_node(*this));
}


function const & branch_node::func() const
{
    return func_;
}


std::size_t branch_node::nbranches() const
{
    return 2;
}


stack * branch_node::branch(std::size_t i)
{
    if (i > 1)
        return 0;
    return i == 0 ? &func_.stk1_ : &func_.stk2_;
}


branch_node::branch_func::branch_func(stack const & stk1, stack const & stk2)
    : function("branch#"),
    stk1_(stk1), stk2_(stk2)
{}


function::size_type branch_node::branch_func::arity() const
{
    return 1;
}


  RealType branch_node::branch_func::eval_priv(vector<double> const &
					       params, RealType x, RealType
					       y, RealType t) const
{
  return evaluate(params[0] ? stk1_ : stk2_, x, y, t);
}

  x_node::x_node() {}

  boost::shared_ptr<node> x_node::clone() const
  {
    return boost::shared_ptr<node>(new x_node(*this));
  }

  y_node::y_node() {}

  boost::shared_ptr<node> y_node::clone() const
  {
    return boost::shared_ptr<node>(new y_node(*this));
  }

  t_node::t_node() {}

  boost::shared_ptr<node> t_node::clone() const
  {
    return boost::shared_ptr<node>(new t_node(*this));
  }

string const name_mangler(string const & name,
                          std::size_t const arity)
{
    ostringstream ss;
    ss << name << '#' << arity;
    return ss.str();
}


namespace {

number_node * get_number_node(stack::iterator it)
{
    number_node * nn = dynamic_cast<number_node *>(&*it);
    BOOST_ASSERT(nn);
    return nn;
}


variable_node * get_variable_node(stack::iterator it)
{
    variable_node * vn = dynamic_cast<variable_node *>(&*it);
    BOOST_ASSERT(vn);
    return vn;
}


std::vector<double> const make_arg_list(stack::iterator const & begin,
                                        stack::iterator const & end)
{
    std::vector<double> vec;
    vec.reserve(std::distance(begin, end));

    stack::iterator it = begin;
    for (; it != end; ++it)
        vec.push_back(get_number_node(it)->value());

    return vec;
}

} // namespace anon


  RealType evaluate(stack stk, RealType &x_value, RealType &y_value, 
		    RealType &t_value )
{

  RealType res = 0.0;

    if (stk.empty())
        return res;

    stack::iterator it = stk.begin();
    while (it != stk.end()) {
        if (dynamic_cast<print_node const *>(&*it))
        {
            stack::iterator val_it = boost::prior(it);
            double const result = get_number_node(val_it)->value();
            it = stk.erase(val_it, boost::next(it));

	    res = result;
            continue;
        }

        if (function_node const * f =
            dynamic_cast<function_node const *>(&*it))
        {
            function const & func = f->func();

            // arity *must* be a signed type or boost::prior will fail.
            stack::difference_type const arity = func.arity();
            stack::iterator data_begin = boost::prior(it, arity);
            RealType const result = func.eval(make_arg_list(data_begin,
							    it),
					      x_value, y_value, t_value);

#ifdef DEBUG
            std::vector<double> args = make_arg_list(data_begin, it);
            std::cerr << "function " << func.name() << '(';
            for (std::size_t i = 0; i != args.size(); ++i) {
                if (i != 0)
                    std::cerr << ',';
                std::cerr << args[i];
            }
            std::cerr << ") = " << result << "\n";
#endif // DEBUG

            it = stk.erase(data_begin, boost::next(it));
            stk.insert(it, new yac::const_value_node(result));
            continue;
        }

        if (x_node const * xx =
            dynamic_cast<x_node const *>(&*it))
        {
	  it = stk.erase(it, boost::next(it));
	  stk.insert(it, new yac::const_value_node(x_value));
            continue;
        }

	if (y_node const * yy =
            dynamic_cast<y_node const *>(&*it))
	  {
	    it = stk.erase(it, boost::next(it));
	    stk.insert(it, new yac::const_value_node(y_value));
            continue;
	  }

	if (t_node const * tt =
            dynamic_cast<t_node const *>(&*it))
	  {
	    it = stk.erase(it, boost::next(it));
	    stk.insert(it, new yac::const_value_node(t_value));
            continue;
	  }

	if (assign_node const * a =
            dynamic_cast<assign_node const *>(&*it))
	  {
            stack::iterator val_it = boost::prior(it);
            stack::iterator var_it = boost::prior(val_it);
            a->assign(get_variable_node(var_it)->var(),
                      get_number_node(val_it)->value());
            it = stk.erase(var_it, boost::next(it));
            continue;
	  }

        if (func_def_node const * fd =
            dynamic_cast<func_def_node const *>(&*it))
        {
            fd->assign();
            it = stk.erase(it);
            continue;
        }

        ++it;
    }

    // The stack contained variable and function definitions only.
    if (stk.empty())
      return res;

    BOOST_ASSERT(stk.size() == 1);
    //return get_number_node(stk.begin())->value();
    return res;
}
  
  

  std::vector<RealType> vectorizedEval(stack stk, 
				       std::vector<RealType> &X,
				       std::vector<RealType> &Y, 
				       RealType &t ) {
    
    if (X.size() != Y.size()) {
      std::cerr << "Error: X and Y must have the same length.\n";
    }

    int n = X.size();
    std::vector<RealType> result(n, 0.0);

#pragma omp parallel for
    for ( int i = 0; i<n; ++i )
      result.at(i) = evaluate( stk, X.at(i), Y.at(i), t );
    
    return result;

  }

} // namespace yac
