// Author(s): Wieger Wesselink. Based on the aterm library by Paul Klint and others.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/atermpp/detail/type_traits_impl.h
/// \brief Some type traits for implementations. 

#ifndef ATERMPP_DETAIL_TYPE_TRAITS_H
#define ATERMPP_DETAIL_TYPE_TRAITS_H

namespace atermpp::detail
{

// Condition for recognising types that represent containers
template < typename T >
struct is_container_impl : public std::false_type
{ };

template < typename T >
struct is_container_impl<std::list<T> > : public std::true_type
{ };

template < typename T >
struct is_container_impl< std::set< T > > : public std::true_type
{ };

template < typename T >
struct is_container_impl< std::multiset< T > > : public std::true_type
{ };

template < typename T >
struct is_container_impl< std::vector< T > > : public std::true_type
{ };

template < bool C, typename Container, typename Value >
struct lazy_check_value_type : public std::false_type
{ };

template < typename Container, typename ValueType >
struct lazy_check_value_type< true, Container, ValueType > : public
    std::is_convertible< typename Container::value_type, ValueType >
{ };

template < typename T >
struct is_set_impl : public std::false_type
  { };

template < typename T >
struct is_set_impl< std::set< T > > : public std::true_type
  { };

template < typename T >
struct is_set_impl< std::multiset< T > > :public std::true_type
  { };

} // namespace atermpp::detail



#endif // ATERMPP_DETAIL_TYPE_TRAITS_H

