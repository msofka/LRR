#ifndef cdcl_macros_h
#define cdcl_macros_h

#include <vcl_typeinfo.h>
#include <vcl_iostream.h>

//:
// \file
// \brief  Macros for cdcl library.
// \author Michal Sofka
// \date   Sep 2006


//  Borrowed from rgrl
//: Macro to define type-related functions of the class
//
//  Please note, type_id() is static, and therefore
//  non-virtual. Always do class::type_id(), not
//  class_instance.type_id(). If the type_id of a class instance is
//  needed, use RTTI (run-time-type-identification), i.e.
//  typeid(class_instance).
//
//  is_type(.) checks if the instance is a type of \a classname
#define cdcl_type_macro( classname ) \
     typedef classname       self; \
     static const vcl_type_info& type_id() \
         { return typeid(self); } \
     virtual bool is_type( const vcl_type_info& type ) const\
         { return (typeid(self) == type)!=0; }
         
         
#endif // end of cdcl_macros.h
