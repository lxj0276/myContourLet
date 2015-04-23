/*
   libit - Library for basic source and channel coding functions
   Copyright (C) 2005-2005 Vivien Chappelier, Herve Jegou

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free
   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/ 

/*
  Basic types, declarations and definitions
  Copyright (C) 2005 Vivien Chappelier, Herve Jegou
*/

#ifndef __it_types_h
#define __it_types_h

#define LIBIT_VERSION(a,b,c) ((a << 16) + ((b) << 8) + (c))
#define LIBIT_VERSION_CODE   LIBIT_VERSION(0,1,0)

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef INT_MAX
#define INT_MAX 2147483647
#define INT_MIN	(-INT_MAX - 1)
#endif

#ifdef WIN32
#define inline __inline /* stoopid VC doesn't know about inline... */
#endif

/* Errors */
#define IT_EINVAL 1    /* invalid parameter */
#define IT_ENOMEM 2    /* not enough memory */
#define IT_ENOENT 3    /* no such entity */

/* Basic types */
typedef unsigned char byte;

static const size_t end = (size_t) (INT_MAX - 1);
static const size_t NULL_INDEX = (size_t) (INT_MAX);

typedef int idx_t;
/* Functions */

/* the type for function arguments */
typedef void * it_args_t;

/* a generic real valued function */
typedef double (* it_function_t)(double, it_args_t);

/* a generic integer valued function */
typedef int (* it_ifunction_t)(int, it_args_t);

/* use these macro to declare functions */
#define it_function(name)					\
it_function_args(name);						\
double __##name(double x, it_function_args(name) *it_this);	\
it_function_t const name = IT_FUNCTION(__##name);		\
double __##name(double x, it_function_args(name) *it_this)

#define it_ifunction(name)					\
it_function_args(name);						\
int __##name(int x, it_function_args(name) *it_this);		\
it_function_t const name = IT_FUNCTION(__##name);		\
int __##name(int x, it_function_args(name) *it_this)

#define it_vfunction(name)					\
it_function_args(name);						\
double __##name(vec v, it_function_args(name) *it_this);	\
it_function_t const name = IT_FUNCTION(__##name);		\
double __##name(vec v, it_function_args(name) *it_this)

/* use this macro to declare function arguments */ 
#define it_function_args(f) struct _##f##_

/* casting to an it_function */
#define IT_FUNCTION(f) ((it_function_t) f)

/* Object management */

/* These are unique random numbers identifying macros object in libit, and used for dynamic type checking. */

typedef unsigned int it_magic_t; /* magics are stored as 32-bit unique identifiers (UID) */

#define IT_MAGIC_it_object_t                  0x81e73650
#define IT_MAGIC_it_quantizer_t               0x8807fc6c 
#define IT_MAGIC_it_scalar_quantizer_t        0x48f58421
#define IT_MAGIC_it_uniform_quantizer_t       0x4a378552
#define IT_MAGIC_it_trellis_coded_quantizer_t 0x0a1d2ed3
#define IT_MAGIC_it_convolutional_code_t      0x93f58bf7
#define IT_MAGIC_it_transform2D_t             0x46e2af32
#define IT_MAGIC_it_wavelet2D_t               0xe5c8ef1a
#define IT_MAGIC_it_transform_t               0xaa709ce4
#define IT_MAGIC_it_wavelet_t                 0x8ba63118
#define IT_MAGIC_it_separable2D_t             0x299270ff
#define IT_MAGIC_it_fourier_t                 0xe0caafee

/* Objects are handled as structure, with methods defined as function pointers. Inheritance is achieved by including
   the parent structure at the beginning of the child structure. Type checking is done at runtime (in debug mode) by
   assigning unique identifiers (magics) to each object type and checking taht casts are valid (matching magics).
   This is the object structure from which all other objects inherit. It store the true type of the object.
   The virtual destructor (which can be empty) is always called when an object is destroyed.
*/

typedef struct _it_object_ {
  it_magic_t type;                                   /* true type of object */
  it_magic_t magic;    /* set to the unique ID corresponding to it_object_t */
  void (* destructor)(struct _it_object_ *it_this); /* virtual destructor */
} it_object_t;


/* Type checking */
#define it_set_magic(object, type_t) \
  do { object->magic = IT_MAGIC_##type_t; IT_OBJECT(object)->type = IT_MAGIC_##type_t; } while(0)

#define it_check_magic(object, type_t) \
  (object->magic == IT_MAGIC_##type_t)

#define it_check_type(object, type_t) \
(IT_OBJECT(object)->type == IT_MAGIC_##type_t)

#if !defined(NDEBUG) && defined(__GNUC__)
#define IT_CAST(type_t, x) 						\
({									\
  type_t *obj = (type_t *) x;						\
  									\
  if(obj->magic != IT_MAGIC_##type_t) {				\
    fprintf(stderr, "*** fatal error ***:%s:%d: object %p is not of type %s\n", __FILE__, __LINE__, (void *) obj, #type_t); \
    abort();                                                            \
  }									\
  obj;									\
})
#else /* !NDEBUG && __GNUC__ */
#define IT_CAST(type_t, x) ((type_t *) (x))
#endif /* !NDEBUG && __GNUC__ */

/* Casts x to an it_object_t, checking type in debug mode. */
#define IT_OBJECT(x) IT_CAST(it_object_t, x)

/* This is used for inheritance, to be placed at the beginning of the child class with type_t set to the parent class.
   It includes the type_t (parent) structure at the beginning of the child structure so that casting it to the
   parent is valid. An additional magic is declared to store the unique identifier of the type of the child structure
   so that dynamic type checking can be done.
*/
#define it_extends(type_t) type_t super; it_magic_t magic

/* Object creation is done by allocating a structure of the appropriate size to store the object. Then the methods
   are linked to the function pointers by the instanciate function, before the pointer is returned to the caller.
   it_newp allows for variable arguments to be passed to the constructor of the object (the instanciate function). 
*/
#define it_new(type_t) type_t##_instanciate((type_t *) malloc(sizeof(type_t)))
#define it_new_va(type_t) type_t##_instanciate((type_t *) malloc(sizeof(type_t)),
#define it_va void *)0
#define it_instanciate(type_t) type_t * type_t##_instanciate(type_t *it_this, ...)
#define it_construct(type_t) type_t##_instanciate((type_t *) it_this)
#define it_construct_va(type_t) type_t##_instanciate((type_t *) it_this,

/* Objects are destructed by calling the virtual destructor functions. The destructor of the it_object is set to the
   'free' function, which is called last and frees the memory allocated for the object. */
#define it_delete(object) IT_OBJECT(object)->destructor(IT_OBJECT(object))

/* Methods can be overloaded by declaring them as it_overload(parent_method) in the child object. */
#define it_overloaded(method) super_##method
#define it_overload(object, type, method, function)			  \
  do {									  \
    object->it_overloaded(method) = IT_CAST(type, it_this)->method; \
    IT_CAST(type, it_this)->method = function;			  \
  } while(0)

/* Instanciation of a generic object. The magic is set to the proper identifier and the destructor is initialized to
   the 'free' function. The pointer to the object in the instanciate function is always called it_this.
*/
static inline it_instanciate(it_object_t)
{

  do {
    it_this->magic = IT_MAGIC_it_object_t;
    IT_CAST(it_object_t, it_this)->type = IT_MAGIC_it_object_t;
  } while(0);


  /*  it_set_magic(it_this, it_object_t); */
  it_this->destructor = (void (*)(it_object_t *)) free;
  return(it_this);
}

/* start variable args */
#define it_new_args_start()          \
  va_list args;                      \
  va_start(args, it_this);           \
  (void) va_arg(args, void *) /* dummy */

/* stop variable args */
#define it_new_args_stop()           \
  va_end(args)

/* get next argument */
#define it_new_args_next(type_t)     \
  va_arg(args, type_t)

#ifdef __cplusplus
}
#endif /* extern "C" */
#endif
