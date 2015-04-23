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
  I/O functions for the libit library
  Copyright (C) 2005 Herve Jegou
*/

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <it/vec.h>
#include <it/mat.h>
#include <it/io.h>
#include <it/types.h>
#include <it/cplx.h>
#include <it/poly.h>


/* Default format for vector (called by format $v)  */
#define it_fmt_max_length (20)
static char it_printf_vec_default_fmt[it_fmt_max_length + 1] = "%.3lf";
static char it_printf_mat_default_fmt[it_fmt_max_length + 1] = "$9.3lf";


/*------------------------------------------------------------------------------*/
void it_set_vec_default_fmt (char *fmt)
{
  assert (strlen (fmt) < it_fmt_max_length);
  strcpy (it_printf_vec_default_fmt, fmt);
}


/*------------------------------------------------------------------------------*/
void it_set_mat_default_fmt (char *fmt)
{
  assert (strlen (fmt) < it_fmt_max_length);
  strcpy (it_printf_vec_default_fmt, fmt);
}


/*------------------------------------------------------------------------------*/
void it_printf (char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  it_vfprintf (stdout, fmt, ap);
  va_end (ap);
}

/*------------------------------------------------------------------------------*/
void it_vprintf (char *fmt, va_list args)
{
  it_vfprintf (stdout, fmt, args);
}

/*------------------------------------------------------------------------------*/
void it_warning (char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
#if defined(__GNUC__)
  fprintf (stderr, "## IT Warning : %s:%d ##  ", __FILE__, __LINE__);
#else
  fprintf (stderr, "## IT Warning ##  ");
#endif
  it_vfprintf (stderr, fmt, ap);
  fprintf (stderr, "\n");
  va_end (ap);
}


/*------------------------------------------------------------------------------*/
void it_error (char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
#if defined(__GNUC__)
  fprintf (stderr, "## IT Error : %s:%d ##  ", __FILE__, __LINE__);
#else
  fprintf (stderr, "## IT Error ##  ");
#endif
  it_vfprintf (stderr, fmt, ap);
  fprintf (stderr, "\n");
  va_end (ap);
  exit (1);
}


/*---------------------------------------------------------------------------*/
void __it_assert (int a, const char *filename, int line, char *message)
{
  if (!a) {

#if defined(__GNUC__)
    fprintf (stderr, "## IT Assert : %s:%d ##  %s\n", filename, line,
	     message);
#else
    fprintf (stderr, "## IT Assert  ##  %s\n", message);
#endif
    abort ();
  }
}


/*------------------------------------------------------------------------------*/
void it_fprintf (FILE * output, char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  it_vfprintf (output, fmt, ap);
  va_end (ap);
}


/*------------------------------------------------------------------------------*/
void it_vfprintf (FILE * output, char *fmt, va_list ap)
{
  char fmt_opt[100];		/* A string for the format and one for the options */
  int  fmt_opt_len;

  /* Declare a variable for each kind of type that may be displayed */
  char c, *s;
  int  d;
  idx_t i;
  double f;
  cplx z;
  void *p;
  vec  P;

  vec  v;
  ivec iv;
  bvec bv;
  cvec cv;
  pvec Pv;

  mat  m;
  imat im;
  bmat bm;
  cmat cm;
  pmat Pm;


  for (; *fmt; fmt++) {

    /*--------------------         common display       ----------------------*/
    switch (*fmt) {
    case '%':
      assert (*(fmt + 1));	/* Unfinished sequence. Format string must continue. */

      if (*(fmt + 1) == '%') {
	fprintf (output, "%%");
	fmt++;
	continue;
      }

      /* Find the length of the segment consisting to the format option */
      fmt_opt_len = strcspn (fmt, "sdxcfgpzP");


      memcpy (fmt_opt, fmt, fmt_opt_len + 1);
      fmt_opt[fmt_opt_len + 1] = '\0';
      fmt = fmt + fmt_opt_len;

      switch (*fmt) {
      case 's':		/* string */
	s = va_arg (ap, char *);
	fprintf (output, fmt_opt, s);
	break;

      case 'd':		/* int */
	d = va_arg (ap, int);
	fprintf (output, fmt_opt, d);
	break;

      case 'x':		/* hexadecimal int */
	d = va_arg (ap, int);
	fprintf (output, "0x");
	fprintf (output, fmt_opt, d);
	break;

      case 'c':		/* char */
	/* need a cast here since va_arg only
	   takes fully promoted types */
	c = (char) va_arg (ap, int);
	fprintf (output, fmt_opt, c);
	break;

      case 'f':
	f = va_arg (ap, double);
	fprintf (output, fmt_opt, f);
	break;

      case 'g':
	f = va_arg (ap, double);
	fprintf (output, fmt_opt, f);
	break;

      case 'z':		/* complex */
	/* This one is a particular case: there are two values to display */
	z = va_arg (ap, cplx);
	fmt_opt[strlen (fmt_opt) - 1] = 'f';

	if (creal (z) == 0 && cimag (z) == 0) {
	  fprintf (output, "0");
	  break;
	}

	if (creal (z) != 0)
	  fprintf (output, fmt_opt, creal (z));

	if (cimag (z) > 0) {
	  fprintf (output, "+");
	  fprintf (output, fmt_opt, cimag (z));
	  fprintf (output, "i");
	}
	else if (cimag (z) < 0) {
	  fprintf (output, fmt_opt, cimag (z));
	  fprintf (output, "i");
	}
	break;

      case 'p':
	p = va_arg (ap, void *);
	fprintf (output, fmt_opt, p);
	break;

      case 'P':		/* polynomial */
	P = va_arg (ap, vec);

	if (vec_length (P) == 0) {
	  fprintf (output, "(0)");
	  break;
	}

	fmt_opt[fmt_opt_len] = 'f';

	fprintf (output, "(");
	for (i = 0; i < vec_length (P) - 1; i++) {
	  fprintf (output, fmt_opt, P[i]);
	  if (i)
	    fprintf (output, "*X^%d", i);
	  fprintf (output, " + ");
	}
	fprintf (output, fmt_opt, P[i]);
	if (i)
	  fprintf (output, "*X^%d", i);
	fprintf (output, ")");
	break;

      default:
	fprintf (stderr, "## it_vfprintf : unrecognized type ## ");
	exit (1);
      }
      break;

      /*--------------------         Vector display       ----------------------*/
    case '$':
      if (*(fmt + 1) == '$') {
	fprintf (output, "$");
	fmt++;
	continue;
      }

      /* Find the length of the segment consisting to the format option */
      fmt_opt_len = strcspn (fmt, "idxbfgvzP");

      memcpy (fmt_opt, fmt, fmt_opt_len + 1);
      fmt_opt[0] = '%';
      fmt_opt[fmt_opt_len + 1] = '\0';
      fmt = fmt + fmt_opt_len;

      switch (*fmt) {
      case 'v':		/* Default representation of the vector. Other format are ignored */
	strcpy (fmt_opt, it_printf_vec_default_fmt);
      case 'f':
      case 'g':
	v = va_arg (ap, vec);
	if (Vec_header (v).element_size != sizeof (double)) {
	  it_warning
	    ("Incompatible formatting type in it_fprintf. Double base type was expected\n");
	  return;
	}
	fprintf (output, "[");
	for (i = 0; i < vec_length (v); i++) {
	  it_fprintf (output, fmt_opt, v[i]);
	  if (i < vec_length (v) - 1)
	    fprintf (output, " ");
	}
	fprintf (output, "]");
	break;

      case 'i':
	fmt_opt[strlen (fmt_opt) - 1] = 'd';
      case 'x':
      case 'd':
	iv = va_arg (ap, ivec);
	if (Vec_header (iv).element_size != sizeof (int)) {
	  it_warning
	    ("Incompatible formatting type in it_fprintf. Integer base type was expected\n");
	  return;
	}
	fprintf (output, "[");
	for (i = 0; i < ivec_length (iv); i++) {
	  it_fprintf (output, fmt_opt, iv[i]);
	  if (i < ivec_length (iv) - 1)
	    fprintf (output, " ");
	}
	fprintf (output, "]");
	break;

      case 'b':
	fmt_opt[strlen (fmt_opt) - 1] = 'd';	/* Displayed as an integer */
      case 'c':
	bv = va_arg (ap, bvec);
	if (Vec_header (bv).element_size != sizeof (byte)) {
	  it_warning
	    ("Incompatible formatting type in it_fprintf. Byte base type was expected\n");
	  return;
	}
	fprintf (output, "[");
	for (i = 0; i < bvec_length (bv); i++) {
	  it_fprintf (output, fmt_opt, bv[i]);
	  if (i < bvec_length (bv) - 1)
	    fprintf (output, " ");
	}
	fprintf (output, "]");
	break;

      case 'z':
	cv = va_arg (ap, cvec);
	if (Vec_header (cv).element_size != sizeof (cplx)) {
	  it_warning
	    ("Incompatible formatting type in it_fprintf. Complex base type was expected\n");
	  return;
	}
	fmt_opt[strlen (fmt_opt) - 1] = 'f';
	fprintf (output, "[");
	for (i = 0; i < cvec_length (cv); i++) {
	  z = cv[i];

	  if (creal (z) == 0 && cimag (z) == 0) {
	    fprintf (output, "0");
	  }
	  else {

	    if (creal (z) != 0)
	      fprintf (output, fmt_opt, creal (z));

	    if (cimag (z) > 0) {
	      fprintf (output, "+");
	      fprintf (output, fmt_opt, cimag (z));
	      fprintf (output, "i");
	    }
	    else if (cimag (z) < 0) {
	      fprintf (output, fmt_opt, cimag (z));
	      fprintf (output, "i");
	    }
	  }
	  if (i < cvec_length (cv) - 1)
	    fprintf (output, " ");
	}
	fprintf (output, "]");
	break;

      case 'P':
	Pv = va_arg (ap, pvec);

	if (Vec_header (Pv).element_size != sizeof (vec)) {
	  it_warning
	    ("Incompatible formatting type in it_fprintf. Polynomial base type was expected\n");
	  return;
	}

	fprintf (output, "[");
	for (i = 0; i < Vec_length (Pv); i++) {
	  P = Pv[i];

	  it_fprintf (output, fmt_opt, Pv[i]);
	  if (i < Vec_length (Pv) - 1)
	    fprintf (output, " ");
	}
	fprintf (output, "]");
	break;

      default:
	fprintf (stderr, "## it_vfprintf : unrecognized type ## ");
	exit (1);
      }
      break;

      /*--------------------         Matrix display       ----------------------*/
    case '#':
      if (*(fmt + 1) == '#') {
	fprintf (output, "#");
	fmt++;
	continue;
      }

      /* Find the length of the segment consisting to the format option */
      fmt_opt_len = strcspn (fmt, "idxbfgmzP");

      memcpy (fmt_opt, fmt, fmt_opt_len + 1);
      fmt_opt[0] = '$';
      fmt_opt[fmt_opt_len + 1] = '\0';
      fmt = fmt + fmt_opt_len;

      switch (*fmt) {
      case 'm':		/* Default representation of the matrix. Other format are ignored */
	strcpy (fmt_opt, it_printf_mat_default_fmt);
      case 'f':
      case 'g':
	m = va_arg (ap, mat);
	fprintf (output, "[");
	for (i = 0; i < mat_height (m); i++) {
	  it_fprintf (output, fmt_opt, m[i]);
	  if (i < mat_height (m) - 1)
	    fprintf (output, "\n ");
	}
	fprintf (output, "]");
	break;

      case 'i':
	fmt_opt[strlen (fmt_opt) - 1] = 'd';
      case 'x':
      case 'd':
	im = va_arg (ap, imat);
	fprintf (output, "[");
	for (i = 0; i < imat_height (im); i++) {
	  it_fprintf (output, fmt_opt, im[i]);
	  if (i < imat_height (im) - 1)
	    fprintf (output, "\n ");
	}
	fprintf (output, "]");
	break;

      case 'b':
	fmt_opt[strlen (fmt_opt) - 1] = 'b';	/* Displayed as integer */
      case 'c':
	bm = va_arg (ap, bmat);
	fprintf (output, "[");
	for (i = 0; i < bmat_height (bm); i++) {
	  it_fprintf (output, fmt_opt, bm[i]);
	  if (i < bmat_height (bm) - 1)
	    fprintf (output, "\n ");
	}
	fprintf (output, "]");
	break;


      case 'z':
	cm = va_arg (ap, cmat);
	fprintf (output, "[");
	for (i = 0; i < cmat_height (cm); i++) {
	  it_fprintf (output, fmt_opt, cm[i]);
	  if (i < cmat_height (cm) - 1)
	    fprintf (output, "\n ");
	}
	fprintf (output, "]");
	break;

      case 'P':
	Pm = va_arg (ap, pmat);
	fprintf (output, "[");
	for (i = 0; i < Mat_height (Pm); i++) {
	  it_fprintf (output, fmt_opt, Pm[i]);
	  if (i < Mat_height (Pm) - 1)
	    fprintf (output, "\n ");
	}
	fprintf (output, "]");
	break;

      default:
	fprintf (stderr, "## it_vfprintf : unrecognized type ## ");
	exit (1);
      }
      break;

      /*---------------         Direct display       ------------------*/
    default:
      fprintf (output, "%c", *fmt);
    }
  }
}


/*---------------------------------------------------------------------*/
char *it_read_double (char *s, double *p_val)
{
  int  nb_char_number = strspn (s, "0123456789.e+-");
  int  nb_char_not_delim = strcspn (s, " ,\t\n];");
  char *buf;

  if (nb_char_number != nb_char_not_delim)
    it_error ("Unable to read a double in string %s\n", s);

  buf = (char *) malloc (nb_char_number + 1);
  buf[nb_char_number] = '\0';
  strncpy (buf, s, nb_char_number);

  sscanf (buf, "%lf", p_val);
  free (buf);
  return s + nb_char_not_delim;
}


/*---------------------------------------------------------------------*/
char *it_read_float (char *s, float *p_val)
{
  int  nb_char_number = strspn (s, "0123456789.e+-");
  int  nb_char_not_delim = strcspn (s, " ,\t\n];");
  char *buf;

  if (nb_char_number != nb_char_not_delim)
    it_error ("Unable to read a float in string %s\n", s);

  buf = (char *) malloc (nb_char_number + 1);
  buf[nb_char_number] = '\0';
  strncpy (buf, s, nb_char_number);

  sscanf (buf, "%f", p_val);
  free (buf);
  return s + nb_char_not_delim;
}


/*---------------------------------------------------------------------*/
char *it_read_cplx (char *s, cplx * p_val)
{
  int  nb_char_number;
  int  nb_char_not_delim;
  int  read_something = 0, r;
  char *buf;
  creal (*p_val) = 0;
  cimag (*p_val) = 0;

  while (*s == ' ' || *s == '\t')
    s++;

  if (*s == '+' || *s == '-') {
    nb_char_number = strspn (s + 1, "0123456789.e") + 1;
    nb_char_not_delim = strcspn (s + 1, " ,\t\n];i+-") + 1;
  }
  else {
    nb_char_number = strspn (s, "0123456789.e");
    nb_char_not_delim = strcspn (s, " ,\t\n];i+-");
  }

  if (nb_char_number != nb_char_not_delim)
    it_error ("Unable to read a complex number in string %s\n", s);

  if (*s == 'i') {		/* Special case */
    nb_char_not_delim = strcspn (s + 1, " ,\t\n];");
    if (nb_char_not_delim)
      it_warning ("Something nasty occured while reading the complex %s\n",
		  s);

    cimag (*p_val) = 1;
    return s + 1;
  }

  if (s[nb_char_not_delim] != 'i') {	/* It is the real part */
    buf = (char *) malloc (nb_char_number + 1);
    buf[nb_char_number] = '\0';
    strncpy (buf, s, nb_char_number);
    sscanf (buf, "%lf", &(creal (*p_val)));
    free (buf);
    s += nb_char_not_delim;
    read_something = 1;

    /* Ready to read the (optional) imaginary part */
    nb_char_number = strspn (s, "0123456789.e+-");
    nb_char_not_delim = strcspn (s, " ,\t\n];i");
  }

  if (s[nb_char_not_delim] == 'i') {	/* It is the complex part */
    buf = (char *) malloc (nb_char_number + 1);
    buf[nb_char_number] = '\0';
    strncpy (buf, s, nb_char_number);
    r = sscanf (buf, "%lf", &(cimag (*p_val)));
    if (r == 0) {		/* Bad return value for sscanf */
      if (*s == 'i' || (*s == '+' && s[1] == 'i'))
	cimag (*p_val) = 1;
      else if (*s == '-' && s[1] == 'i')
	cimag (*p_val) = -1;
      else
	it_warning
	  ("Something wrong happened while reading the complex. Unsafe result.\n");
    }
    free (buf);
    s += nb_char_not_delim + 1;
    read_something = 1;
  }

  it_assert (read_something, "Invalid complex");
  return s;
}


/*---------------------------------------------------------------------*/
char *it_read_int (char *s, int *p_val)
{
  int  nb_char_number = strspn (s, "0123456789+-");
  int  nb_char_not_delim = strcspn (s, " ,\t\n];");
  char *buf;

  if (nb_char_number != nb_char_not_delim)
    it_error ("Bad integer numbers in string %s\n", s);

  buf = (char *) malloc (nb_char_number + 1);
  buf[nb_char_number] = '\0';
  strncpy (buf, s, nb_char_number);

  sscanf (buf, "%d", p_val);
  free (buf);
  return s + nb_char_not_delim;
}


/*---------------------------------------------------------------------*/
char *it_read_vec (char *s, vec * v)
{
  char *s_start, *s_end;
  char *s_first_b = strpbrk (s, "[");
  char *s_last_b = strpbrk (s, "]");
  char *s_first_number = strpbrk (s, "0123456789.e+-");
  char *s_last_semicol = strpbrk (s, "\n];");
  double d;

  if (s_first_b < s_last_semicol && s_last_b > s_last_semicol && s_first_b)
    return (NULL);

  if (s_first_b == NULL) {
    if (s_last_semicol == NULL)
      s_last_semicol = s + strlen (s);

    s_start = s_first_number;
    s_end = s_last_semicol;
  }
  else {
    s_start = strpbrk (s_first_b, "0123456789.e+-");
    s_end = s_last_b;
  }

  if (s_first_number == NULL)
    s_start = s_end;

  vec_set_length (*v, 0);

  for (s = s_start; s && s < s_end;) {
    s = it_read_double (s, &d);
    vec_push (*v, d);
    s = strpbrk (s, "0123456789.e+-");
  }
  return s_end + 1;
}


/*---------------------------------------------------------------------*/
char *it_read_fvec (char *s, vec * v)
{
  char *s_start, *s_end;
  char *s_first_b = strpbrk (s, "[");
  char *s_last_b = strpbrk (s, "]");
  char *s_first_number = strpbrk (s, "0123456789.e+-");
  char *s_last_semicol = strpbrk (s, "\n];");
  float d;

  if (s_first_b < s_last_semicol && s_last_b > s_last_semicol && s_first_b)
    return (NULL);

  if (s_first_b == NULL) {
    if (s_last_semicol == NULL)
      s_last_semicol = s + strlen (s);

    s_start = s_first_number;
    s_end = s_last_semicol;
  }
  else {
    s_start = strpbrk (s_first_b, "0123456789.e+-");
    s_end = s_last_b;
  }

  if (s_first_number == NULL)
    s_start = s_end;

  vec_set_length (*v, 0);

  for (s = s_start; s && s < s_end;) {
    s = it_read_float (s, &d);
    vec_push (*v, d);
    s = strpbrk (s, "0123456789.e+-");
  }
  return s_end + 1;
}


/*---------------------------------------------------------------------*/
char *it_read_cvec (char *s, cvec * v)
{
  char *s_start, *s_end;
  char *s_first_b = strpbrk (s, "[");
  char *s_last_b = strpbrk (s, "]");
  char *s_first_number = strpbrk (s, "0123456789.e+-i");
  char *s_last_semicol = strpbrk (s, "\n];");
  cplx c;

  if (s_first_b < s_last_semicol && s_last_b > s_last_semicol && s_first_b)
    return (NULL);

  if (s_first_b == NULL) {
    if (s_last_semicol == NULL)
      s_last_semicol = s + strlen (s);

    s_start = s_first_number;
    s_end = s_last_semicol;
  }
  else {
    s_start = strpbrk (s_first_b, "0123456789.e+-i");
    s_end = s_last_b;
  }

  if (s_first_number == NULL)
    s_start = s_end;

  cvec_set_length (*v, 0);

  for (s = s_start; s && s < s_end;) {
    s = it_read_cplx (s, &c);
    cvec_push (*v, c);
    s = strpbrk (s, "0123456789.e+-i");
  }
  return s_end + 1;
}


/*---------------------------------------------------------------------*/
char *it_read_ivec (char *s, ivec * v)
{
  char *s_start, *s_end;
  char *s_first_b = strpbrk (s, "[");
  char *s_last_b = strpbrk (s, "]");
  char *s_first_number = strpbrk (s, "0123456789+-");
  char *s_last_semicol = strpbrk (s, "\n];");
  int  d;

  if (s_first_b < s_last_semicol && s_last_b > s_last_semicol && s_first_b)
    return (NULL);

  if (s_first_b == NULL) {
    if (s_last_semicol == NULL)
      s_last_semicol = s + strlen (s);

    s_start = s_first_number;
    s_end = s_last_semicol;
  }
  else {
    s_start = strpbrk (s_first_b, "0123456789+-");
    s_end = s_last_b;
  }

  if (s_first_number == NULL)
    s_start = s_end;

  ivec_set_length (*v, 0);

  for (s = s_start; s && s < s_end;) {
    s = it_read_int (s, &d);
    ivec_push (*v, d);
    s = strpbrk (s, "0123456789+-");
  }
  return s_end + 1;
}


/*---------------------------------------------------------------------*/
char *it_read_bvec (char *s, bvec * v)
{
  char *s_start, *s_end;
  char *s_first_b = strpbrk (s, "[");
  char *s_last_b = strpbrk (s, "]");
  char *s_first_number = strpbrk (s, "+-0123456789");
  char *s_last_semicol = strpbrk (s, "\n];");
  int  d;

  if (s_first_b < s_last_semicol && s_last_b > s_last_semicol && s_first_b)
    return (NULL);

  if (s_first_b == NULL) {
    if (s_last_semicol == NULL)
      s_last_semicol = s + strlen (s);

    s_start = s_first_number;
    s_end = s_last_semicol;
  }
  else {
    s_start = strpbrk (s_first_b, "+-0123456789");
    s_end = s_last_b;
  }

  if (s_first_number == NULL)
    s_start = s_end;

  bvec_set_length (*v, 0);

  for (s = s_start; s && s < s_end;) {

    s = it_read_int (s, &d);
    if (d < 0 || d >= (1 << (8 * sizeof (byte))))
      it_warning ("Invalid byte value : %d. Casting to %d...\n", d, (byte) d);

    bvec_push (*v, (byte) d);
    s = strpbrk (s, "+-0123456789");
  }
  return s_end + 1;
}


/*---------------------------------------------------------------------*/
char *it_read_mat (char *s, mat * m)
{
  char *s_start, *s_end;
  vec  v;
  int  c, l;

  /* find the opening bracket */
  s_start = strchr (s, '[');
  /* find the closing bracket */
  s_end = NULL;
  if (s_start) {
    l = strlen (s);
    c = 0;
    for (s_end = s_start; s_end < s + l; s_end++) {
      if (*s_end == '[')
	c++;
      if (*s_end == ']')
	c--;
      if (c == 0)
	break;
    }
    if (s_end >= s + l || *s_end != ']')
      s_end = NULL;
  }

  if (!s_start || !s_end)
    it_error ("Unable to read a matrix of double in string %s\n", s);

  mat_set_height (*m, 0);

  /* read the vectors */
  s = s_start + 1;
  while (s && s < s_end) {
    v = vec_new (10);
    s = it_read_vec (s, &v);
    if (s)
      mat_push_row (*m, v);
    else
      vec_delete (v);
  }

  return s_end + 1;
}


/*---------------------------------------------------------------------*/
char *it_read_fmat (char *s, mat * m)
{
  char *s_start, *s_end;
  vec  v;
  int  c, l;

  /* find the opening bracket */
  s_start = strchr (s, '[');
  /* find the closing bracket */
  s_end = NULL;
  if (s_start) {
    l = strlen (s);
    c = 0;
    for (s_end = s_start; s_end < s + l; s_end++) {
      if (*s_end == '[')
	c++;
      if (*s_end == ']')
	c--;
      if (c == 0)
	break;
    }
    if (s_end >= s + l || *s_end != ']')
      s_end = NULL;
  }

  if (!s_start || !s_end)
    it_error ("Unable to read a matrix of float in string %s\n", s);

  mat_set_height (*m, 0);

  /* read the vectors */
  s = s_start + 1;
  while (s && s < s_end) {
    v = vec_new (10);
    s = it_read_vec (s, &v);
    if (s)
      mat_push_row (*m, v);
    else
      vec_delete (v);
  }

  return s_end + 1;
}


/*---------------------------------------------------------------------*/
char *it_read_cmat (char *s, cmat * m)
{
  char *s_start, *s_end;
  cvec v;
  int  c, l;

  /* find the opening bracket */
  s_start = strchr (s, '[');
  /* find the closing bracket */
  s_end = NULL;
  if (s_start) {
    l = strlen (s);
    c = 0;
    for (s_end = s_start; s_end < s + l; s_end++) {
      if (*s_end == '[')
	c++;
      if (*s_end == ']')
	c--;
      if (c == 0)
	break;
    }
    if (s_end >= s + l || *s_end != ']')
      s_end = NULL;
  }

  if (!s_start || !s_end)
    it_error ("Unable to read a matrix of double in string %s\n", s);

  cmat_set_height (*m, 0);

  /* read the vectors */
  s = s_start + 1;
  while (s && s < s_end) {
    v = cvec_new (10);
    s = it_read_cvec (s, &v);
    if (s)
      cmat_push_row (*m, v);
    else
      cvec_delete (v);
  }

  return s_end + 1;
}

/*---------------------------------------------------------------------*/
char *it_read_imat (char *s, imat * m)
{
  char *s_start, *s_end;
  ivec v;
  int  c, l;

  /* find the opening bracket */
  s_start = strchr (s, '[');
  /* find the closing bracket */
  s_end = NULL;
  if (s_start) {
    l = strlen (s);
    c = 0;
    for (s_end = s_start; s_end < s + l; s_end++) {
      if (*s_end == '[')
	c++;
      if (*s_end == ']')
	c--;
      if (c == 0)
	break;
    }
    if (s_end >= s + l || *s_end != ']')
      s_end = NULL;
  }

  if (!s_start || !s_end)
    it_error ("Unable to read a matrix of double in string %s\n", s);

  imat_set_height (*m, 0);

  /* read the vectors */
  s = s_start + 1;
  while (s && s < s_end) {
    v = ivec_new (10);
    s = it_read_ivec (s, &v);
    if (s)
      imat_push_row (*m, v);
    else
      ivec_delete (v);
  }

  return s_end + 1;
}

/*---------------------------------------------------------------------*/
char *it_read_bmat (char *s, bmat * m)
{
  char *s_start, *s_end;
  bvec v;
  int  c, l;

  /* find the opening bracket */
  s_start = strchr (s, '[');
  /* find the closing bracket */
  s_end = NULL;
  if (s_start) {
    l = strlen (s);
    c = 0;
    for (s_end = s_start; s_end < s + l; s_end++) {
      if (*s_end == '[')
	c++;
      if (*s_end == ']')
	c--;
      if (c == 0)
	break;
    }
    if (s_end >= s + l || *s_end != ']')
      s_end = NULL;
  }

  if (!s_start || !s_end)
    it_error ("Unable to read a matrix of double in string %s\n", s);

  bmat_set_height (*m, 0);

  /* read the vectors */
  s = s_start + 1;
  while (s && s < s_end) {
    v = bvec_new (10);
    s = it_read_bvec (s, &v);
    if (s)
      bmat_push_row (*m, v);
    else
      bvec_delete (v);
  }

  return s_end + 1;
}

/*---------------------------------------------------------------------*/
vec vec_new_string (char *s)
{
  vec  v = vec_new_alloc (0, 10);
  if (!it_read_vec (s, &v))
    it_error ("Unable to read a vector of double in string %s\n", s);

  return v;
}


/*---------------------------------------------------------------------*/
vec fvec_new_string (char *s)
{
  vec  v = vec_new_alloc (0, 10);
  if (!it_read_fvec (s, &v))
    it_error ("Unable to read a vector of float in string %s\n", s);

  return v;
}


/*---------------------------------------------------------------------*/
ivec ivec_new_string (char *s)
{
  ivec v = ivec_new_alloc (0, 10);
  if (!it_read_ivec (s, &v))
    it_error ("Unable to read a vector of integers in string %s\n", s);

  return v;
}


/*---------------------------------------------------------------------*/
bvec bvec_new_string (char *s)
{
  bvec v = bvec_new_alloc (0, 10);
  if (!it_read_bvec (s, &v))
    it_error ("Unable to read a vector of bytes in string %s\n", s);

  return v;
}


/*---------------------------------------------------------------------*/
cvec cvec_new_string (char *s)
{
  cvec v = cvec_new_alloc (0, 10);
  if (!it_read_cvec (s, &v))
    it_error ("Unable to read a vector of complexes in string %s\n", s);

  return v;
}


/*---------------------------------------------------------------------*/
mat mat_new_string (char *s)
{
  mat  m = mat_new_alloc (0, 0, 10, 10);
  if (!it_read_mat (s, &m))
    it_error ("Unable to read a matrix of double in string %s\n", s);

  return m;
}


/*---------------------------------------------------------------------*/
mat fmat_new_string (char *s)
{
  mat  m = mat_new_alloc (0, 0, 10, 10);
  if (!it_read_fmat (s, &m))
    it_error ("Unable to read a matrix of double in string %s\n", s);

  return m;
}


/*---------------------------------------------------------------------*/
cmat cmat_new_string (char *s)
{
  cmat m = cmat_new_alloc (0, 0, 10, 10);
  if (!it_read_cmat (s, &m))
    it_error ("Unable to read a matrix of double in string %s\n", s);

  return m;
}

/*---------------------------------------------------------------------*/
imat imat_new_string (char *s)
{
  imat m = imat_new_alloc (0, 0, 10, 10);
  if (!it_read_imat (s, &m))
    it_error ("Unable to read a matrix of double in string %s\n", s);

  return m;
}

/*---------------------------------------------------------------------*/
bmat bmat_new_string (char *s)
{
  bmat m = bmat_new_alloc (0, 0, 10, 10);
  if (!it_read_bmat (s, &m))
    it_error ("Unable to read a matrix of double in string %s\n", s);

  return m;
}


/*---------------------------------------------------------------------*/
/*   PNM related functions                                             */
/*---------------------------------------------------------------------*/

/* Suppress the additional white characters and concatenate the comments 
 to the string comments                                                */
static void pnm_read_comments (FILE * F, char *comments, int length)
{
  int  r;
  int  comment_pos = 0;

  do
    r = fgetc (F);
  while (r == (int) '\n' || r == (int) ' ' || r == (int) '\t');

  while (r == (int) '#')
    /* Comments, read characters until the end of line                 */
    do {
      r = fgetc (F);
      if (comment_pos < length - 1)
	comments[comment_pos++] = (char) r;
    }
    while (r != (int) '\n' && !feof (F));

  if (r >= '0' && r <= '9')
    ungetc (r, F);

  if (length)
    comments[comment_pos] = '\0';
}


/* Read/Write the header for the pnm file format */
static int pnm_read_header (FILE * file, char *p_pnm_type, int *p_width,
			    int *p_height, int *p_max_val, char *comments,
			    int length)
{

  int  r;

  /* Set default values if the parsing would fails before the end */
  *p_pnm_type = 0;
  *p_width = 0;
  *p_height = 0;
  *p_max_val = 0;

  /* Read the 'P' identifier required for a pnm file */
  r = fgetc (file);
  if (r != (int) 'P') {
    it_warning ("Invalid format file: not a pnm file\n");
    return 0;
  }

  *p_pnm_type = fgetc (file);
  if (*p_pnm_type < '1' || *p_pnm_type > '6') {
    it_warning ("Unknown pnm type");
    return 0;
  }

  /* Read the comments */
  pnm_read_comments (file, comments, length);

  fscanf (file, "%d", p_width);
  fscanf (file, "%d", p_height);


  /* Maximal value does not exist in PBM files */
  if (*p_pnm_type == '2' || *p_pnm_type == '3' ||
      *p_pnm_type == '5' || *p_pnm_type == '6')
    fscanf (file, "%d", p_max_val);
  else
    *p_max_val = 1;

  /* According to the pnm specification, the maximal value should not
     be greater than 65536 and lower than 0                             */
  if (*p_max_val >= 65536 || *p_max_val <= 0) {
    it_warning ("Invalid maximum number in pnm header\n");
    return 0;
  }

  /* Eat the last whitespace */
  r = fgetc (file);

  /* For type P5 and P6, the value have to be lower than 255            */
  if ((*p_pnm_type == '5' || *p_pnm_type == '6') && *p_max_val > 255) {
    it_warning ("Invalid maximum number in pnm header\n");
    return 0;
  }
  return 1;
}


static int pnm_write_header (FILE * file, char type, int width, int height,
			     int max_val, const char *comments)
{
  fprintf (file, "P%c\n#%s\n%d %d\n", type, comments, width, height);
  if (type == '2' || type == '3' || type == '5' || type == '6')
    fprintf (file, "%d\n", max_val);
  return 1;
}


/*----------------------------------------------------------------------*/
char pnm_type (const char *filename)
{
  FILE *F = fopen (filename, "r+b");
  char pnm_type;
  int  width, height, max_val;

  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return 0;
  }

  pnm_read_header (F, &pnm_type, &width, &height, &max_val, NULL, 0);
  fclose (F);
  return pnm_type;
}


/*----------------------------------------------------------------------*/
int pnm_info (const char *filename, char *p_pnm_type, int *p_width,
	      int *p_height, int *p_max_val, char *comments, int length)
{
  FILE *F = fopen (filename, "r+b");
  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return 0;
  }

  pnm_read_header (F, p_pnm_type, p_width, p_height, p_max_val, comments,
		   length);
  fclose (F);
  return 1;
}


/*----------------------------------------------------------------------*/
mat mat_pgm_read (const char *filename)
{
  FILE *F = fopen (filename, "r+b");
  char type;
  int  width, height, max_val, i, j;
  mat  m;

  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return 0;
  }

  /* Return a void matrix */
  if (!pnm_read_header (F, &type, &width, &height, &max_val, NULL, 0))
    return mat_new (0, 0);

  m = mat_new (height, width);
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++)
      m[i][j] = (double) fgetc (F);

  fclose (F);
  return m;
}


/*----------------------------------------------------------------------*/
imat imat_pgm_read (const char *filename)
{
  FILE *F = fopen (filename, "r+b");
  char type;
  int  width, height, max_val, i, j;
  imat m;

  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return 0;
  }

  /* Return a void matrix */
  if (!pnm_read_header (F, &type, &width, &height, &max_val, NULL, 0))
    return imat_new (0, 0);

  m = imat_new (height, width);
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++)
      m[i][j] = (int) fgetc (F);

  fclose (F);
  return m;
}


/*----------------------------------------------------------------------*/
/* Write a matrix of double as a pgm image file                         */
int mat_pgm_write (const char *filename, mat m)
{
  idx_t i, j;
  double v;
  FILE *F = fopen (filename, "w+b");
  pnm_write_header (F, '5', mat_width (m), mat_height (m), 255,
		    "Generated by libit");

  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++) {
      v = m[i][j] + 0.5;
      if (v < 0)
	v = 0;
      if (v > 255)
	v = 255;
      fputc ((int) v, F);
    }

  fclose (F);
  return 1;
}


/*----------------------------------------------------------------------*/
/* Write a matrix of integers as a pgm file                             */
int imat_pgm_write (const char *filename, imat m)
{
  idx_t i, j;
  FILE *F = fopen (filename, "w+b");
  pnm_write_header (F, '5', imat_width (m), imat_height (m), 255,
		    "Generated by libit");

  for (i = 0; i < imat_height (m); i++)
    for (j = 0; j < imat_width (m); j++)
      fputc (m[i][j], F);

  fclose (F);
  return 1;
}


/*---------------------------------------------------------------------*/
/*   WAV related functions                                             */
/*---------------------------------------------------------------------*/

static byte wav_read8 (FILE * file)
{
  int  c;
  if ((c = fgetc (file)) == EOF) {
    it_warning ("unexpected end of file\n");
    return (0);
  }
  return ((byte) c);
}


static unsigned short wav_read16 (FILE * file)
{
  unsigned short v = 0;

  v |= wav_read8 (file);
  v |= wav_read8 (file) << 8;

  return (v);
}


static unsigned long wav_read32 (FILE * file)
{
  unsigned long v = 0;

  v |= wav_read8 (file);
  v |= wav_read8 (file) << 8;
  v |= wav_read8 (file) << 16;
  v |= wav_read8 (file) << 24;

  return (v);
}


static int wav_read_sample (FILE * file, int depth)
{
  int  v = 0;

  v |= wav_read8 (file);
  if (depth == 8)
    return (v - 128);		/* what a stupid file format */

  if (depth > 8)
    v |= wav_read8 (file) << 8;
  if (depth > 16)
    v |= wav_read8 (file) << 8;
  if (depth > 24)
    v |= wav_read8 (file) << 8;

  /* sign extend */
  v <<= 32 - (((depth + 7) >> 3) << 3);
  v >>= 32 - depth;

  return (v);
}


static void wav_write8 (FILE * file, byte c)
{
  fputc (c, file);
}

static void wav_write16 (FILE * file, unsigned short c)
{
  wav_write8 (file, (byte) c);
  wav_write8 (file, (byte) (c >> 8));
}


static void wav_write32 (FILE * file, unsigned long c)
{
  wav_write8 (file, (byte) c);
  wav_write8 (file, (byte) (c >> 8));
  wav_write8 (file, (byte) (c >> 16));
  wav_write8 (file, (byte) (c >> 24));
}


static void wav_write_sample (FILE * file, int depth, int v)
{
  if (depth == 8)
    v += 128;			/* what a stupid file format */

  /* align to the next byte msb */
  v <<= (((depth + 7) >> 3) << 3) - depth;

  wav_write8 (file, (byte) v);
  if (depth > 8)
    wav_write8 (file, (byte) (v >> 8));
  if (depth > 16)
    wav_write8 (file, (byte) (v >> 16));
  if (depth > 24)
    wav_write8 (file, (byte) (v >> 24));
}


#define FOURCC(a,b,c,d) ((int)(a) | ((int)(b) << 8) | ((int)(c) << 16) | ((int)(d) << 24) )

/* Read/Write the header for the wav file format */
static int wav_read_header (FILE * file, int *p_channels, int *p_srate,
			    int *p_depth, int *p_length)
{
  int  r;

  /* Set default values if the parsing would fails before the end */
  *p_channels = 0;
  *p_srate = 0;
  *p_depth = 0;
  *p_length = 0;

  /* Read the 'RIFF' identifier required for a wave file */
  r = wav_read32 (file);
  if (r != FOURCC ('R', 'I', 'F', 'F')) {
    it_warning ("Invalid format file: not a wav file\n");
    return 0;
  }

  /* Read the size of the following data */
  r = wav_read32 (file);

  /* Read the 'WAVE' identifier required for a wave file */
  r = wav_read32 (file);
  if (r != FOURCC ('W', 'A', 'V', 'E')) {
    it_warning ("Invalid format file: not a wav file\n");
    return 0;
  }

  /* Read the 'fmt ' identifier */
  r = wav_read32 (file);
  if (r != FOURCC ('f', 'm', 't', ' ')) {
    it_warning ("Invalid format file: not a wav file\n");
    return 0;
  }

  /* Read the size of the following data */
  r = wav_read32 (file);

  /* Read the format */
  r = wav_read16 (file);
  if (r != 1) {
    it_warning ("WAV: unsupported compressed file\n");
    return 0;
  }

  /* Read the number of channels */
  r = wav_read16 (file);
  if (r <= 0) {
    it_warning ("Invalid format file: invalid number of channels\n");
    return 0;
  }
  *p_channels = r;

  /* Read the sampling rate */
  r = wav_read32 (file);
  if (r <= 0) {
    it_warning ("Invalid format file: invalid sampling rate\n");
    return 0;
  }
  *p_srate = r;

  /* Read the average number of bytes per second (ignore) */
  r = wav_read32 (file);

  /* Read the block alignment (ignore) */
  r = wav_read16 (file);

  /* Read the number of bits per sample */
  r = wav_read16 (file);
  if (r <= 0) {
    it_warning ("Invalid format file: invalid number of bits per sample\n");
    return 0;
  }
  *p_depth = r;

  /* Read the 'data' identifier */
  r = wav_read32 (file);
  if (r != FOURCC ('d', 'a', 't', 'a')) {
    it_warning ("Invalid format file: not a wav file\n");
    return 0;
  }

  /* Read the size of the data chunk */
  r = wav_read32 (file);
  *p_length = r / ((*p_depth / 8) * *p_channels);

  return 1;
}


static void wav_write_header (FILE * file, int channels, int srate,
			      int depth, int length)
{
  int  data_size = length * ((depth + 7) / 8) * channels;

  /* Write the 'RIFF' identifier required for a wave file */
  wav_write32 (file, FOURCC ('R', 'I', 'F', 'F'));

  /* Read the size of the following data */
  wav_write32 (file, data_size + 36);

  /* Write the 'WAVE' identifier */
  wav_write32 (file, FOURCC ('W', 'A', 'V', 'E'));

  /* Write the 'fmt ' identifier */
  wav_write32 (file, FOURCC ('f', 'm', 't', ' '));

  /* Write the size of the following data */
  wav_write32 (file, 16);

  /* Write the format (uncompressed = 1) */
  wav_write16 (file, 1);

  /* Write the number of channels */
  wav_write16 (file, (unsigned short) channels);

  /* Write the sampling rate */
  wav_write32 (file, srate);

  /* Write the average number of bytes per second */
  wav_write32 (file, srate * channels * ((depth + 7) / 8));

  /* Write the block alignment */
  wav_write16 (file, (unsigned short) (channels * ((depth + 7) / 8)));

  /* Write the number of bits per sample */
  wav_write16 (file, (unsigned short) depth);

  /* Read the 'data' identifier */
  wav_write32 (file, FOURCC ('d', 'a', 't', 'a'));

  /* Write the size of the data chunk */
  wav_write32 (file, data_size);
}


/*----------------------------------------------------------------------*/
int wav_info (const char *filename, int *p_channels, int *p_srate,
	      int *p_depth, int *p_length)
{
  FILE *F = fopen (filename, "r+b");
  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return 0;
  }

  wav_read_header (F, p_channels, p_srate, p_depth, p_length);
  fclose (F);
  return 1;
}


/*----------------------------------------------------------------------*/
mat mat_wav_read (const char *filename)
{
  FILE *F = fopen (filename, "r+b");
  int  channels, srate, depth, length;
  int  i, j;
  mat  m;

  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return 0;
  }

  /* Return a void matrix */
  if (!wav_read_header (F, &channels, &srate, &depth, &length))
    return mat_new (0, 0);

  m = mat_new (channels, length);
  for (j = 0; j < length; j++)
    for (i = 0; i < channels; i++)
      m[i][j] = (double) wav_read_sample (F, depth);

  fclose (F);
  return m;
}


/*----------------------------------------------------------------------*/
imat imat_wav_read (const char *filename)
{
  FILE *F = fopen (filename, "r+b");
  int  channels, srate, depth, length;
  int  i, j;
  imat m;

  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return 0;
  }

  /* Return a void matrix */
  if (!wav_read_header (F, &channels, &srate, &depth, &length))
    return imat_new (0, 0);

  m = imat_new (channels, length);
  for (j = 0; j < length; j++)
    for (i = 0; i < channels; i++)
      m[i][j] = wav_read_sample (F, depth);

  fclose (F);
  return m;
}


/*----------------------------------------------------------------------*/
void mat_wav_write (const char *filename, mat m, int srate, int depth)
{
  FILE *F = fopen (filename, "w+b");
  int  channels, length;
  int  i, j;
  double s;
  double min, max;

  channels = mat_height (m);
  length = mat_width (m);

  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return;
  }

  wav_write_header (F, channels, srate, depth, length);

  max = (1 << (depth - 1)) - 1;
  min = -max - 1;

  for (j = 0; j < length; j++)
    for (i = 0; i < channels; i++) {
      s = m[i][j];
      if (s > 0)
	s += 0.5;
      if (s < 0)
	s -= 0.5;
      if (s < min)
	s = min;
      if (s > max)
	s = max;
      wav_write_sample (F, depth, (int) s);
    }

  fclose (F);
}

/*----------------------------------------------------------------------*/
void imat_wav_write (const char *filename, imat m, int srate, int depth)
{
  FILE *F = fopen (filename, "w+b");
  int  channels, length;
  int  i, j;

  channels = imat_height (m);
  length = imat_width (m);

  if (!F) {
    it_printf ("Unable to open file %s\n", filename);
    return;
  }

  wav_write_header (F, channels, srate, depth, length);

  for (j = 0; j < length; j++)
    for (i = 0; i < channels; i++)
      wav_write_sample (F, depth, m[i][j]);

  fclose (F);
}


/*----------------------------------------------------------------------*/
/* Write and read matrix in a packed format                             */

void vec_fwrite (FILE * stream, vec v)
{
  int  n = vec_length (v), ret;
  ret = fwrite (&n, sizeof (n), 1, stream);
  assert (ret == 1);

  if (n) {
    ret = fwrite (v, sizeof (*v), n, stream);
    assert (ret == n);
  }
}


void fvec_fwrite (FILE * stream, vec v)
{
  int  n = vec_length (v), ret;
  float ftmp;
  ret = fwrite (&n, sizeof (n), 1, stream);
  assert (ret == 1);

  if (n) {
    int i;
    for (i = 0 ; i < n ; i++) {
      ftmp = (float) v[i];
      ret = fwrite (&ftmp, sizeof (ftmp), 1, stream);
      assert (ret == 1);
    }
  }
}


/*----------------------------------------------------------------------*/
void bvec_fwrite (FILE * stream, bvec v)
{
  int  n = bvec_length (v), ret;
  ret = fwrite (&n, sizeof (n), 1, stream);
  assert (ret == 1);

  if (n) {
    ret = fwrite (v, sizeof (*v), n, stream);
    assert (ret == n);
  }
}


/*----------------------------------------------------------------------*/
void ivec_fwrite (FILE * stream, ivec v)
{
  int  n = ivec_length (v), ret;
  ret = fwrite (&n, sizeof (n), 1, stream);
  assert (ret == 1);

  if (n) {
    ret = fwrite (v, sizeof (*v), n, stream);
    assert (ret == n);
  }
}


/*----------------------------------------------------------------------*/
void cvec_fwrite (FILE * stream, cvec v)
{
  int  n = cvec_length (v), ret;
  ret = fwrite (&n, sizeof (n), 1, stream);
  assert (ret == 1);

  if (n) {
    ret = fwrite (v, sizeof (*v), n, stream);
    assert (ret == n);
  }
}



/*----------------------------------------------------------------------*/
int vec_fread (FILE * stream, vec v)
{
  int  n, ret;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return 0;

  if (n != vec_length (v))
    return 0;

  if (n) {
    ret = fread (v, sizeof (*v), n, stream);
    if (ret != n)
      return 0;
  }

  return 1;
}


/*----------------------------------------------------------------------*/
int fvec_fread (FILE * stream, vec v)
{
  int  n, ret;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return 0;

  if (n != vec_length (v))
    return 0;

  if (n) {
    int i;
    float ftmp;
    for (i = 0 ; i < n ; i++) {
      ret = fread (&ftmp, sizeof (ftmp), 1, stream);
      if (ret != 1)
	return 0;
      v[i] = ftmp;
    }
  }

  return 1;
}



/*----------------------------------------------------------------------*/
int ivec_fread (FILE * stream, ivec v)
{
  int  n, ret;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return 0;

  if (n != ivec_length (v))
    return 0;

  if (n) {
    ret = fread (v, sizeof (*v), n, stream);
  }

  return 1;
}


/*----------------------------------------------------------------------*/
int bvec_fread (FILE * stream, bvec v)
{
  int  n, ret;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return 0;

  if (n != bvec_length (v))
    return 0;

  if (n) {
    ret = fread (v, sizeof (*v), n, stream);
    if (ret != n)
      return 0;
  }

  return 1;
}


/*----------------------------------------------------------------------*/
int cvec_fread (FILE * stream, cvec v)
{
  int  n, ret;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return 0;

  if (n != cvec_length (v))
    return 0;

  if (n) {
    ret = fread (v, sizeof (*v), n, stream);
    if (ret != n)
      return 0;
  }

  return 1;
}



/*----------------------------------------------------------------------*/
vec vec_new_fread (FILE * stream)
{
  int  n, ret;
  vec  v;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return NULL;

  v = vec_new (n);
  if (n) {
    ret = fread (v, sizeof (*v), n, stream);
    if (ret != n)
      return NULL;
  }
  return v;
}


/*----------------------------------------------------------------------*/
vec fvec_new_fread (FILE * stream)
{
  int  n, ret;
  vec  v;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return NULL;

  v = vec_new (n);
  if (n) {
    int i;
    float ftmp;
    for (i = 0 ; i < n ; i++) {
      ret = fread (&ftmp, sizeof (ftmp), 1, stream);
      if (ret != 1)
	return NULL;
      v[i] = ftmp;
    }
  }
  return v;
}


/*----------------------------------------------------------------------*/
bvec bvec_new_fread (FILE * stream)
{
  int  n, ret;
  bvec v;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return NULL;

  v = bvec_new (n);
  if (n) {
    ret = fread (v, sizeof (*v), n, stream);
    if (ret != n)
      return NULL;
  }
  return v;

}


/*----------------------------------------------------------------------*/
ivec ivec_new_fread (FILE * stream)
{
  int  n, ret;
  ivec v;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return NULL;

  v = ivec_new (n);

  if (n) {
    ret = fread (v, sizeof (*v), n, stream);

    if (ret != n)
      return NULL;
  }
  return v;
}


/*----------------------------------------------------------------------*/
cvec cvec_new_fread (FILE * stream)
{
  int  n, ret;
  cvec v;
  ret = fread (&n, sizeof (n), 1, stream);
  if (ret != 1)
    return NULL;

  v = cvec_new (n);
  ret = fread (v, sizeof (*v), n, stream);
  if (ret != n)
    return NULL;
  return v;
}


/*----------------------------------------------------------------------*/
void mat_fwrite (FILE * stream, mat m)
{
  int  w = mat_width (m);
  int  h = mat_height (m), i, ret;

  /* Write the dimension of the matrix (first the height) */
  ret = fwrite (&h, sizeof (h), 1, stream);
  assert (ret == 1);

  fwrite (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  if (w && h)
    for (i = 0; i < h;  i++) {
      ret = fwrite (m[i], sizeof (**m), w, stream);
      assert (ret == w);
    }
}


/*----------------------------------------------------------------------*/
void fmat_fwrite (FILE * stream, mat m)
{
  int  w = mat_width (m);
  int  h = mat_height (m), j, i, ret;

  /* Write the dimension of the matrix (first the height) */
  ret = fwrite (&h, sizeof (h), 1, stream);
  assert (ret == 1);

  fwrite (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  if (w && h)
    for (i = 0; i < h; i++) {
      ret = 0;
      for (j = 0 ; j < w ; j++) {
	float tmp = (float) m[i][j];
	ret += fwrite (&tmp, sizeof (tmp), 1, stream);
      }
      assert (ret == w);
    }
}


/*----------------------------------------------------------------------*/
void bmat_fwrite (FILE * stream, bmat m)
{
  int  w = bmat_width (m);
  int  h = bmat_height (m), i, ret;

  /* Write the dimension of the matrix (first the height) */
  ret = fwrite (&h, sizeof (h), 1, stream);
  assert (ret == 1);
  ret = fwrite (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  if (w && h)
    for (i = 0; i < h; i++) {
      ret = fwrite (m[i], sizeof (**m), w, stream);
      assert (ret == w);
    }
}


/*----------------------------------------------------------------------*/
void imat_fwrite (FILE * stream, imat m)
{
  int  w = imat_width (m);
  int  h = imat_height (m), i, ret;

  /* Write the dimension of the matrix (first the height) */
  ret = fwrite (&h, sizeof (h), 1, stream);
  assert (ret == 1);
  ret = fwrite (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  if (w && h)
    for (i = 0; i < h; i++) {
      ret = fwrite (m[i], sizeof (**m), w, stream);
      assert (ret == w);
    }
}


/*----------------------------------------------------------------------*/
void cmat_fwrite (FILE * stream, cmat m)
{
  int  w = cmat_width (m);
  int  h = cmat_height (m), i, ret;

  /* Write the dimension of the matrix (first the height) */
  ret = fwrite (&h, sizeof (h), 1, stream);
  assert (ret == 1);

  ret = fwrite (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  if (w && h)
    for (i = 0; i < h; i++) {
      ret = fwrite (m[i], sizeof (**m), w, stream);
      assert (ret == w);
    }
}



/*----------------------------------------------------------------------*/
int vec_sread (void * buffer, vec v)
{
  int n;

  n = * ((int *) buffer);
  buffer += sizeof (int);  /* Position the buffer at the beginning of the date */
  
  vec_init (v, (double *) buffer, n);

  return sizeof (int) + n * sizeof (double);
}


/*----------------------------------------------------------------------*/
int ivec_sread (void * buffer, ivec v)
{
  int n;

  n = * ((int *) buffer);
  buffer += sizeof (int);  /* Position the buffer at the beginning of the date */
  
  ivec_init (v, (int *) buffer, n);

  return sizeof (int) + n * sizeof (int);
}


/*----------------------------------------------------------------------*/
mat mat_new_fread (FILE * stream)
{
  int  w, h, i, ret;
  mat  m;
  ret = fread (&h, sizeof (h), 1, stream);
  assert (ret == 1);
  ret = fread (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  m = mat_new (h, w);

  if (w && h)
    for (i = 0; i < h; i++) {
      ret = fread (m[i], sizeof (**m), w, stream);
      assert (ret == w);
    }
  return m;
}


/*----------------------------------------------------------------------*/
bmat bmat_new_fread (FILE * stream)
{
  int  w, h, i, ret;
  bmat m;
  ret = fread (&h, sizeof (h), 1, stream);
  assert (ret == 1);
  ret = fread (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  m = bmat_new (h, w);

  if (w && h)
    for (i = 0; i < h; i++) {
      ret = fread (m[i], sizeof (**m), w, stream);
      assert (ret == w);
    }
  return m;
}


/*----------------------------------------------------------------------*/
imat imat_new_fread (FILE * stream)
{
  int  w, h, i, ret;
  imat m;
  ret = fread (&h, sizeof (h), 1, stream);
  assert (ret == 1);
  ret = fread (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  m = imat_new (h, w);

  if (w && h)
    for (i = 0; i < h; i++) {
      ret = fread (m[i], sizeof (**m), w, stream);
      assert (ret == w);
    }
  return m;
}


/*----------------------------------------------------------------------*/
cmat cmat_new_fread (FILE * stream)
{
  int  w, h, i, ret;
  cmat m;
  ret = fread (&h, sizeof (h), 1, stream);
  assert (ret == 1);
  ret = fread (&w, sizeof (w), 1, stream);
  assert (ret == 1);

  m = cmat_new (h, w);

  if (w && h)
    for (i = 0; i < h; i++) {
      ret = fread (m[i], sizeof (**m), w, stream);
      assert (ret == w);
    }
  return m;
}


/*----------------------------------------------------------------------*/
bvec bvec_file_read_bits (const char *filename, int nb_max)
{
  FILE *F = fopen (filename, "r+b");
  byte b = 0;
  bvec v, buf;
  int  i, return_code;

  if (F == NULL)
    it_error ("Unable to open file %s\n", filename);

  /* The whole file has to be read. We have to find the length of this file */
  if (nb_max <= 0) {
    fseek (F, 0, SEEK_END);	/* Reach the end of the file */
    nb_max = ftell (F) * 8;	/* vector length = Length of the file * 8 */
    rewind (F);			/* Come back to the start of the file */
  }

  /* Read a vector of characters from the file */
  buf = bvec_new ((nb_max + 7) / 8);
  return_code = fread (buf, 1, bvec_length (buf), F);
  fclose (F);

  /* Verify that we have read the correct number of bits, otherwise quit */
  if (return_code != (nb_max + 7) / 8)
    it_error ("Incorrect number of bits read from file %s\n", filename);

  v = bvec_new (nb_max);

  for (i = 0; i < nb_max; i++) {
    if (i % 8 == 0)
      b = buf[i / 8];

    v[i] = (b >> (7 - (i % 8))) & 1;
  }

  bvec_delete (buf);
  return v;
}


/*----------------------------------------------------------------------*/
void bvec_file_write_bits (const char *filename, bvec v)
{
  FILE *F = fopen (filename, "w+b");
  bvec buf = bvec_new_zeros ((bvec_length (v) + 7) / 8);
  int  i, return_code;

  if (F == NULL)
    it_error ("Unable to create the file %s\n", filename);

  for (i = 0; i < bvec_length (v); i++) {

    /* Note that v must be a vector of bits */
    buf[i / 8] ^= (v[i] & 1) << (7 - i % 8);
  }

  return_code = fwrite (buf, bvec_length (buf), 1, F);

  printf ("bvec_length(v)   = %d\n", bvec_length (v));
  printf ("bvec_length(buf) = %d\n", bvec_length (buf));
  printf ("return_code      = %d\n", return_code);

  if (!return_code)
    it_error ("Could not write a bit vector in file %d\n", filename);
  fclose (F);
  bvec_delete (buf);
}


