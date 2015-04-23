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
  Parser for Parameter files and command line arguments
  Copyright (C) 2005 Herve Jegou
*/

#include <it/io.h>
#include <it/parser.h>
#include <stdlib.h>

/* Where the display the output of the parser */
#define parser_output stderr

/*----------------------------------------------------------------------*/
parser_t *_parser_add_string (parser_t * p, char *s)
{
  int  line_len;
  char *next_endl;
  char *line;

  if (Vec_length (p)) {
    /* continue an existing parser multiline */
    line = Vec_head (p);
  }
  else {
    /* initialize the parser by creating a new parser multiline */
    Vec_push (p, (char *) malloc (sizeof (char)));
    line = Vec_head (p);
    line[0] = 0;
  }

  while (s != NULL) {

    /* skip spaces */
    while ((*s) && (*s == ' ' || *s == '\t' || *s == '\n' || *s == ';'))
      s++;
    if (*s == '\0')
      break;

    /* Search the segment to copy and check for commentaries */
    next_endl = strpbrk (s, ";\n#");

    if (next_endl == NULL) {
      line_len = strlen (s);
      next_endl = s + strlen (s);
    }
    else {
      if (*next_endl == ';')
	next_endl++;

      line_len = (int) (next_endl - s) / sizeof (char);
      if (*next_endl == '#')
	next_endl = strchr (s, '\n');	/* Search the true end of line */
    }

    /* skip empty lines */
    if (!line_len) {
      s = next_endl;
      continue;
    }

    /* Add the text line to the parser */
    if (strchr (s, '=') && strchr (s, '=') < next_endl) {
      /* create a new multiline */
      Vec_push (p, (char *) malloc ((line_len + 2) * sizeof (char)));
      line = Vec_head (p);
      strncpy (line, s, line_len);
      line[line_len + 1] = 0;
      line[line_len] = '\n';
    }
    else {
      /* add to an existing multiline */
      line = Vec_head (p);
      line = (char *) realloc (line,
			       (strlen (line) + line_len +
				2) * sizeof (char));
      strncat (line, s, line_len);
      line[strlen (line) + 1] = 0;
      line[strlen (line)] = '\n';
    }
    Vec_head (p) = line;

    /* go to the next text line */
    s = next_endl;
  }

  return p;
}


/*----------------------------------------------------------------------*/
parser_t *_parser_add_params (parser_t * p, int argc, char **argv)
{
  int  i;
  char *line;
  for (i = 1; i < argc; i++) {
    line = (char *) malloc (strlen (argv[i]) + 1);
    strcpy (line, argv[i]);
    Vec_push (p, line);
  }
  return p;
}


/*----------------------------------------------------------------------*/
parser_t *_parser_add_file (parser_t * p, const char *filename)
{
  FILE *F;
  char *Fstring;
  int  Fsize;

  F = fopen (filename, "r+b");
  if (F == NULL) {
    it_warning ("Unable to open parameter file %s\n", filename);
    return p;
  }

  fseek (F, 0, SEEK_END);
  Fsize = ftell (F);
  fseek (F, 0, SEEK_SET);

  Fstring = (char *) malloc (Fsize + 1);
  fread (Fstring, Fsize, 1, F);

  /* replace ; with \n to handle multiple variables per line */

  p = _parser_add_string (p, Fstring);
  free (Fstring);
  return p;
}


/*----------------------------------------------------------------------*/
parser_t *parser_init (int argc, char **argv,
		       const char *filename, char *cmdline)
{
  parser_t *p = Vec_new (char *, 0);

  if (argv != NULL && argc > 0)
    p = _parser_add_params (p, argc, argv);

  if (filename != NULL)
    p = _parser_add_file (p, filename);

  if (cmdline != NULL)
    p = _parser_add_string (p, cmdline);

  return p;
}


/*----------------------------------------------------------------------*/
void parser_delete (parser_t * p)
{
  idx_t i;
  assert (p);
  for (i = 0; i < Vec_length (p); i++)
    free (p[i]);
  Vec_delete (p);
}

/*----------------------------------------------------------------------*/
void parser_print (parser_t * p)
{
  idx_t line;
  for (line = 0; line < Vec_length (p); line++)
    printf ("%s", p[line]);

}


/*----------------------------------------------------------------------*/
/* Return the line of the parser which contains the variable name       */
char *parser_var_line (const parser_t * p, const char *varname)
{
  idx_t line;
  int  found = -1;
  unsigned int n;

  for (line = 0; line < Vec_length (p); line++) {
    char const *end = strpbrk (p[line], " \t\n=");
    if (end)
      n = end - p[line];
    else
      n = strlen (p[line]);
    if (n == strlen (varname) && strncmp (p[line], varname, n) == 0) {
      found = line;
      break;
    }
  }

  if (found == -1)
    return NULL;
  else
    return p[line];
}


/*----------------------------------------------------------------------*/
int parser_exists (const parser_t * p, const char *varname)
{
  if (parser_var_line (p, varname) == NULL)
    return 0;
  return 1;
}


/*----------------------------------------------------------------------*/
int parser_get_int (const parser_t * p, const char *varname)
{
  int  return_code;
  int  r;
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;
    s++;

    return_code = sscanf (s, "%d", &r);
    if (return_code != 1)
      break;

    /* Search if the value has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      fprintf (parser_output, "%s = %d\n", varname, r);

    return r;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
int parser_get_int_verbose (const parser_t * p, const char *varname,
			    int verbose)
{
  int  return_code;
  int  r;
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;
    s++;

    return_code = sscanf (s, "%d", &r);
    if (return_code != 1)
      break;

    /* Search if the value has to be displayed on the standard output */
    if (verbose)		/* Yes */
      fprintf (parser_output, "%s = %d\n", varname, r);

    return r;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
double parser_get_double (const parser_t * p, const char *varname)
{
  double r;
  int  return_code;
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;
    s++;

    return_code = sscanf (s, "%lg", &r);
    if (return_code != 1)
      break;

    /* Search if the value has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      fprintf (parser_output, "%s = %g\n", varname, r);

    return r;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0.;
}


/*----------------------------------------------------------------------*/
double parser_get_double_verbose (const parser_t * p, const char *varname,
				  int verbose)
{
  double r;
  int  return_code;
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;
    s++;

    return_code = sscanf (s, "%lg", &r);
    if (return_code != 1)
      break;

    /* Search if the value has to be displayed on the standard output */
    if (verbose)		/* Yes */
      fprintf (parser_output, "%s = %g\n", varname, r);

    return r;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0.;
}


/*----------------------------------------------------------------------*/
byte parser_get_byte (const parser_t * p, const char *varname)
{
  int  return_code;
  int  r;
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;
    s++;

    return_code = sscanf (s, "%d", &r);
    if (return_code != 1)
      break;

    if (r > 255 || r < 0)
      it_warning ("Parameter byte " "%s" " has been rounded\n", varname);

    /* Search if the value has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      fprintf (parser_output, "%s = %d\n", varname, (int) (byte) r);

    return (byte) r;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
byte parser_get_byte_verbose (const parser_t * p, const char *varname,
			      int verbose)
{
  int  return_code;
  int  r;
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;
    s++;

    return_code = sscanf (s, "%d", &r);
    if (return_code != 1)
      break;

    if (r > 255 || r < 0)
      it_warning ("Parameter byte " "%s" " has been rounded\n", varname);

    /* Search if the value has to be displayed on the standard output */
    if (verbose)		/* Yes */
      fprintf (parser_output, "%s = %d\n", varname, (int) (byte) r);

    return (byte) r;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
char *parser_get_string (const parser_t * p, const char *varname)
{
  int  rs_len;
  int  doublequote = 1;
  char *rs;
  char *line, *s, *s_end, *s_noquote;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    /* Two cases: string is between quota or not */
    s_noquote = s;
    s = strpbrk (s, "\"");

    if (s == NULL) {		/* Assume at this point that there is no double-quote */
      doublequote = 0;
      s = s_noquote + 1;
    }
    else
      s++;

    if (doublequote) {
      s_end = strpbrk (s, "\"");
      if (s_end == NULL)
	break;
    }
    else {			/* No double quote. Stop on space, tab or \n */
      s_end = strpbrk (s, " \t\n");
      if (s_end == NULL)
	s_end = s + strlen (s);
    }

    rs_len = (int) (s_end - s);

    /* Violent malloc */
    rs = (char *) malloc (rs_len + 1);
    strncpy (rs, s, rs_len);
    rs[rs_len] = '\0';

    /* Search if the value has to be displayed on the standard output */
    if (strpbrk (s_end, ";") == NULL)	/* Yes */
      fprintf (parser_output, "%s = %s\n", varname, rs);

    return rs;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
char *parser_get_string_verbose (const parser_t * p, const char *varname,
				 int verbose)
{
  int  rs_len;
  char *rs;
  char *line, *s, *s_end, *s_noquote;
  int  doublequote = 1;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    /* Two cases: string is between quota or not */
    s_noquote = s;
    s = strpbrk (s, "\"");

    if (s == NULL) {		/* Assume at this point that there is no double-quote */
      doublequote = 0;
      s = s_noquote + 1;
    }
    else
      s++;

    if (doublequote) {
      s_end = strpbrk (s, "\"");
      if (s_end == NULL)
	break;
    }
    else {			/* No double quote. Stop on space, tab or \n */
      s_end = strpbrk (s, " \t\n");
      if (s_end == NULL)
	s_end = s + strlen (s);
    }

    rs_len = (int) (s_end - s);

    /* Violent malloc */
    rs = (char *) malloc (rs_len + 1);
    strncpy (rs, s, rs_len);
    rs[rs_len] = '\0';

    /* Search if the value has to be displayed on the standard output */
    if (verbose)		/* Yes */
      fprintf (parser_output, "%s = %s\n", varname, rs);

    return rs;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
cplx parser_get_cplx (const parser_t * p, const char *varname)
{
  cplx c;
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_cplx (s, &c);

    /* Search if the vector has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = %z\n", varname, c);

    return c;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return cplx_0;
}


/*----------------------------------------------------------------------*/
cplx parser_get_cplx_verbose (const parser_t * p, const char *varname,
			      int verbose)
{
  cplx c;
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_cplx (s, &c);

    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = %z\n", varname, c);

    return c;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return cplx_0;
}


/*----------------------------------------------------------------------*/
vec parser_get_vec (const parser_t * p, const char *varname)
{
  vec  v = vec_new_alloc (0, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_vec (s, &v);
    if (!s)
      it_error ("Unable to read a vector of double in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = $v\n", varname, v);

    return v;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
vec parser_get_vec_verbose (const parser_t * p, const char *varname,
			    int verbose)
{
  vec  v = vec_new_alloc (0, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_vec (s, &v);
    if (!s)
      it_error ("Unable to read a vector of double in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = $v\n", varname, v);

    return v;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
cvec parser_get_cvec (const parser_t * p, const char *varname)
{
  cvec v = cvec_new_alloc (0, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_cvec (s, &v);
    if (!s)
      it_error ("Unable to read a vector of complexes in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = $z\n", varname, v);

    return v;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return v;
}


/*----------------------------------------------------------------------*/
cvec parser_get_cvec_verbose (const parser_t * p, const char *varname,
			      int verbose)
{
  cvec v = cvec_new_alloc (0, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_cvec (s, &v);
    if (!s)
      it_error ("Unable to read a vector of complexes in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = $z\n", varname, v);

    return v;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return v;
}


/*----------------------------------------------------------------------*/
ivec parser_get_ivec (const parser_t * p, const char *varname)
{
  ivec v = ivec_new_alloc (0, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_ivec (s, &v);
    if (!s)
      it_error ("Unable to read a vector of integers in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = $d\n", varname, v);

    return v;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
ivec parser_get_ivec_verbose (const parser_t * p, const char *varname,
			      int verbose)
{
  ivec v = ivec_new_alloc (0, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_ivec (s, &v);
    if (!s)
      it_error ("Unable to read a vector of integers in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = $d\n", varname, v);

    return v;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
bvec parser_get_bvec (const parser_t * p, const char *varname)
{
  bvec v = bvec_new_alloc (0, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_bvec (s, &v);
    if (!s)
      it_error ("Unable to read a vector of bytes in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = $b\n", varname, v);

    return v;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
bvec parser_get_bvec_verbose (const parser_t * p, const char *varname,
			      int verbose)
{
  bvec v = bvec_new_alloc (0, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;

    s = it_read_bvec (s, &v);
    if (!s)
      it_error ("Unable to read a vector of bytes in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = $b\n", varname, v);

    return v;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
mat parser_get_mat (const parser_t * p, const char *varname)
{
  mat  m = mat_new_alloc (0, 0, 10, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;


    s = it_read_mat (s, &m);
    if (!s)
      it_error ("Unable to read a matrix of double in string %s\n", s);

    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = #m\n", varname, m);

    return m;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
mat parser_get_mat_verbose (const parser_t * p, const char *varname,
			    int verbose)
{
  mat  m = mat_new_alloc (0, 0, 10, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;


    s = it_read_mat (s, &m);
    if (!s)
      it_error ("Unable to read a matrix of double in string %s\n", s);

    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = #m\n", varname, m);

    return m;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
cmat parser_get_cmat (const parser_t * p, const char *varname)
{
  cmat m = cmat_new_alloc (0, 0, 10, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;


    s = it_read_cmat (s, &m);
    if (!s)
      it_error ("Unable to read a matrix of complexes in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = #.3z\n", varname, m);

    return m;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
cmat parser_get_cmat_verbose (const parser_t * p, const char *varname,
			      int verbose)
{
  cmat m = cmat_new_alloc (0, 0, 10, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;


    s = it_read_cmat (s, &m);
    if (!s)
      it_error ("Unable to read a matrix of complexes in string %s\n", s);

    /* Search if the matrix has to be displayed on the standard output */
    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = #.3z\n", varname, m);

    return m;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
imat parser_get_imat (const parser_t * p, const char *varname)
{
  imat m = imat_new_alloc (0, 0, 10, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;


    s = it_read_imat (s, &m);
    if (!s)
      it_error ("Unable to read a matrix of double in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = #d\n", varname, m);

    return m;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
imat parser_get_imat_verbose (const parser_t * p, const char *varname,
			      int verbose)
{
  imat m = imat_new_alloc (0, 0, 10, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;


    s = it_read_imat (s, &m);
    if (!s)
      it_error ("Unable to read a matrix of double in string %s\n", s);

    /* Search if the matrix has to be displayed on the standard output */
    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = #d\n", varname, m);

    return m;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
bmat parser_get_bmat (const parser_t * p, const char *varname)
{
  bmat m = bmat_new_alloc (0, 0, 10, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;


    s = it_read_bmat (s, &m);
    if (!s)
      it_error ("Unable to read a matrix of double in string %s\n", s);

    /* Search if the vector has to be displayed on the standard output */
    if (strpbrk (s, ";") == NULL)	/* Yes */
      it_fprintf (parser_output, "%s = #d\n", varname, m);

    return m;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}


/*----------------------------------------------------------------------*/
bmat parser_get_bmat_verbose (const parser_t * p, const char *varname,
			      int verbose)
{
  bmat m = bmat_new_alloc (0, 0, 10, 10);
  char *line, *s;
  line = parser_var_line (p, varname);

  while (1) {
    if (line == NULL)
      break;

    s = strpbrk (line, "=");
    if (s == NULL)
      break;


    s = it_read_bmat (s, &m);
    if (!s)
      it_error ("Unable to read a matrix of double in string %s\n", s);

    /* Search if the matrix has to be displayed on the standard output */
    if (verbose)		/* Yes */
      it_fprintf (parser_output, "%s = #d\n", varname, m);

    return m;
  }

  it_warning ("Unable to retrieve variable " "%s" "\n", varname);
  return 0;
}
