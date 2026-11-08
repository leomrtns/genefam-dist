#ifndef TREESIGNAL_NEWICK_VALIDATION_H
#define TREESIGNAL_NEWICK_VALIDATION_H

#include <Python.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* The historical C reader reports malformed input through biomcmc_error(),
 * which exits the process. Validate its supported Newick subset before calling
 * it from Python, so invalid user input instead raises a Python exception. */
typedef struct {
  const char *next;
  PyObject *leaves;
  int nonroot_nodes, nonroot_lengths, has_lengths;
} NewickValidation;

static int
newick_validation_error (const char *message)
{
  PyErr_SetString (PyExc_ValueError, message);
  return 0;
}

static void
newick_validation_skip_space (NewickValidation *state)
{
  while (isspace ((unsigned char) *state->next)) state->next++;
}

static int
newick_validation_label (NewickValidation *state, int leaf)
{
  const char *start;
  PyObject *label;
  int present;

  newick_validation_skip_space (state);
  start = state->next;
  while (*state->next && !isspace ((unsigned char) *state->next) &&
         !strchr ("(),:;", *state->next)) {
    if (strchr ("'\"[]", *state->next))
      return newick_validation_error ("quoted labels and Newick comments are not supported by the native reader");
    state->next++;
  }
  if (!leaf) return 1;
  if (state->next == start) return newick_validation_error ("each leaf must have a nonempty label");
  label = PyUnicode_DecodeUTF8 (start, state->next - start, "strict");
  if (!label) return 0;
  present = PySet_Contains (state->leaves, label);
  if (present < 0) { Py_DECREF (label); return 0; }
  if (present) {
    Py_DECREF (label);
    return newick_validation_error ("leaf labels must be unique within each tree");
  }
  present = PySet_Add (state->leaves, label);
  Py_DECREF (label);
  return present == 0;
}

static int
newick_validation_subtree (NewickValidation *state, int depth)
{
  char *end;
  double length;
  int children = 0, has_length = 0;

  if (depth > 512) return newick_validation_error ("Newick nesting exceeds the native reader's supported depth (512)");
  newick_validation_skip_space (state);
  if (*state->next == '(') {
    state->next++;
    do {
      if (children) state->next++;
      if (!newick_validation_subtree (state, depth + 1)) return 0;
      children++;
      newick_validation_skip_space (state);
    } while (*state->next == ',');
    if (*state->next != ')') return newick_validation_error ("unbalanced or malformed Newick tree");
    if ((children != 2) && !((depth == 0) && (children == 3)))
      return newick_validation_error ("the native reader requires binary trees (a root trifurcation is also supported)");
    state->next++;
    if (!newick_validation_label (state, 0)) return 0;
  }
  else if (!newick_validation_label (state, 1)) return 0;

  newick_validation_skip_space (state);
  if (*state->next == ':') {
    state->next++;
    newick_validation_skip_space (state);
    length = strtod (state->next, &end);
    if (end == state->next || !isfinite (length) || length < 0.)
      return newick_validation_error ("branch lengths must be finite nonnegative numbers");
    state->next = end;
    has_length = state->has_lengths = 1;
    newick_validation_skip_space (state);
  }
  if (depth) {
    state->nonroot_nodes++;
    state->nonroot_lengths += has_length;
  }
  if (*state->next && !strchr (",);", *state->next))
    return newick_validation_error ("unexpected text in Newick tree (labels cannot contain spaces)");
  return 1;
}

static int
validate_newick_input (const char *text)
{
  NewickValidation state;
  PyObject *first_leaves = NULL;
  int equal;

  state.next = text;
  state.leaves = NULL;
  newick_validation_skip_space (&state);
  if (!*state.next) return newick_validation_error ("Newick input must contain at least one tree");
  while (*state.next) {
    state.nonroot_nodes = state.nonroot_lengths = state.has_lengths = 0;
    state.leaves = PySet_New (NULL);
    if (!state.leaves) goto failure;
    if (*state.next != '(') {
      newick_validation_error ("each Newick tree must start with '('");
      goto failure;
    }
    if (!newick_validation_subtree (&state, 0)) goto failure;
    if (PySet_Size (state.leaves) < 3) {
      newick_validation_error ("the native reader requires at least three leaves");
      goto failure;
    }
    if (state.has_lengths && state.nonroot_lengths != state.nonroot_nodes) {
      newick_validation_error ("specify a branch length on every non-root edge, or omit all branch lengths");
      goto failure;
    }
    if (first_leaves) {
      equal = PyObject_RichCompareBool (state.leaves, first_leaves, Py_EQ);
      if (equal <= 0) {
        if (!equal) newick_validation_error ("all trees in a TreeSpace must have the same leaf labels");
        goto failure;
      }
    }
    else { first_leaves = state.leaves; Py_INCREF (first_leaves); }
    Py_CLEAR (state.leaves);
    if (*state.next == ';') state.next++;
    else if (*state.next) {
      newick_validation_error ("separate Newick trees with semicolons");
      goto failure;
    }
    newick_validation_skip_space (&state);
  }
  Py_XDECREF (first_leaves);
  return 1;

failure:
  Py_XDECREF (state.leaves);
  Py_XDECREF (first_leaves);
  return 0;
}

#endif
