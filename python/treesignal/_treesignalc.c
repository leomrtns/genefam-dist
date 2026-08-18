#include <Python.h>
#include <stdlib.h>

#include "genefam_dist.h"

static PyObject *TreesignalcError;

// TODO: replace three versions with one that accepts option about scaling
static PyObject *
treesignalc_fromtrees (PyObject *self, PyObject *args)
{
  const char *gtree_str, *splist_str;
  PyObject *res_tuple;
  double *res_doublevector=NULL; /* output with distances, allocated by library function and freed here */
  int i, n_res = -1;

  if (!PyArg_ParseTuple(args, "ss", &gtree_str, &splist_str)) return NULL;

  // printf ("I got [%s] and [%s] \n", gtree_str, splist_str); // DEBUG
  n_res = genefam_module_treesignal_fromtrees (gtree_str, splist_str, &res_doublevector);
  if (n_res < 1) { PyErr_SetString(TreesignalcError, "gene and species trees can't be compared"); return NULL; }

  res_tuple = PyTuple_New(n_res);
  for (i = 0; i < n_res; i++) PyTuple_SetItem(res_tuple, i, PyFloat_FromDouble (res_doublevector[i]));
  if (res_doublevector) free (res_doublevector);

  return res_tuple;
}

static PyObject *
treesignalc_fromtrees_rescale (PyObject *self, PyObject *args)
{
  const char *gtree_str, *splist_str;
  PyObject *res_tuple;
  double *res_doublevector=NULL; /* output with distances, allocated by library function and freed here */
  int i, n_res = -1;

  if (!PyArg_ParseTuple(args, "ss", &gtree_str, &splist_str)) return NULL;
  // printf ("I got [%s] and [%s] \n", gtree_str, splist_str); // DEBUG
  n_res = genefam_module_treesignal_fromtrees_rescale (gtree_str, splist_str, &res_doublevector);
  if (n_res < 1) { PyErr_SetString(TreesignalcError, "gene and species trees can't be compared"); return NULL; }

  res_tuple = PyTuple_New(n_res);
  for (i = 0; i < n_res; i++) PyTuple_SetItem(res_tuple, i, PyFloat_FromDouble (res_doublevector[i]));
  if (res_doublevector) free (res_doublevector);

  return res_tuple;
}

static PyObject *
treesignalc_fromtrees_pvalue (PyObject *self, PyObject *args)
{
  const char *gtree_str, *splist_str;
  PyObject *res_tuple;
  double *res_doublevector=NULL; /* output with distances, allocated by library function and freed here */
  int i, n_replicates = 1000, n_res = -1;

  if (!PyArg_ParseTuple(args, "ss|i", &gtree_str, &splist_str, &n_replicates)) return NULL;
  if (n_replicates < 10) n_replicates = 10;
  //printf ("I got [%s] and [%s] \n n = %d\n", gtree_str, splist_str, n_replicates); // DEBUG
  n_res = genefam_module_treesignal_fromtrees_pvalue (gtree_str, splist_str, n_replicates, &res_doublevector);
  if (n_res < 1) { PyErr_SetString(TreesignalcError, "gene and species trees can't be compared"); return NULL; }

  res_tuple = PyTuple_New(n_res);
  for (i = 0; i < n_res; i++) PyTuple_SetItem(res_tuple, i, PyFloat_FromDouble (res_doublevector[i]));
  if (res_doublevector) free (res_doublevector);

  return res_tuple;
}

static PyObject *
treesignalc_randomise_trees_with_spr (PyObject *self, PyObject *args)
{
  char *output_trees = NULL;
  const char *splist_str;
  PyObject *res_string;
  int n_copies = 2, n_spr = 1;  

  if (!PyArg_ParseTuple(args, "s|ii", &splist_str, &n_copies, &n_spr)) return NULL;
  output_trees = genefam_module_randomise_trees_with_spr (splist_str, n_copies, n_spr);
  if (!output_trees) { PyErr_SetString(TreesignalcError, "Could not expand trees with SPR neighbours"); return NULL; }
  
  res_string = PyUnicode_FromString ((const char *)output_trees);
  if (output_trees) free (output_trees);
  return res_string;
}

static PyObject *
treesignalc_generate_spr_trees (PyObject *self, PyObject *args)
{
  char *output_trees = NULL;
  PyObject *res_string;
  int n_leaves, n_iter, n_spr;  

  if (!PyArg_ParseTuple(args,  "iii", &n_leaves, &n_iter, &n_spr))  return NULL; 
  output_trees = genefam_module_generate_spr_trees (n_leaves, n_iter, n_spr); 
  if (!output_trees) { PyErr_SetString(TreesignalcError, "Could not create chain of SPR trees"); return NULL; }
  
  res_string = PyUnicode_FromString ((const char *)output_trees);
  if (output_trees) free (output_trees);
  return res_string;
}

PyMODINIT_FUNC
PyInit__treesignalc(void) /* it has to be named PyInit_<module name in python> */
{
  PyObject *m;
  static PyMethodDef TreesignalcMethods[] = {
     {"fromtrees", (PyCFunction) treesignalc_fromtrees, METH_VARARGS, 
      "given a set of sptrees and a gene tree, calculates a set of distances."},
     {"fromtrees_rescale", (PyCFunction) treesignalc_fromtrees_rescale, METH_VARARGS, 
      "given a set of sptrees and a gene tree, calculates a set of distances and normalise them through division by theoretical upper bounds."},
     {"fromtrees_pvalue", (PyCFunction) treesignalc_fromtrees_pvalue, METH_VARARGS, 
      "given sptrees and a genetree, returns concatenated vectors of p-values (% trees more similar) and distances rescaled to 0-1 using empirical bounds."},
     {"randomise_trees_with_spr", (PyCFunction) treesignalc_randomise_trees_with_spr, METH_VARARGS, 
      "expands set of give trees by generating SPR neighbours."},
     {"generate_spr_trees", (PyCFunction) treesignalc_generate_spr_trees, METH_VARARGS, 
      "creates a chain of trees with given SPR step."},
     {NULL, NULL, 0, NULL}        /* Sentinel */
  };
  PyDoc_STRVAR(treesignalc__doc__,"lowlevel functions in C for treesignal module");

  static struct PyModuleDef treesignalcmodule = { PyModuleDef_HEAD_INIT, "_treesignalc", treesignalc__doc__, -1, TreesignalcMethods};

  m = PyModule_Create(&treesignalcmodule);
  if (m == NULL) return NULL;

  TreesignalcError = PyErr_NewException("treesignal._treesignalc.error", NULL, NULL);
  Py_INCREF(TreesignalcError);
  PyModule_AddObject(m, "error", TreesignalcError);
  return m;
}
