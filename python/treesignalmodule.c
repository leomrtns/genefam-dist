#include <Python.h>  // timestamp 2016.09.22 
#include <genefam/genefam_dist.h>

static PyObject *TreesignalcError;

static PyObject *
treesignalc_fromtrees (PyObject *self, PyObject *args)
{
  const char *gtree_str, *splist_str;
  PyObject *arg1, *arg2, *res_tuple;
  double *res_doublevector=NULL; /* output with distances, allocated by library function and freed here */
  int i, n_res = -1;

  if (!PyArg_ParseTuple(args,  "UU", &arg1, &arg2))  return NULL; 
  gtree_str  = PyBytes_AsString(PyUnicode_AsUTF8String(arg1));
  splist_str = PyBytes_AsString(PyUnicode_AsUTF8String(arg2));

  // printf ("I got [%s] and [%s] \n", gtree_str, splist_str); // DEBUG
  n_res = genefam_module_treesignal_fromtrees (gtree_str, splist_str, &res_doublevector);
  printf ("%d\n", n_res);
  if (n_res < 2) { PyErr_SetString(TreesignalcError, "Could not find set of species trees"); return NULL; }

  res_tuple = PyTuple_New(n_res);
  for (i = 0; i < n_res; i++) PyTuple_SetItem(res_tuple, i, PyFloat_FromDouble (res_doublevector[i]));
  if (res_doublevector) free (res_doublevector);

  return res_tuple;
}

static PyObject *
treesignalc_fromtrees_pvalue (PyObject *self, PyObject *args)
{
  const char *gtree_str, *splist_str;
  PyObject *arg1, *arg2, *res_tuple;
  double *res_doublevector=NULL; /* output with distances, allocated by library function and freed here */
  int i, n_replicates, n_res = -1;

  if (!PyArg_ParseTuple(args,  "UUi", &arg1, &arg2, &n_replicates))  return NULL; 
  gtree_str  = PyBytes_AsString(PyUnicode_AsUTF8String(arg1));
  splist_str = PyBytes_AsString(PyUnicode_AsUTF8String(arg2));
  if (n_replicates < 10) n_replicates = 10;
  // printf ("I got [%s] and [%s] \n", gtree_str, splist_str); // DEBUG
  n_res = genefam_module_treesignal_fromtrees_pvalue (gtree_str, splist_str, n_replicates, &res_doublevector);
  if (n_res < 2) { PyErr_SetString(TreesignalcError, "Could not find set of species trees"); return NULL; }

  res_tuple = PyTuple_New(n_res);
  for (i = 0; i < n_res; i++) PyTuple_SetItem(res_tuple, i, PyFloat_FromDouble (res_doublevector[i]));
  if (res_doublevector) free (res_doublevector);

  return res_tuple;
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
  //if (output_trees) free (output_trees);
  return res_string;
}

PyMODINIT_FUNC
PyInit___treesignalc(void)
{
  PyObject *m;
  static PyMethodDef TreesignalcMethods[] = {
     {"fromtrees", (PyCFunction) treesignalc_fromtrees, METH_VARARGS, "calculates distances given sptrees."},
     {"fromtrees_pvalue", (PyCFunction) treesignalc_fromtrees_pvalue, METH_VARARGS, "calculates p-values from distances given sptrees."},
     {"generate_spr_trees", (PyCFunction) treesignalc_generate_spr_trees, METH_VARARGS, "creates a chain of trees with given SPR step."},
     {NULL, NULL, 0, NULL}        /* Sentinel */
  };
  PyDoc_STRVAR(treesignalc__doc__,"lowlevel functions in C for treesignal module");

  static struct PyModuleDef treesignalcmodule = { PyModuleDef_HEAD_INIT, "treesignalc", treesignalc__doc__, -1, TreesignalcMethods};

  m = PyModule_Create(&treesignalcmodule);
  if (m == NULL) return NULL;

  TreesignalcError = PyErr_NewException("__treesignalc.error", NULL, NULL);
  Py_INCREF(TreesignalcError);
  PyModule_AddObject(m, "error", TreesignalcError);
  return m;
}

int
main (int argc, char *argv[])
{
  wchar_t *program = Py_DecodeLocale(argv[0], NULL);
  if (program == NULL) { fprintf(stderr, "Fatal error: cannot decode argv[0]\n"); exit(1); }
  
  PyImport_AppendInittab("__treesignalc", PyInit___treesignalc); /* Add a built-in module, before Py_Initialize */
  Py_SetProgramName(program); /* Pass argv[0] to the Python interpreter */
  Py_Initialize(); /* Initialize the Python interpreter.  Required. */
 
  PyMem_RawFree(program);
  return 0;
} 
