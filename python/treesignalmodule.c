#include <Python.h>  // timestamp 2016.09.22 
// #include <genefam_dist.h>  // does not work externally (since genefam_dist.h calls other hidden headers)

static PyObject *TreesignalError;

static PyObject *
treesignal_fromtrees (PyObject *self, PyObject *args)
{
  const char *gtree_str, *splist_str;
  PyObject *arg1, *arg2, *res_tuple;

  if (!PyArg_ParseTuple(args,  "UU", &arg1, &arg2))  return NULL; 
  gtree_str  = PyBytes_AsString(PyUnicode_AsUTF8String(arg1));
  splist_str = PyBytes_AsString(PyUnicode_AsUTF8String(arg2));

  printf ("I got [%s] and [%s] \n", gtree_str, splist_str);
  // PyErr_SetString(TreesignalError, "Problem reading the string with species trees");

  res_tuple = PyTuple_New(2);
  PyTuple_SetItem(res_tuple, 0, PyFloat_FromDouble(0.));
  PyTuple_SetItem(res_tuple, 1, PyFloat_FromDouble(1.));

  return res_tuple;
}

PyMODINIT_FUNC
PyInit_treesignal(void)
{
  PyObject *m;
  static PyMethodDef TreesignalMethods[] = {
     {"fromtrees", (PyCFunction) treesignal_fromtrees, METH_VARARGS, "calculates distances given sptrees."},
     {NULL, NULL, 0, NULL}        /* Sentinel */
  };
  //static struct PyModuleDef treesignalmodule = { PyModuleDef_HEAD_INIT, "treesignal", treesignal_doc, -1, TreesignalMethods};
  static struct PyModuleDef treesignalmodule = { PyModuleDef_HEAD_INIT, "treesignal", NULL, -1, TreesignalMethods};

  m = PyModule_Create(&treesignalmodule);
  if (m == NULL) return NULL;

  TreesignalError = PyErr_NewException("treesignal.error", NULL, NULL);
  Py_INCREF(TreesignalError);
  PyModule_AddObject(m, "error", TreesignalError);
  return m;
}

int
main (int argc, char *argv[])
{
  wchar_t *program = Py_DecodeLocale(argv[0], NULL);
  if (program == NULL) { fprintf(stderr, "Fatal error: cannot decode argv[0]\n"); exit(1); }
  
  PyImport_AppendInittab("treesignal", PyInit_treesignal); /* Add a built-in module, before Py_Initialize */
  Py_SetProgramName(program); /* Pass argv[0] to the Python interpreter */
  Py_Initialize(); /* Initialize the Python interpreter.  Required. */
 
  PyMem_RawFree(program);
  return 0;
} 
