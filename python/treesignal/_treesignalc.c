#include <Python.h>
#include <stdlib.h>

#include "genefam_dist.h"

typedef struct {
  PyObject_HEAD
  topology_space space;
  int rooted;
} TreeSpaceObject;

static PyObject *TreesignalcError;
static PyTypeObject TreeSpaceType;

static PyObject *
distance_vector_to_tuple (double *values, int length)
{
  PyObject *item, *result;
  int i;

  if (length < 1) {
    if (values) free (values);
    PyErr_SetString (TreesignalcError, "gene and species trees cannot be compared");
    return NULL;
  }
  result = PyTuple_New (length);
  if (!result) {
    free (values);
    return NULL;
  }
  for (i = 0; i < length; i++) {
    item = PyFloat_FromDouble (values[i]);
    if (!item) {
      free (values);
      Py_DECREF (result);
      return NULL;
    }
    PyTuple_SET_ITEM (result, i, item);
  }
  free (values);
  return result;
}

static PyObject *
TreeSpace_new (PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
  TreeSpaceObject *self = (TreeSpaceObject*) type->tp_alloc (type, 0);
  if (self) {
    self->space = NULL;
    self->rooted = 0;
  }
  return (PyObject*) self;
}

static int
TreeSpace_init (TreeSpaceObject *self, PyObject *args, PyObject *kwargs)
{
  static char *keywords[] = {"newick", "rooted", NULL};
  const char *newick;
  int rooted = 0;
  topology_space parsed;

  if (!PyArg_ParseTupleAndKeywords (args, kwargs, "s|p:TreeSpace", keywords, &newick, &rooted)) return -1;
  parsed = genefam_module_parse_newick_trees (newick, rooted != 0);
  if (!parsed) {
    PyErr_SetString (TreesignalcError, "could not parse the Newick tree data");
    return -1;
  }
  if (self->space) del_topology_space (self->space);
  self->space = parsed;
  self->rooted = rooted;
  return 0;
}

static void
TreeSpace_dealloc (TreeSpaceObject *self)
{
  if (self->space) del_topology_space (self->space);
  Py_TYPE (self)->tp_free ((PyObject*) self);
}

static PyObject *
TreeSpace_repr (TreeSpaceObject *self)
{
  if (!self->space) return PyUnicode_FromString ("<treesignal._treesignalc.TreeSpace (uninitialized)>");
  return PyUnicode_FromFormat (
      "<treesignal._treesignalc.TreeSpace trees=%d distinct=%d rooted=%s>",
      self->space->ntrees,
      self->space->ndistinct,
      self->rooted ? "True" : "False");
}

static PyObject *
TreeSpace_get_ntrees (TreeSpaceObject *self, void *closure)
{
  if (!self->space) return PyLong_FromLong (0);
  return PyLong_FromLong (self->space->ntrees);
}

static PyObject *
TreeSpace_get_ndistinct (TreeSpaceObject *self, void *closure)
{
  if (!self->space) return PyLong_FromLong (0);
  return PyLong_FromLong (self->space->ndistinct);
}

static PyObject *
TreeSpace_get_rooted (TreeSpaceObject *self, void *closure)
{
  return PyBool_FromLong (self->rooted);
}

static PyGetSetDef TreeSpace_getset[] = {
  {"ntrees", (getter) TreeSpace_get_ntrees, NULL, "Number of input trees.", NULL},
  {"ndistinct", (getter) TreeSpace_get_ndistinct, NULL, "Number of distinct native topologies.", NULL},
  {"rooted", (getter) TreeSpace_get_rooted, NULL, "Whether root location distinguishes topologies.", NULL},
  {NULL, NULL, NULL, NULL, NULL}
};

PyDoc_STRVAR (
    TreeSpace_doc,
    "TreeSpace(newick, rooted=False)\n--\n\n"
    "Owning native representation of one or more trees. Distance functions use "
    "the first gene tree and all distinct species trees.\n\n"
    "Newick data is parsed once during construction. The underlying C topology_space "
    "is retained until this Python object is destroyed, so object-based distance "
    "functions avoid repeated string serialization and parsing.");

static PyTypeObject TreeSpaceType = {
  PyVarObject_HEAD_INIT (NULL, 0)
  .tp_name = "treesignal._treesignalc.TreeSpace",
  .tp_basicsize = sizeof (TreeSpaceObject),
  .tp_dealloc = (destructor) TreeSpace_dealloc,
  .tp_repr = (reprfunc) TreeSpace_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = TreeSpace_doc,
  .tp_getset = TreeSpace_getset,
  .tp_init = (initproc) TreeSpace_init,
  .tp_new = TreeSpace_new,
};

static int
parse_tree_spaces (PyObject *args, TreeSpaceObject **gene, TreeSpaceObject **species, const char *function_name)
{
  if (!PyArg_ParseTuple (args, "O!O!", &TreeSpaceType, gene, &TreeSpaceType, species)) return 0;
  if ((!(*gene)->space) || (!(*species)->space)) {
    PyErr_Format (TreesignalcError, "%s received an uninitialized TreeSpace", function_name);
    return 0;
  }
  return 1;
}

/* Native object API. These functions receive owning TreeSpace Python objects
 * and reuse their already-parsed C topology_space structures. The C library
 * creates calculation-specific copies where existing algorithms mutate tree
 * or reconciliation state, leaving the Python-held topologies reusable. */
static PyObject *
treesignalc_from_tree_spaces (PyObject *self, PyObject *args)
{
  TreeSpaceObject *gene, *species;
  double *values = NULL;
  int length;

  if (!parse_tree_spaces (args, &gene, &species, "from_tree_spaces")) return NULL;
  length = genefam_module_treesignal_from_topology_spaces (gene->space, species->space, &values);
  return distance_vector_to_tuple (values, length);
}

static PyObject *
treesignalc_from_tree_spaces_rescale (PyObject *self, PyObject *args)
{
  TreeSpaceObject *gene, *species;
  double *values = NULL;
  int length;

  if (!parse_tree_spaces (args, &gene, &species, "from_tree_spaces_rescale")) return NULL;
  length = genefam_module_treesignal_from_topology_spaces_rescale (gene->space, species->space, &values);
  return distance_vector_to_tuple (values, length);
}

static PyObject *
treesignalc_from_tree_spaces_pvalue (PyObject *self, PyObject *args)
{
  TreeSpaceObject *gene, *species;
  double *values = NULL;
  int n_replicates = 1000, length;

  if (!PyArg_ParseTuple (args, "O!O!|i", &TreeSpaceType, &gene, &TreeSpaceType, &species, &n_replicates)) return NULL;
  if ((!gene->space) || (!species->space)) {
    PyErr_SetString (TreesignalcError, "from_tree_spaces_pvalue received an uninitialized TreeSpace");
    return NULL;
  }
  if (n_replicates < 10) n_replicates = 10;
  length = genefam_module_treesignal_from_topology_spaces_pvalue (
      gene->space, species->space, n_replicates, &values);
  return distance_vector_to_tuple (values, length);
}

/* Legacy string-based API. These names are retained for compatibility. Each
 * call parses its Newick strings into temporary C structures; new code should
 * construct TreeSpace objects once and use from_tree_spaces* above. */
static PyObject *
treesignalc_fromtrees_string (PyObject *self, PyObject *args)
{
  const char *gene_newick, *species_newick;
  double *values = NULL;
  int length;

  if (!PyArg_ParseTuple (args, "ss", &gene_newick, &species_newick)) return NULL;
  length = genefam_module_treesignal_fromtrees (gene_newick, species_newick, &values);
  return distance_vector_to_tuple (values, length);
}

static PyObject *
treesignalc_fromtrees_rescale_string (PyObject *self, PyObject *args)
{
  const char *gene_newick, *species_newick;
  double *values = NULL;
  int length;

  if (!PyArg_ParseTuple (args, "ss", &gene_newick, &species_newick)) return NULL;
  length = genefam_module_treesignal_fromtrees_rescale (gene_newick, species_newick, &values);
  return distance_vector_to_tuple (values, length);
}

static PyObject *
treesignalc_fromtrees_pvalue_string (PyObject *self, PyObject *args)
{
  const char *gene_newick, *species_newick;
  double *values = NULL;
  int n_replicates = 1000, length;

  if (!PyArg_ParseTuple (args, "ss|i", &gene_newick, &species_newick, &n_replicates)) return NULL;
  if (n_replicates < 10) n_replicates = 10;
  length = genefam_module_treesignal_fromtrees_pvalue (
      gene_newick, species_newick, n_replicates, &values);
  return distance_vector_to_tuple (values, length);
}

static PyObject *
treesignalc_randomise_trees_with_spr_string (PyObject *self, PyObject *args)
{
  char *output_trees = NULL;
  const char *species_newick;
  PyObject *result;
  int n_copies = 2, n_spr = 1;

  if (!PyArg_ParseTuple (args, "s|ii", &species_newick, &n_copies, &n_spr)) return NULL;
  output_trees = genefam_module_randomise_trees_with_spr (species_newick, n_copies, n_spr);
  if (!output_trees) {
    PyErr_SetString (TreesignalcError, "could not expand trees with SPR neighbours");
    return NULL;
  }
  result = PyUnicode_FromString (output_trees);
  free (output_trees);
  return result;
}

static PyObject *
treesignalc_generate_spr_trees_string (PyObject *self, PyObject *args)
{
  char *output_trees = NULL;
  PyObject *result;
  int n_leaves, n_iter, n_spr;

  if (!PyArg_ParseTuple (args, "iii", &n_leaves, &n_iter, &n_spr)) return NULL;
  output_trees = genefam_module_generate_spr_trees (n_leaves, n_iter, n_spr);
  if (!output_trees) {
    PyErr_SetString (TreesignalcError, "could not create a chain of SPR trees");
    return NULL;
  }
  result = PyUnicode_FromString (output_trees);
  free (output_trees);
  return result;
}

static PyMethodDef TreesignalcMethods[] = {
  {"from_tree_spaces", treesignalc_from_tree_spaces, METH_VARARGS,
   "Calculate distances using reusable native TreeSpace objects."},
  {"from_tree_spaces_rescale", treesignalc_from_tree_spaces_rescale, METH_VARARGS,
   "Calculate theoretically rescaled distances using reusable native TreeSpace objects."},
  {"from_tree_spaces_pvalue", treesignalc_from_tree_spaces_pvalue, METH_VARARGS,
   "Calculate p-values and empirically rescaled distances using reusable native TreeSpace objects."},
  {"fromtrees", treesignalc_fromtrees_string, METH_VARARGS,
   "Legacy string-based API: parse Newick strings and calculate distances."},
  {"fromtrees_rescale", treesignalc_fromtrees_rescale_string, METH_VARARGS,
   "Legacy string-based API: parse Newick strings and calculate rescaled distances."},
  {"fromtrees_pvalue", treesignalc_fromtrees_pvalue_string, METH_VARARGS,
   "Legacy string-based API: parse Newick strings and calculate p-values."},
  {"randomise_trees_with_spr", treesignalc_randomise_trees_with_spr_string, METH_VARARGS,
   "Legacy string-based API: return Newick strings for SPR neighbours."},
  {"generate_spr_trees", treesignalc_generate_spr_trees_string, METH_VARARGS,
   "Legacy string-based API: return a Newick string containing an SPR tree chain."},
  {NULL, NULL, 0, NULL}
};

PyDoc_STRVAR (
    treesignalc_doc,
    "Native tree structures and low-level distance functions for treesignal.\n\n"
    "Prefer TreeSpace and from_tree_spaces* for new code. The fromtrees* names "
    "are retained as legacy string-based compatibility functions.");

static struct PyModuleDef treesignalcmodule = {
  PyModuleDef_HEAD_INIT,
  "_treesignalc",
  treesignalc_doc,
  -1,
  TreesignalcMethods
};

PyMODINIT_FUNC
PyInit__treesignalc (void)
{
  PyObject *module;

  if (PyType_Ready (&TreeSpaceType) < 0) return NULL;
  module = PyModule_Create (&treesignalcmodule);
  if (!module) return NULL;

  TreesignalcError = PyErr_NewException ("treesignal._treesignalc.error", NULL, NULL);
  if (!TreesignalcError || PyModule_AddObject (module, "error", TreesignalcError) < 0) {
    Py_XDECREF (TreesignalcError);
    Py_DECREF (module);
    return NULL;
  }

  Py_INCREF (&TreeSpaceType);
  if (PyModule_AddObject (module, "TreeSpace", (PyObject*) &TreeSpaceType) < 0) {
    Py_DECREF (&TreeSpaceType);
    Py_DECREF (module);
    return NULL;
  }
  return module;
}
