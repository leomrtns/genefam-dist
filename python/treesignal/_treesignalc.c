#include <Python.h>
#include <stdlib.h>

#include "genefam_dist.h"
#include "newick_validation.h"

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
  if (!validate_newick_input (newick)) return -1;
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

static PyObject *
TreeSpace_get_taxa (TreeSpaceObject *self, void *closure)
{
  PyObject *result, *label;
  int i, n = self->space ? self->space->taxlabel->nstrings : 0;
  result = PyTuple_New (n);
  if (!result) return NULL;
  for (i = 0; i < n; i++) {
    label = PyUnicode_FromString (self->space->taxlabel->string[i]);
    if (!label) { Py_DECREF (result); return NULL; }
    PyTuple_SET_ITEM (result, i, label);
  }
  return result;
}

static PyObject *
TreeSpace_get_has_branch_lengths (TreeSpaceObject *self, void *closure)
{
  int i;
  if (!self->space) Py_RETURN_FALSE;
  for (i = 0; i < self->space->ndistinct; i++) if (!self->space->distinct[i]->blength) Py_RETURN_FALSE;
  Py_RETURN_TRUE;
}

/* Python owns these temporary strings. Avoid the fixed-size sprintf buffers in
 * the historical C serializer, and retain full double precision for lengths. */
static PyObject *
node_to_newick (topology tree, topol_node node)
{
  PyObject *left, *right, *base, *result;
  char length[40];
  if (node->internal) {
    left = node_to_newick (tree, node->left);
    if (!left) return NULL;
    right = node_to_newick (tree, node->right);
    if (!right) { Py_DECREF (left); return NULL; }
    base = PyUnicode_FromFormat ("(%U,%U)", left, right);
    Py_DECREF (left); Py_DECREF (right);
  }
  else base = PyUnicode_FromString (tree->taxlabel->string[node->id]);
  if (!base || !tree->blength || node == tree->root) return base;
  snprintf (length, sizeof (length), "%.17g", tree->blength[node->id]);
  result = PyUnicode_FromFormat ("%U:%s", base, length);
  Py_DECREF (base);
  return result;
}

static PyObject *
TreeSpace_get_newick (TreeSpaceObject *self, void *closure)
{
  PyObject *list, *item, *tree, *separator, *result;
  int i;
  if (!self->space) return PyUnicode_FromString ("");
  list = PyList_New (self->space->ndistinct);
  if (!list) return NULL;
  for (i = 0; i < self->space->ndistinct; i++) {
    tree = node_to_newick (self->space->distinct[i], self->space->distinct[i]->root);
    if (!tree) { Py_DECREF (list); return NULL; }
    item = PyUnicode_FromFormat ("%U;", tree);
    Py_DECREF (tree);
    if (!item) { Py_DECREF (list); return NULL; }
    PyList_SET_ITEM (list, i, item);
  }
  separator = PyUnicode_FromString ("\n");
  if (!separator) { Py_DECREF (list); return NULL; }
  result = PyUnicode_Join (separator, list);
  Py_DECREF (separator); Py_DECREF (list);
  return result;
}

static PyObject *
TreeSpace_patristic_matrix (TreeSpaceObject *self, PyObject *args, PyObject *kwargs)
{
  static char *keywords[] = {"index", NULL};
  int index = 0, i, j;
  distance_matrix dist;
  topology tree;
  PyObject *result, *row, *value;
  if (!PyArg_ParseTupleAndKeywords (args, kwargs, "|i:patristic_matrix", keywords, &index)) return NULL;
  if (!self->space || index < 0 || index >= self->space->ndistinct) {
    PyErr_SetString (PyExc_IndexError, "tree index out of range"); return NULL;
  }
  tree = self->space->distinct[index];
  if (!tree->blength) { PyErr_SetString (PyExc_ValueError, "patristic distances require complete branch lengths"); return NULL; }
  dist = new_distance_matrix_for_topology (tree->nleaves);
  fill_distance_matrix_from_topology (dist, tree, tree->blength, true);
  result = PyTuple_New (tree->nleaves);
  if (!result) { del_distance_matrix (dist); return NULL; }
  for (i = 0; i < tree->nleaves; i++) {
    row = PyTuple_New (tree->nleaves);
    if (!row) { Py_DECREF (result); del_distance_matrix (dist); return NULL; }
    PyTuple_SET_ITEM (result, i, row);
    for (j = 0; j < tree->nleaves; j++) {
      if (!isfinite (i < j ? dist->d[i][j] : dist->d[j][i])) {
        Py_DECREF (result); del_distance_matrix (dist);
        PyErr_SetString (PyExc_ValueError, "patristic distances exceed the finite numeric range"); return NULL;
      }
      value = PyFloat_FromDouble (i < j ? dist->d[i][j] : dist->d[j][i]);
      if (!value) { Py_DECREF (result); del_distance_matrix (dist); return NULL; }
      PyTuple_SET_ITEM (row, j, value);
    }
  }
  del_distance_matrix (dist);
  return result;
}

static PyMethodDef TreeSpace_methods[] = {
  {"patristic_matrix", (PyCFunction) TreeSpace_patristic_matrix, METH_VARARGS | METH_KEYWORDS,
   "Return the raw patristic matrix in taxa order for the selected reference."},
  {NULL, NULL, 0, NULL}
};

static PyGetSetDef TreeSpace_getset[] = {
  {"ntrees", (getter) TreeSpace_get_ntrees, NULL, "Number of input trees.", NULL},
  {"ndistinct", (getter) TreeSpace_get_ndistinct, NULL, "Number of retained trees (lengthless topologies are deduplicated).", NULL},
  {"rooted", (getter) TreeSpace_get_rooted, NULL, "Whether root location distinguishes topologies.", NULL},
  {"taxa", (getter) TreeSpace_get_taxa, NULL, "Leaf labels in native matrix order.", NULL},
  {"has_branch_lengths", (getter) TreeSpace_get_has_branch_lengths, NULL, "Whether every retained tree has lengths.", NULL},
  {"newick", (getter) TreeSpace_get_newick, NULL, "Retained trees serialized with their raw lengths.", NULL},
  {NULL, NULL, NULL, NULL, NULL}
};

PyDoc_STRVAR (
    TreeSpace_doc,
    "TreeSpace(newick, rooted=False)\n--\n\n"
    "Owning native representation of one or more trees. Distance functions use "
    "the first gene tree and all retained species trees. Length-bearing trees "
    "retain input order; only lengthless topologies are deduplicated.\n\n"
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
  .tp_methods = TreeSpace_methods,
  .tp_init = (initproc) TreeSpace_init,
  .tp_new = TreeSpace_new,
};

static int
check_legacy_mapping (topology_space gene, topology_space species)
{
  int i, j;
  for (i = 0; i < gene->taxlabel->nstrings; i++) {
    for (j = 0; j < species->taxlabel->nstrings; j++)
      if (strcasestr (gene->taxlabel->string[i], species->taxlabel->string[j])) break;
    if (j == species->taxlabel->nstrings) {
      PyErr_Format (PyExc_ValueError, "gene leaf %s does not contain a reference species label", gene->taxlabel->string[i]);
      return 0;
    }
  }
  return 1;
}

static int
parse_tree_spaces (PyObject *args, TreeSpaceObject **gene, TreeSpaceObject **species, const char *function_name)
{
  if (!PyArg_ParseTuple (args, "O!O!", &TreeSpaceType, gene, &TreeSpaceType, species)) return 0;
  if ((!(*gene)->space) || (!(*species)->space)) {
    PyErr_Format (TreesignalcError, "%s received an uninitialized TreeSpace", function_name);
    return 0;
  }
  return check_legacy_mapping ((*gene)->space, (*species)->space);
}

static PyObject *
treesignalc_from_tree_spaces_lengths (PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *keywords[] = {"gene", "species", "normalize", NULL};
  TreeSpaceObject *gene, *species;
  double *values = NULL;
  int normalize = 0, length;
  if (!PyArg_ParseTupleAndKeywords (args, kwargs, "O!O!|p:from_tree_spaces_lengths", keywords,
                                   &TreeSpaceType, &gene, &TreeSpaceType, &species, &normalize)) return NULL;
  if (!gene->space || !species->space) {
    PyErr_SetString (PyExc_ValueError, "received an uninitialized TreeSpace"); return NULL;
  }
  length = genefam_module_treesignal_lengths (gene->space, species->space, normalize != 0, &values);
  if (length < 0) {
    PyErr_SetString (PyExc_ValueError, "derived branch-length distances exceed the finite numeric range"); return NULL;
  }
  if (!length) {
    PyErr_SetString (PyExc_ValueError, "length distances require complete branch lengths and at least two common leaf labels");
    return NULL;
  }
  for (int i = 0; i < length; i++) if (!isfinite (values[i])) {
    free (values);
    PyErr_SetString (PyExc_ValueError, "derived branch-length distances exceed the finite numeric range");
    return NULL;
  }
  return distance_vector_to_tuple (values, length);
}

static PyObject *
treesignalc_generate_tree_space (PyObject *self, PyObject *args, PyObject *kwargs)
{
  static char *keywords[] = {"n_species", "chain_size", "n_spr", "seed", NULL};
  int n_species, chain_size, n_spr = 1;
  unsigned long long seed = 0;
  PyObject *seed_object = NULL;
  TreeSpaceObject *result;
  if (!PyArg_ParseTupleAndKeywords (args, kwargs, "ii|iO:generate_tree_space", keywords,
                                   &n_species, &chain_size, &n_spr, &seed_object)) return NULL;
  if (seed_object) {
    seed = PyLong_AsUnsignedLongLong (seed_object);
    if (PyErr_Occurred ()) return NULL;
  }
  if (n_species < 4 || chain_size < 1 || n_spr < 1) {
    PyErr_SetString (PyExc_ValueError, "generation requires at least 4 species, 1 tree and 1 SPR move"); return NULL;
  }
  result = (TreeSpaceObject*) TreeSpace_new (&TreeSpaceType, NULL, NULL);
  if (!result) return NULL;
  result->space = genefam_module_generate_tree_space (n_species, chain_size, n_spr, seed);
  result->rooted = 1;
  return (PyObject*) result;
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
  if (!check_legacy_mapping (gene->space, species->space)) return NULL;
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
  if (!validate_newick_input (gene_newick) || !validate_newick_input (species_newick)) return NULL;
  topology_space gene = genefam_module_parse_newick_trees (gene_newick, false);
  topology_space species = genefam_module_parse_newick_trees (species_newick, true);
  if (!check_legacy_mapping (gene, species)) { del_topology_space (gene); del_topology_space (species); return NULL; }
  length = genefam_module_treesignal_from_topology_spaces (gene, species, &values);
  del_topology_space (gene); del_topology_space (species);
  return distance_vector_to_tuple (values, length);
}

static PyObject *
treesignalc_fromtrees_rescale_string (PyObject *self, PyObject *args)
{
  const char *gene_newick, *species_newick;
  double *values = NULL;
  int length;

  if (!PyArg_ParseTuple (args, "ss", &gene_newick, &species_newick)) return NULL;
  if (!validate_newick_input (gene_newick) || !validate_newick_input (species_newick)) return NULL;
  topology_space gene = genefam_module_parse_newick_trees (gene_newick, false);
  topology_space species = genefam_module_parse_newick_trees (species_newick, true);
  if (!check_legacy_mapping (gene, species)) { del_topology_space (gene); del_topology_space (species); return NULL; }
  length = genefam_module_treesignal_from_topology_spaces_rescale (gene, species, &values);
  del_topology_space (gene); del_topology_space (species);
  return distance_vector_to_tuple (values, length);
}

static PyObject *
treesignalc_fromtrees_pvalue_string (PyObject *self, PyObject *args)
{
  const char *gene_newick, *species_newick;
  double *values = NULL;
  int n_replicates = 1000, length;

  if (!PyArg_ParseTuple (args, "ss|i", &gene_newick, &species_newick, &n_replicates)) return NULL;
  if (!validate_newick_input (gene_newick) || !validate_newick_input (species_newick)) return NULL;
  if (n_replicates < 10) n_replicates = 10;
  topology_space gene = genefam_module_parse_newick_trees (gene_newick, false);
  topology_space species = genefam_module_parse_newick_trees (species_newick, true);
  if (!check_legacy_mapping (gene, species)) { del_topology_space (gene); del_topology_space (species); return NULL; }
  length = genefam_module_treesignal_from_topology_spaces_pvalue (gene, species, n_replicates, &values);
  del_topology_space (gene); del_topology_space (species);
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
  if (!validate_newick_input (species_newick)) return NULL;
  if (n_copies < 0 || n_spr < 0) { PyErr_SetString (PyExc_ValueError, "copies and SPR moves cannot be negative"); return NULL; }
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
  if (n_spr < 0) { PyErr_SetString (PyExc_ValueError, "SPR moves cannot be negative"); return NULL; }
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
  {"from_tree_spaces_lengths", (PyCFunction) treesignalc_from_tree_spaces_lengths, METH_VARARGS | METH_KEYWORDS,
   "Patristic RMS, weighted RF and KF on common taxa; normalize by matrix mean and induced total length respectively."},
  {"generate_tree_space", (PyCFunction) treesignalc_generate_tree_space, METH_VARARGS | METH_KEYWORDS,
   "Generate a native SPR chain with unit unrooted edge lengths and T1...Tn labels."},
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

  TreesignalcError = PyErr_NewException ("treesignal._treesignalc.error", PyExc_ValueError, NULL);
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
