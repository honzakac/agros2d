// This file is part of Hermes2D
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, see <http://www.gnu.prg/licenses/>.

#include <string.h>
#include "mesh.h"
#include "h2d_reader.h"
#include <map>
#include "hash.h"
#include "mesh_parser.h"

extern unsigned g_mesh_seq;

H2DReader::H2DReader()
{
}

H2DReader::~H2DReader()
{
}

//// load //////////////////////////////////////////////////////////////////////////////////////////

/*************** OLD OLD OLD ***************************************************************/

#define eof_error error("Premature end of file")

char* get_line(FILE* f)
{
  static char line[1000];

  // read one line, skipping empty ones and those starting with '*' or '#'
  while (fgets(line, 1000-1, f) != NULL)
  {
    char* p = line;
    while (*p && (unsigned) *p <= ' ') p++;

    if (*p && *p != '*' && *p != '#') return line;
  }

  return NULL;
}


Nurbs* H2DReader::load_nurbs_old(Mesh *mesh, FILE* f, Node** en, int &p1, int &p2)
{
  int i;
  char* line;
  Nurbs* nurbs = new Nurbs;

  // read the end point indices
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d %d", &p1, &p2) != 2) error("Error reading curved boundary data (end point indices)");
  *en = mesh->peek_edge_node(p1, p2);
  if (*en == NULL) error("Error reading curved boundary data (edge %d-%d does not exist)", p1, p2);

  // degree of curved edge
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &(nurbs->degree)) != 1)
    error("Error reading curved boundary data for edge %d-%d (degree)", p1, p2);

  // create a circular arc if degree == 0
  bool circle = (nurbs->degree == 0);
  if (circle) nurbs->degree = 2;
  nurbs->arc = circle;

  // load control points of curved edge
  int inner = 1, outer;
  if (!circle)
  {
    if ((line = get_line(f)) == NULL) eof_error;
    if (sscanf(line, "%d", &inner) != 1)
      error("Error reading curved boundary data for edge %d-%d (# of control points)", p1, p2);
  }
  nurbs->np = inner + 2;

  // edge endpoints are also control points, with weight 1.0
  nurbs->pt = new double3[nurbs->np];
  nurbs->pt[0][0] = mesh->nodes[p1].x;
  nurbs->pt[0][1] = mesh->nodes[p1].y;
  nurbs->pt[0][2] = 1.0;
  nurbs->pt[inner+1][0] = mesh->nodes[p2].x;
  nurbs->pt[inner+1][1] = mesh->nodes[p2].y;
  nurbs->pt[inner+1][2] = 1.0;

  if (!circle)
  {
    // load inner control points
    for (i = 1; i <= inner; i++)
    {
      if ((line = get_line(f)) == NULL) eof_error;
      if (sscanf(line, "%lf %lf %lf", &(nurbs->pt[i][0]), &(nurbs->pt[i][1]), &(nurbs->pt[i][2])) != 3)
        error("Error reading curved boundary data for edge %d-%d (control points)", p1, p2);
    }
  }
  else
  {
    // load the arc angle
    if ((line = get_line(f)) == NULL) eof_error;
    double angle;
    if (sscanf(line, "%lf", &angle) != 1)
      error("Error reading curved boundary data for edge %d-%d (arc angle)", p1, p2);
    nurbs->angle = angle;
    angle = (180.0 - angle) / 180.0 * M_PI;

    // generate one control point
    double x = 1.0 / tan(angle * 0.5);
    nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
    nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
    nurbs->pt[1][2] = cos((M_PI - angle) * 0.5);
  }

  // load the number of knot vector points
  inner = 0;
  if (!circle)
  {
    if ((line = get_line(f)) == NULL) eof_error;
    if (sscanf(line, "%d", &(inner)) != 1)
      error("Error reading curved boundary data for edge %d-%d (# of knot vector points)", p1, p2);
  }

  nurbs->nk = nurbs->degree + nurbs->np + 1;
  outer = nurbs->nk - inner;
  if ((outer & 1) == 1)
    error("Error reading curved boundary data for edge %d-%d (wrong number of knot points)", p1, p2);

  // knot vector is completed by 0.0 on the left and by 1.0 on the right
  nurbs->kv = new double[nurbs->nk];
  for (i = 0; i < outer/2; i++)
    nurbs->kv[i] = 0.0;
  for (i = outer/2; i < inner + outer/2; i++)
  {
    if ((line = get_line(f)) == NULL) eof_error;
    if (sscanf(line, "%lf", &(nurbs->kv[i])) != 1)
      error("Error reading curved boundary data for edge %d-%d (knot points)", p1, p2);
  }
  for (i = outer/2 + inner; i < nurbs->nk; i++)
    nurbs->kv[i] = 1.0;

  nurbs->ref = 0;
  return nurbs;
}

void H2DReader::load_old(const char* filename, Mesh *mesh)
{
  // open the mesh file
  FILE* f = fopen(filename, "r");
  if (f == NULL) error("Could not open the mesh file %s", filename);
  this->load_stream(f, mesh);
}

void H2DReader::load_str(char* mesh_str, Mesh *mesh)
{
  // open the mesh file
  FILE* f = fmemopen(mesh_str, strlen(mesh_str), "r");
  if (f == NULL) error("Could not create the read buffer");
  this->load_stream(f, mesh);
}

/*
   Loads the mesh from a stream.
*/
void H2DReader::load_stream(FILE *f, Mesh *mesh)
{
  int i, j, k, n, maj, min;
  char* line;

  // check file version
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d %d", &maj, &min) != 2) error("Could not read file version");
  if (maj > 1) error("Unsupported file version");

  // read the number of vertices
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &n) != 1) error("Could not read the number of vertices");

  // free all current data
  mesh->free();

  // create a hash table large enough
  int size = HashTable::H2D_DEFAULT_HASH_SIZE;
  while (size < 8*n) size *= 2;
  mesh->init(size);

  // load vertices: create top-level vertex nodes
  for (i = 0; i < n; i++)
  {
    Node* node = mesh->nodes.add();
    assert(node->id == i);
    node->ref = TOP_LEVEL_REF;
    node->type = H2D_TYPE_VERTEX;
    node->bnd = 0;
    node->p1 = node->p2 = -1;
    node->next_hash = NULL;

    if ((line = get_line(f)) == NULL) eof_error;
    if (sscanf(line, "%lf %lf", &node->x, &node->y) != 2) error("Error reading vertex data");
  }
  mesh->ntopvert = n;

  // read the number of elements
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &n) != 1) error("Could not read the number of elements");

  // load elements
  for (i = 0; i < n; i++)
  {
    if ((line = get_line(f)) == NULL) eof_error;

    int ret, idx[5];
    if ((ret = sscanf(line, "%d %d %d %d %d", idx, idx+1, idx+2, idx+3, idx+4)) != 4 && ret != 5)
      error("Error reading elements");

    for (j = 0; j < ret-1; j++)
      if (idx[j] < 0 || idx[j] >= mesh->ntopvert)
        error("Error reading elements: node %d does not exist", idx[j]);

    Node *v0 = &mesh->nodes[idx[0]], *v1 = &mesh->nodes[idx[1]], *v2 = &mesh->nodes[idx[2]];
    if (ret == 4)
    {
      check_triangle(i, v0, v1, v2);
      mesh->create_triangle(idx[3], v0, v1, v2, NULL);
    }
    else
    {
      Node *v3 = &mesh->nodes[idx[3]];
      check_quad(i, v0, v1, v2, v3);
      mesh->create_quad(idx[4], v0, v1, v2, v3, NULL);
    }
  }
  mesh->nbase = mesh->nactive = n;

  // read the number of boundary data
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &n) != 1) error("Could not read the number of boundary markers\n");

  // load boundary data
  Node* en;
  for (i = 0; i < n; i++)
  {
    if ((line = get_line(f)) == NULL) eof_error;

    int v1, v2, marker;
    if (sscanf(line, "%d %d %d", &v1, &v2, &marker) != 3) error("Error reading boundary marker data");

    en = mesh->peek_edge_node(v1, v2);
    if (en == NULL) error("Boundary data error (edge %d-%d does not exist)", v1, v2);
    en->marker = marker;

    if (marker > 0)
    {
      mesh->nodes[v1].bnd = 1;
      mesh->nodes[v2].bnd = 1;
      en->bnd = 1;
    }
  }

  // check that all boundary edges have a marker assigned
  for_all_edge_nodes(en, mesh)
    if (en->ref < 2 && en->marker == 0)
      warn("Boundary edge node does not have a boundary marker");

  // read the number of curved edges
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &n) != 1) error("Could not read the number of curved edges");

  // load curved edges
  for (i = 0; i < n; i++)
  {
    // load the control points, knot vector, etc.
    Node* en;
    int p1, p2;
    Nurbs* nurbs = load_nurbs_old(mesh, f, &en, p1, p2);

    // assign the nurbs to the elements sharing the edge node
    for (k = 0; k < 2; k++)
    {
      Element* e = en->elem[k];
      if (e == NULL) continue;

      if (e->cm == NULL)
      {
        e->cm = new CurvMap;
        memset(e->cm, 0, sizeof(CurvMap));
        e->cm->toplevel = 1;
        e->cm->order = 4;
      }

      int idx = -1;
      for (j = 0; j < e->nvert; j++)
        if (e->en[j] == en) { idx = j; break; }
      assert(idx >= 0);

      if (e->vn[idx]->id == p1)
      {
        e->cm->nurbs[idx] = nurbs;
        nurbs->ref++;
      }
      else
      {
        Nurbs* nurbs_rev = mesh->reverse_nurbs(nurbs);
        e->cm->nurbs[idx] = nurbs_rev;
        nurbs_rev->ref++;
      }
    }
    if (!nurbs->ref) delete nurbs;
  }

  // read the number of initial refinements
  //if () eof_error;
  //if (sscanf(line, "%d", &n) != 1) error("could not read the number of initial refinements");

  if ((line = get_line(f)) == NULL ||
      sscanf(line, "%d", &n) != 1)
  {
    warn("Could not read the number of initial refinements");
  }
  else
  {
    // perform initial refinements
    for (i = 0; i < n; i++)
    {
      if ((line = get_line(f)) == NULL) eof_error;
      int id, ref;
      if (sscanf(line, "%d %d", &id, &ref) != 2)
        error("Error reading initial refinement data");
      mesh->refine_element(id, ref);
    }
  }
  mesh->ninitial = mesh->elements.get_num_items();

  // update refmap coeffs of curvilinear elements
  Element* e;
  for_all_elements(e, mesh)
    if (e->cm != NULL)
      e->cm->update_refmap_coeffs(e);

  fclose(f);
  mesh->seq = g_mesh_seq++;
}

/*************** OLD OLD OLD ***************************************************************/


//// load_nurbs ////////////////////////////////////////////////////////////////////////////////////

Nurbs* H2DReader::load_nurbs(Mesh *mesh, MItem* curve, int id, Node** en, int &p1, int &p2)
{
  int i;
  Nurbs* nurbs = new Nurbs;

  if (curve == NULL || curve->n < 0 || curve->n != 3 && curve->n != 5)
    error("Invalid curve #%d.", id);
  bool circle = (curve->n == 3);
  nurbs->arc = circle;

  // read the end point indices
  MItem* edge = curve->list;
  if (edge->n >= 0 || !H2D_IS_INT(edge->val))
    error("Curve #%d: invalid edge definition.", id);
  p1 = (int) edge->val;
  edge = edge->next;

  if (edge->n >= 0 || !H2D_IS_INT(edge->val))
    error("Curve #%d: invalid edge definition.", id);
  p2 = (int) edge->val;
  edge = edge->next;

  *en = mesh->peek_edge_node(p1, p2);
  if (*en == NULL)
    error("Curve #%d: edge %d-%d does not exist.", id, p1, p2);

  // degree of curved edge
  MItem* deg = edge;
  nurbs->degree = 2;
  if (!circle)
  {
    if (deg == NULL || deg->n >= 0 || !H2D_IS_INT(deg->val) || deg->val < 0 || deg->val == 1)
      error("Curve #%d: invalid degee.", id);
    nurbs->degree = (int) deg->val;
  }

  // get the number of control points
  MItem* pts = deg->next;
  int inner = 1, outer;
  if (!circle)
  {
    if (pts == NULL || pts->n < 0)
      error("Curve #%d: control points not defined.", id);
    inner = pts->n;
  }
  nurbs->np = inner + 2;

  // edge endpoints are also control points, with weight 1.0
  nurbs->pt = new double3[nurbs->np];
  nurbs->pt[0][0] = mesh->nodes[p1].x;
  nurbs->pt[0][1] = mesh->nodes[p1].y;
  nurbs->pt[0][2] = 1.0;
  nurbs->pt[inner+1][0] = mesh->nodes[p2].x;
  nurbs->pt[inner+1][1] = mesh->nodes[p2].y;
  nurbs->pt[inner+1][2] = 1.0;

  if (!circle)
  {
    // read inner control points
    MItem* it = pts->list;
    for (i = 1; i <= inner; i++, it = it->next)
    {
      if (!mesh_parser_get_doubles(it, 3, &(nurbs->pt[i][0]), &(nurbs->pt[i][1]), &(nurbs->pt[i][2])))
        error("Curve #%d: invalid control point #%d.", id, i-1);
    }
  }
  else
  {
    // read the arc angle
    MItem* angle = deg;
    if (angle == NULL || angle->n >= 0)
      error("Curve #%d: invalid arc angle.", id);
    double a = (180.0 - angle->val) / 180.0 * M_PI;
    nurbs->angle = angle->val;

    // generate one control point
    double x = 1.0 / tan(a * 0.5);
    nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
    nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
    nurbs->pt[1][2] = cos((M_PI - a) * 0.5);
  }

  // get the number of knot vector points
  inner = 0;
  MItem* knot=pts;
  if (!circle && knot != NULL)
  {
    knot = pts->next;
    if (knot->n < 0) error("Curve #%d: invalid knot vector.", id);
    inner = knot->n;
  }

  nurbs->nk = nurbs->degree + nurbs->np + 1;
  outer = nurbs->nk - inner;
  if ((outer & 1) == 1)
    error("Curve #%d: incorrect number of knot points.", id);

  // knot vector is completed by 0.0 on the left and by 1.0 on the right
  nurbs->kv = new double[nurbs->nk];
  for (i = 0; i < outer/2; i++)
    nurbs->kv[i] = 0.0;
  if (inner) {
    MItem* it = knot->list;
    for (i = outer/2; i < inner + outer/2; i++, it = it->next) {
      if (it->n >= 0) error("Curve #%d: invalid knot vector item.", id);
      nurbs->kv[i] = it->val;
    }
  }
  for (i = outer/2 + inner; i < nurbs->nk; i++)
    nurbs->kv[i] = 1.0;

  nurbs->ref = 0;
  return nurbs;
}


bool H2DReader::load(const char *filename, Mesh *mesh)
{
  int i, j, k, n;
  Node* en;
  bool debug = false;

  // open the mesh file
  FILE* f = fopen(filename, "r");
  if (f == NULL) error("Could not open the mesh file %s", filename);

  mesh->free();

  // parse the file
  mesh_parser_init(f, filename);
  mesh_parser_run(debug);
  fclose(f);

  //// vertices ////////////////////////////////////////////////////////////////

  MSymbol* sym = mesh_parser_find_symbol("vertices");
  if (sym == NULL) error("File %s: 'vertices' not found.", filename);
  n = sym->data->n;
  if (n < 0) error("File %s: 'vertices' must be a list.", filename);
  if (n < 2) error("File %s: invalid number of vertices.", filename);

  // create a hash table large enough
  int size = HashTable::H2D_DEFAULT_HASH_SIZE;
  while (size < 8*n) size *= 2;
  mesh->init(size);

  // create top-level vertex nodes
  MItem* pair = sym->data->list;
  for (i = 0; i < n; i++, pair = pair->next)
  {
    Node* node = mesh->nodes.add();
    assert(node->id == i);
    node->ref = TOP_LEVEL_REF;
    node->type = H2D_TYPE_VERTEX;
    node->bnd = 0;
    node->p1 = node->p2 = -1;
    node->next_hash = NULL;

    if (!mesh_parser_get_doubles(pair, 2, &node->x, &node->y))
      error("File %s: invalid vertex #%d.", filename, i);
  }
  mesh->ntopvert = n;

  //// elements ////////////////////////////////////////////////////////////////

  sym = mesh_parser_find_symbol("elements");
  if (sym == NULL) error("File %s: 'elements' not found.", filename);
  n = sym->data->n;
  if (n < 0) error("File %s: 'elements' must be a list.", filename);
  if (n < 1) error("File %s: no elements defined.", filename);

  // create elements
  MItem* elem = sym->data->list;
  mesh->nactive = 0;
  for (i = 0; i < n; i++, elem = elem->next)
  {
    // read and check vertex indices
    int nv = elem->n, idx[5];
    if (!nv) { mesh->elements.skip_slot(); continue; }
    if (nv < 4 || nv > 5)
      error("File %s: element #%d: wrong number of vertex indices.", filename, i);
    if (!mesh_parser_get_ints(elem, nv, &idx[0], &idx[1], &idx[2], &idx[3], &idx[4]))
      error("File %s: invalid definition of element #%d.", filename, i);
    for (j = 0; j < nv-1; j++)
      if (idx[j] < 0 || idx[j] >= mesh->ntopvert)
        error("File %s: error creating element #%d: vertex #%d does not exist.", filename, i, idx[j]);

    // create triangle/quad
    Node *v0 = &mesh->nodes[idx[0]], *v1 = &mesh->nodes[idx[1]], *v2 = &mesh->nodes[idx[2]];
    if (nv == 4)
    {
      check_triangle(i, v0, v1, v2);
      mesh->create_triangle(idx[3], v0, v1, v2, NULL);
    }
    else
    {
      Node *v3 = &mesh->nodes[idx[3]];
      check_quad(i, v0, v1, v2, v3);
      mesh->create_quad(idx[4], v0, v1, v2, v3, NULL);
    }
    mesh->nactive++;
  }
  mesh->nbase = n;

  //// boundaries //////////////////////////////////////////////////////////////

  sym = mesh_parser_find_symbol("boundaries");
  if (sym != NULL)
  {
    n = sym->data->n;
    if (n < 0) error("File %s: 'boundaries' must be a list.", filename);

    // read boundary data
    MItem* triple = sym->data->list;
    for (i = 0; i < n; i++, triple = triple->next)
    {
      int v1, v2, marker;
      if (!mesh_parser_get_ints(triple, 3, &v1, &v2, &marker))
        error("File %s: invalid boundary data #%d.", filename, i);

      en = mesh->peek_edge_node(v1, v2);
      if (en == NULL)
        error("File %s: boundary data #%d: edge %d-%d does not exist", filename, i, v1, v2);
      en->marker = marker;

      if (marker > 0)
      {
        mesh->nodes[v1].bnd = 1;
        mesh->nodes[v2].bnd = 1;
        en->bnd = 1;
      }
    }
  }

  // check that all boundary edges have a marker assigned
  for_all_edge_nodes(en, mesh)
    if (en->ref < 2 && en->marker == 0)
      warn("Boundary edge node does not have a boundary marker");

  //// curves //////////////////////////////////////////////////////////////////

  sym = mesh_parser_find_symbol("curves");
  if (sym != NULL)
  {
    n = sym->data->n;
    if (n < 0) error("File %s: 'curves' must be a list.", filename);

    // load curved edges
    MItem* curve = sym->data->list;
    for (i = 0; i < n; i++, curve = curve->next)
    {
      // load the control points, knot vector, etc.
      Node* en;
      int p1, p2;
      Nurbs* nurbs = load_nurbs(mesh, curve, i, &en, p1, p2);

      // assign the nurbs to the elements sharing the edge node
      for (k = 0; k < 2; k++)
      {
        Element* e = en->elem[k];
        if (e == NULL) continue;

        if (e->cm == NULL)
        {
          e->cm = new CurvMap;
          memset(e->cm, 0, sizeof(CurvMap));
          e->cm->toplevel = 1;
          e->cm->order = 4;
        }

        int idx = -1;
        for (j = 0; j < e->nvert; j++)
          if (e->en[j] == en) { idx = j; break; }
        assert(idx >= 0);

        if (e->vn[idx]->id == p1)
        {
          e->cm->nurbs[idx] = nurbs;
          nurbs->ref++;
        }
        else
        {
          Nurbs* nurbs_rev = mesh->reverse_nurbs(nurbs);
          e->cm->nurbs[idx] = nurbs_rev;
          nurbs_rev->ref++;
        }
      }
      if (!nurbs->ref) delete nurbs;
    }
  }

  // update refmap coeffs of curvilinear elements
  Element* e;
  for_all_elements(e, mesh)
    if (e->cm != NULL)
      e->cm->update_refmap_coeffs(e);

  //// refinements /////////////////////////////////////////////////////////////

  sym = mesh_parser_find_symbol("refinements");
  if (sym != NULL)
  {
    n = sym->data->n;
    if (n < 0) error("File %s: 'refinements' must be a list.", filename);

    // perform initial refinements
    MItem* pair = sym->data->list;
    for (i = 0; i < n; i++, pair = pair->next)
    {
      int id, ref;
      if (!mesh_parser_get_ints(pair, 2, &id, &ref))
        error("File %s: invalid refinement #%d.", filename, i);
      mesh->refine_element(id, ref);
    }
  }
  mesh->ninitial = mesh->elements.get_num_items();

  mesh_parser_free();
  mesh->seq = g_mesh_seq++;

  return true;
}

//// save ////////////////////////////////////////////////////////////////////////////////////

void H2DReader::save_refinements(Mesh *mesh, FILE* f, Element* e, int id, bool& first)
{
  if (e->active) return;
  fprintf(f, first ? "refinements =\n{\n" : ",\n"); first = false;
  if (e->bsplit())
  {
    fprintf(f, "  { %d, 0 }", id);
    int sid = mesh->seq; mesh->seq += 4;
    for (int i = 0; i < 4; i++)
      save_refinements(mesh, f, e->sons[i], sid+i, first);
  }
  else if (e->hsplit())
  {
    fprintf(f, "  { %d, 1 }", id);
    int sid = mesh->seq; mesh->seq += 2;
    save_refinements(mesh, f, e->sons[0], sid, first);
    save_refinements(mesh, f, e->sons[1], sid+1, first);
  }
  else
  {
    fprintf(f, "  { %d, 2 }", id);
    int sid = mesh->seq; mesh->seq += 2;
    save_refinements(mesh, f, e->sons[2], sid, first);
    save_refinements(mesh, f, e->sons[3], sid+1, first);
  }
}


void H2DReader::save_nurbs(Mesh *mesh, FILE* f, int p1, int p2, Nurbs* nurbs)
{
  if (nurbs->arc)
  {
    fprintf(f, "  { %d, %d, %.16g }", p1, p2, nurbs->angle);
  }
  else
  {
    int inner = nurbs->np - 2;
    int outer = nurbs->nk - inner;
    fprintf(f, "  { %d, %d, %d, { ", p1, p2, nurbs->degree);
    for (int i = 1; i < nurbs->np-1; i++)
      fprintf(f, "{ %.16g, %.16g, %.16g }%s ",
                 nurbs->pt[i][0], nurbs->pt[i][1], nurbs->pt[i][2],
                 i < nurbs->np-2 ? "," : "");

    fprintf(f, "}, { ", nurbs->nk - 2*(nurbs->degree+1));
    int max = nurbs->nk - (nurbs->degree+1);
    for (int i = nurbs->degree+1; i < max; i++)
      fprintf(f, "%.16g%s", nurbs->kv[i], i < max-1 ? "," : "");
    fprintf(f, "} }");
  }
}


static bool is_twin_nurbs(Element* e, int i)
{
  // on internal edges, where there are two Nurbs', we only save one of them
  return e->cm->nurbs[i]->twin && e->en[i]->ref == 2;
}

bool H2DReader::save(const char* filename, Mesh *mesh)
{
  int i, mrk;
  Element* e;

  // open output file
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not create mesh file.");
  //fprintf(f, "# hermes2d saved mesh\n\n");

  // save vertices
  fprintf(f, "vertices =\n{\n");
  for (i = 0; i < mesh->ntopvert; i++)
    fprintf(f, "  { %.16g, %.16g }%s\n", mesh->nodes[i].x, mesh->nodes[i].y, (i < mesh->ntopvert-1 ? "," : ""));

  // save elements
  fprintf(f, "}\n\nelements =\n{");
  bool first = true;
  for (i = 0; i < mesh->get_num_base_elements(); i++)
  {
    const char* nl = first ? "\n" : ",\n";  first = false;
    e = mesh->get_element_fast(i);
    if (!e->used)
      fprintf(f, "%s  { }", nl);
    else if (e->is_triangle())
      fprintf(f, "%s  { %d, %d, %d, %d }", nl, e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, e->marker);
    else
      fprintf(f, "%s  { %d, %d, %d, %d, %d }", nl, e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, e->vn[3]->id, e->marker);
  }

  // save boundary markers
  fprintf(f, "\n}\n\nboundaries =\n{");
  first = true;
  for_all_base_elements(e, mesh)
    for (i = 0; i < e->nvert; i++)
      if ((mrk = mesh->get_base_edge_node(e, i)->marker)) {
        const char* nl = first ? "\n" : ",\n";  first = false;
        fprintf(f, "%s  { %d, %d, %d }", nl, e->vn[i]->id, e->vn[e->next_vert(i)]->id, mrk);
      }
  fprintf(f, "\n}\n\n");

  // save curved edges
  first = true;
  for_all_base_elements(e, mesh)
    if (e->is_curved())
      for (i = 0; i < e->nvert; i++)
        if (e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i)) {
          fprintf(f, first ? "curves =\n{\n" : ",\n");  first = false;
          save_nurbs(mesh, f, e->vn[i]->id, e->vn[e->next_vert(i)]->id, e->cm->nurbs[i]);
        }
  if (!first) fprintf(f, "\n}\n\n");

  // save refinements
  unsigned temp = mesh->seq;
  mesh->seq = mesh->nbase;
  first = true;
  for_all_base_elements(e, mesh)
    save_refinements(mesh, f, e, e->id, first);
  if (!first) fprintf(f, "\n}\n\n");

  mesh->seq = temp;
  fclose(f);

  return true;
}
