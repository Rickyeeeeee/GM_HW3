#ifndef __LEARNPLY_H__
#define __LEARNPLY_H__

#include "ply.h"
#include <Eigen/Dense>
#include <vector>

/* forward declarations */
class Triangle;
class Vertex;
class Edge;
class Corner;

class Vertex {
public:
  int index;
  Eigen::Vector3d pos;
  std::vector<Triangle *> tris;
  std::vector<Corner *> corners;

  // Attributes
  Eigen::Vector3f color;
  Eigen::Vector3d normal;
  bool selected;

  // Curvature
  double area;
  double gaussCurvature;
  double meanCurvature;
  double maxPrincCurvature;
  double minPrincCurvature;
  // Frames
  Eigen::Matrix2d tensor;
  Eigen::Vector2d princDir2D[2];
  Eigen::Vector3d princDir3D[2];

public:
  Vertex(Eigen::Vector3d &p)
      : index(-1), pos(p), color(Eigen::Vector3f(0.58f, 0.72f, 0.75f)),
        normal(Eigen::Vector3d(0.0, 0.0, 1.0)), selected(false),
        area(0.0), gaussCurvature(0.0), meanCurvature(0.0), maxPrincCurvature(0.0), minPrincCurvature(0.0) {}
  Vertex(double xx, double yy, double zz)
      : index(-1), pos(Eigen::Vector3d(xx, yy, zz)), color(Eigen::Vector3f(0.58f, 0.72f, 0.75f)),
        normal(Eigen::Vector3d(0.0, 0.0, 1.0)), selected(false),
        area(0.0), gaussCurvature(0.0), meanCurvature(0.0), maxPrincCurvature(0.0), minPrincCurvature(0.0) {}

  int ntris() { return (int)tris.size(); }
};

class Edge {
public:
  int index;

  Vertex *verts[2];

  std::vector<Triangle *> tris;
  std::vector<Corner *> corners;

  // Attributes
  bool selected;
  double length;
  double dihedral_angle;

public:
  Edge() : index(-1), verts{NULL, NULL}, length(0.0), selected(false), dihedral_angle(0.0) {}
  int ntris() { return (int)tris.size(); }
};

class Triangle {
public:
  int index;

  Vertex *verts[3];
  Edge   *edges[3];
  Corner *corners[3];

  // Attributes
  bool selected;
  double area;
  Eigen::Vector3d normal;

public:
  Triangle()
      : index(-1), verts{NULL, NULL, NULL}, edges{NULL, NULL, NULL}, corners{NULL, NULL, NULL},
        area(0.0), normal(Eigen::Vector3d(0.0, 0.0, 1.0)), selected(false) {}
  int nverts() { return 3; }
};

class Corner {
public:
  Vertex *vertex;
  Edge *edge;
  Triangle *tri;

  Corner *prev;
  Corner *next;
  Corner *oppsite;

  // Attributes
  double interior_angle;

public:
  Corner()
      : vertex(NULL), edge(NULL), tri(NULL), prev(NULL), next(NULL), oppsite(NULL), interior_angle(0.0){}
};

class Polyhedron {
public:
  std::vector<Triangle *> tlist; /* list of triangles */
  std::vector<Vertex *> vlist;   /* list of vertices */
  std::vector<Edge *> elist;     /* list of edges */
  std::vector<Corner *> clist;   /* list of corners */

public:
  Polyhedron();
  Polyhedron(std::vector<Vertex *> &verts, std::vector<Triangle *> &tris, bool re_index = true);
  static Polyhedron *createPolyhedron(const std::string &filename);

  //
  int ntris() { return (int)tlist.size(); }
  int nverts() { return (int)vlist.size(); }
  int nedges() { return (int)elist.size(); }
  int ncorners() { return (int)clist.size(); }

  // initialization and finalization
  void initialize();
  void finalize();

  void recreate_corners();

private:
  void create_edge(Vertex *, Vertex *);
  void create_edges();
  void create_corners();
  void set_vertex_to_tri_ptrs();
};

#endif /* __LEARNPLY_H__ */
