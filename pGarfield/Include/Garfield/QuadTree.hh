#ifndef G_QUAD_TREE_H
#define G_QUAD_TREE_H

#include <cstddef>
#include <tuple>
#include <vector>

namespace Garfield {

/// Quadtree search.

class QuadTree {
 public:
  /// Constructor
  QuadTree(const double x0, const double y0, const double hx, const double hy);

  /// Destructor
  ~QuadTree();

  /// Insert a mesh node (a vertex/point) to the tree.
  void InsertMeshNode(const double x, const double y, const int index);

  /// Insert a mesh element with given bounding box and index to the tree.
  void InsertMeshElement(const double bb[4], const int index);

  /// Get all elements linked to a block corresponding to the given point.
  std::vector<int> GetElementsInBlock(const double x, const double y) const;

 private:
  // Centre of this tree node.
  double m_x0, m_y0;
  // Half-width in x and y of this tree node.
  double m_hx, m_hy;
  // Bounding box of the tree node.
  double m_xmin, m_ymin, m_xmax, m_ymax;

  // Pointers to child quadrants.
  QuadTree* children[4];  

  // Children follow a predictable pattern to make accesses simple.
  // Here, - means less than 'origin' in that dimension, + means greater than.
  // child:	0 1 2 3 
  // x:     - - + +
  // y:     - + - +

  std::vector<std::tuple<float, float, int> > nodes;
  std::vector<int> elements;

  static const size_t BlockCapacity = 10;

  // Check if the given box overlaps with this tree node.
  bool DoesBoxOverlap(const double bb[4]) const;

  int GetQuadrant(const double x, const double y) const;

  // Check if this tree node is a leaf or intermediate node.
  bool IsLeafNode() const;

  // Get a block containing the input point
  const QuadTree* GetBlockFromPoint(const double x, const double y) const;

  // A helper function used by the function above.
  // Called recursively on the child nodes.
  const QuadTree* GetBlockFromPointHelper(const double x, const double y) const;
};
}

#endif
