#ifndef TETRAHEDRAL_TREE_H
#define TETRAHEDRAL_TREE_H

#include <cstddef>
#include <vector>

namespace Garfield {

// TODO: replace this class with ROOT's TVector3 class

struct Vec3 {
  float x = 0., y = 0., z = 0.;

  Vec3() {}
  Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}

  Vec3 operator+(const Vec3& r) const {
    return Vec3(x + r.x, y + r.y, z + r.z);
  }

  Vec3 operator-(const Vec3& r) const {
    return Vec3(x - r.x, y - r.y, z - r.z);
  }

  Vec3& operator+=(const Vec3& r) {
    x += r.x;
    y += r.y;
    z += r.z;
    return *this;
  }

  Vec3& operator-=(const Vec3& r) {
    x -= r.x;
    y -= r.y;
    z -= r.z;
    return *this;
  }

  Vec3 operator*(float r) const { return Vec3(x * r, y * r, z * r); }

  Vec3 operator/(float r) const { return Vec3(x / r, y / r, z / r); }
};

/**

\brief Helper class for searches in field maps.

This class stores the mesh nodes and elements in an Octree data
structure to optimize the element search operations

Author: Ali Sheharyar

Organization: Texas A&M University at Qatar

*/
class TetrahedralTree {
 public:
  /// Constructor
  TetrahedralTree(const Vec3& origin, const Vec3& halfDimension);

  /// Destructor
  ~TetrahedralTree();

  /// Insert a mesh node (a vertex/point) to the tree.
  void InsertMeshNode(Vec3 point, const int index);

  /// Insert a mesh element with given bounding box and index to the tree.
  void InsertMeshElement(const double bb[6], const int index);

  /// Get all elements linked to a block corresponding to the given point.
  std::vector<int> GetElementsInBlock(const Vec3& point) const;

 private:
  // Physical centre of this tree node.
  Vec3 m_origin;
  // Half the width/height/depth of this tree node. 
  Vec3 m_halfDimension;
  // Storing min and max points for convenience
  Vec3 m_min, m_max;  

  // The tree has up to eight children and can additionally store
  // a list of mesh nodes and mesh elements.
  // Pointers to child octants.
  TetrahedralTree* children[8];  

  // Children follow a predictable pattern to make accesses simple.
  // Here, - means less than 'origin' in that dimension, + means greater than.
  // child:	0 1 2 3 4 5 6 7
  // x:     - - - - + + + +
  // y:     - - + + - - + +
  // z:     - + - + - + - +

  std::vector<std::pair<Vec3, int> > nodes;
  std::vector<int> elements;

  static const size_t BlockCapacity = 10;

  // Check if the given box overlaps with this tree node.
  bool DoesBoxOverlap(const double bb[6]) const;

  int GetOctantContainingPoint(const Vec3& point) const;

  // Check if this tree node is a leaf or intermediate node.
  bool IsLeafNode() const;

  // Get a block containing the input point
  const TetrahedralTree* GetBlockFromPoint(const Vec3& point) const;

  // A helper function used by the function above.
  // Called recursively on the child nodes.
  const TetrahedralTree* GetBlockFromPointHelper(const Vec3& point) const;
};
}

#endif
