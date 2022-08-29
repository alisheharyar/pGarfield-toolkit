#include "Garfield/TetrahedralTree.hh"
#include <iostream>

namespace Garfield {

/**
TetrahedralTree.cc
This class stores the mesh nodes and elements in an Octree data
structure to optimize the element search operations

Author: Ali Sheharyar
Organization: Texas A&M University at Qatar
*/
TetrahedralTree::TetrahedralTree(const Vec3& origin, const Vec3& halfDimension)
    : m_origin(origin), m_halfDimension(halfDimension) {
  m_min.x = origin.x - halfDimension.x;
  m_min.y = origin.y - halfDimension.y;
  m_min.z = origin.z - halfDimension.z;
  m_max.x = origin.x + halfDimension.x;
  m_max.y = origin.y + halfDimension.y;
  m_max.z = origin.z + halfDimension.z;

  // Initially, there are no children
  for (int i = 0; i < 8; ++i) children[i] = nullptr;
}

TetrahedralTree::~TetrahedralTree() {
  // Recursively destroy octants
  for (int i = 0; i < 8; ++i) delete children[i];
}

// Check if a box overlaps with this node
bool TetrahedralTree::DoesBoxOverlap(const double bb[6]) const {
  if (m_max.x < bb[0] || m_max.y < bb[1] || m_max.z < bb[2]) return false;
  if (m_min.x > bb[3] || m_min.y > bb[4] || m_min.z > bb[5]) return false;
  return true;
}

// Determine which octant of the tree would contain 'point'
int TetrahedralTree::GetOctantContainingPoint(const Vec3& point) const {
  int oct = 0;
  if (point.x >= m_origin.x) oct |= 4;
  if (point.y >= m_origin.y) oct |= 2;
  if (point.z >= m_origin.z) oct |= 1;
  return oct;
}

bool TetrahedralTree::IsLeafNode() const {
  // We are a leaf if we have no children. Since we either have none, or
  // all eight, it is sufficient to just check the first.
  return children[0] == nullptr;
}

void TetrahedralTree::InsertMeshNode(Vec3 point, const int index) {
  // Check if it is a leaf node.
  if (!IsLeafNode()) {
    // We are at an interior node. Insert recursively into the
    // appropriate child octant.
    int octant = GetOctantContainingPoint(point);
    children[octant]->InsertMeshNode(point, index);
    return;
  }
  
  // Add the new point if the block is not full.
  if (nodes.size() < BlockCapacity) {
    nodes.push_back(std::make_pair(point, index));
    return;
  } 
  // Block is full, so we need to partition it.
  // Split the current node and create new empty trees for each child octant.
  for (int i = 0; i < 8; ++i) {
    // Compute new bounding box for this child
    Vec3 newOrigin = m_origin;
    newOrigin.x += m_halfDimension.x * (i & 4 ? .5f : -.5f);
    newOrigin.y += m_halfDimension.y * (i & 2 ? .5f : -.5f);
    newOrigin.z += m_halfDimension.z * (i & 1 ? .5f : -.5f);
    children[i] = new TetrahedralTree(newOrigin, m_halfDimension * .5f);
  }

  // Move the mesh nodes from the partitioned node (now marked as interior) to
  // its children.
  while (!nodes.empty()) {
    auto node = nodes.back();
    nodes.pop_back();
    const int oct = GetOctantContainingPoint(node.first);
    children[oct]->InsertMeshNode(node.first, node.second);
  }
  // Insert the new point in the appropriate octant.
  children[GetOctantContainingPoint(point)]->InsertMeshNode(point, index);
}

void TetrahedralTree::InsertMeshElement(const double bb[6], const int index) {
  if (IsLeafNode()) {
    // Add the element to the list of this octant.
    elements.push_back(index);
    return;
  } 
  // Check which children overlap with the element's bounding box.
  for (int i = 0; i < 8; i++) {
    if (!children[i]->DoesBoxOverlap(bb)) continue;
    children[i]->InsertMeshElement(bb, index);
  }
}

// It returns the list of tetrahedrons that intersects in a bounding box (Octree
// block) that contains the
// point passed as input.
std::vector<int> TetrahedralTree::GetElementsInBlock(const Vec3& point) const {
  const TetrahedralTree* octreeNode = GetBlockFromPoint(point);

  if (octreeNode) {
    return octreeNode->elements;
  }

  return std::vector<int>();
}

// check if the point is inside the domain.
// This function is only executed at root to ensure that input point is inside
// the mesh's bounding box
// If we don't check this, the case when root is leaf node itself will return
// wrong block
const TetrahedralTree* TetrahedralTree::GetBlockFromPoint(
    const Vec3& point) const {
  if (!(m_min.x <= point.x && point.x <= m_max.x &&
        m_min.y <= point.y && point.y <= m_max.y &&
        m_min.z <= point.z && point.z <= m_max.z))
    return nullptr;

  return GetBlockFromPointHelper(point);
}

const TetrahedralTree* TetrahedralTree::GetBlockFromPointHelper(
    const Vec3& point) const {
  // If we're at a leaf node, it means, the point is inside this block
  if (IsLeafNode()) return this;
  // We are at the interior node, so check which child octant contains the
  // point
  int octant = GetOctantContainingPoint(point);
  return children[octant]->GetBlockFromPointHelper(point);
}
}
