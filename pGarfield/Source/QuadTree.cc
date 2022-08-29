#include "Garfield/QuadTree.hh"
#include <iostream>

namespace Garfield {

QuadTree::QuadTree(const double x0, const double y0, 
                   const double hx, const double hy) 
    : m_x0(x0), m_y0(y0), m_hx(hx), m_hy(hy) {
  m_xmin = x0 - hx;
  m_ymin = y0 - hy;
  m_xmax = x0 + hx;
  m_ymax = y0 + hy;

  // Initially, there are no children.
  for (int i = 0; i < 4; ++i) children[i] = nullptr;
}

QuadTree::~QuadTree() {
  for (int i = 0; i < 4; ++i) delete children[i];
}

bool QuadTree::DoesBoxOverlap(const double bb[4]) const {
  if (m_xmax < bb[0] || m_ymax < bb[1]) return false;
  if (m_xmin > bb[2] || m_ymin > bb[3]) return false;
  return true;
}

int QuadTree::GetQuadrant(const double x, const double y) const {
  int quad = 0;
  if (x >= m_x0) quad |= 2;
  if (y >= m_y0) quad |= 1;
  return quad;
}

bool QuadTree::IsLeafNode() const {
  // We are a leaf if we have no children. Since we either have none or
  // all, it is sufficient to just check the first.
  return children[0] == nullptr;
}

void QuadTree::InsertMeshNode(const double x, const double y, const int index) {
  // Check if it is a leaf node.
  if (!IsLeafNode()) {
    // We are at an interior node. 
    // Insert recursively into appropriate child quadrant.
    int quad = GetQuadrant(x, y);
    children[quad]->InsertMeshNode(x, y, index);
    return;
  }
  
  // Add the new point if the block is not full.
  if (nodes.size() < BlockCapacity) {
    nodes.push_back(std::make_tuple(x, y, index));
    return;
  } 
  // Block is full, so we need to partition it.
  // Split the current node and create new empty trees for each child.
  for (int i = 0; i < 4; ++i) {
    // Compute new bounding box for this child.
    const double xi = m_x0 + m_hx * (i & 2 ? 0.5 : -0.5);
    const double yi = m_y0 + m_hy * (i & 1 ? 0.5 : -0.5);
    children[i] = new QuadTree(xi, yi, 0.5 * m_hx, 0.5 * m_hy);
  }

  // Move the mesh nodes from the partitioned node (now marked as interior) to
  // its children.
  while (!nodes.empty()) {
    auto node = nodes.back();
    nodes.pop_back();
    const double xn = std::get<0>(node);
    const double yn = std::get<1>(node);
    const int quad = GetQuadrant(xn, yn);
    children[quad]->InsertMeshNode(xn, yn, std::get<2>(node));
  }
  // Insert the new point in the appropriate octant.
  children[GetQuadrant(x, y)]->InsertMeshNode(x, y, index);
}

void QuadTree::InsertMeshElement(const double bb[4], const int index) {
  if (IsLeafNode()) {
    // Add the element to the list of this quadrant.
    elements.push_back(index);
    return;
  } 
  // Check which children overlap with the element's bounding box.
  for (int i = 0; i < 4; ++i) {
    if (!children[i]->DoesBoxOverlap(bb)) continue;
    children[i]->InsertMeshElement(bb, index);
  }
}

std::vector<int> QuadTree::GetElementsInBlock(const double x, 
                                              const double y) const {
  const auto node = GetBlockFromPoint(x, y);
  if (node) return node->elements;
  return std::vector<int>();
}

const QuadTree* QuadTree::GetBlockFromPoint(
    const double x, const double y) const {
  if (x < m_xmin || x > m_xmax || y < m_ymin || y > m_ymax) {
    return nullptr;
  }
  return GetBlockFromPointHelper(x, y);
}

const QuadTree* QuadTree::GetBlockFromPointHelper(
    const double x, const double y) const {
  // If we're at a leaf node, it means, the point is inside this block.
  if (IsLeafNode()) return this;
  // We are at the interior node, so check which child contains the point.
  int quad = GetQuadrant(x, y);
  return children[quad]->GetBlockFromPointHelper(x, y);
}
}
