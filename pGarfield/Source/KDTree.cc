//
// (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
//
// Licensed under the Academic Free License version 1.1 found in file LICENSE
// with additional provisions in that same file.

#include <algorithm> 
#include <limits>
#include <iostream>

#include "Garfield/KDTree.hh"

namespace {

double squared(const double x) { return x * x; }

double dis_from_bnd(const double x, const double amin, const double amax) {
  if (x > amax) {
    return(x-amax); 
  } else if (x < amin)
    return (amin-x);
  else
    return 0.0;
}

}

namespace Garfield {

inline bool operator<(const KDTreeResult& e1, const KDTreeResult& e2) {
  return (e1.dis < e2.dis);
}

// Constructor
KDTree::KDTree(KDTreeArray& data_in)
  : m_data(data_in) {

  const size_t n = data_in.size(); 
  if (!data_in.empty()) {
    m_dim = data_in[0].size();
  } 

  m_ind.resize(n);
  for (size_t i = 0; i < n; i++) m_ind[i] = i; 
  // Build the tree.
  m_root = build_tree_for_range(0, n - 1, 0); 
}

// Destructor
KDTree::~KDTree() {
  delete m_root;
}

KDTreeNode* KDTree::build_tree_for_range(int l, int u, KDTreeNode* parent) {

  if (u < l) return nullptr;
  KDTreeNode* node = new KDTreeNode(m_dim);
  if ((u - l) <= bucketsize) {
    // Create a terminal node. Always compute true bounding box. 
    for (size_t i = 0; i < m_dim; i++) {
      node->box[i] = spread_in_coordinate(i, l, u);
    }
    node->cut_dim = 0; 
    node->cut_val = 0.0;
    node->m_l = l;
    node->m_u = u;
    node->left = node->right = nullptr;
  } else {
    // Compute an APPROXIMATE bounding box for this node.
    // If parent == nullptr, then this is the root node, and 
    // we compute for all dimensions.
    // Otherwise, we copy the bounding box from the parent for
    // all coordinates except for the parent's cut dimension. 
    // That, we recompute ourself.
    int c = -1;
    double maxspread = 0.0;
    for (size_t i = 0; i < m_dim; i++) {
      if (!parent || (parent->cut_dim == i)) {
        node->box[i] = spread_in_coordinate(i, l, u);
      } else {
        node->box[i] = parent->box[i];
      }
      const double spread = node->box[i][1] - node->box[i][0];
      if (spread > maxspread) {
        maxspread = spread;
        c = i; 
      }
    }

    // Now, c is the identity of which coordinate has the greatest spread.
    double sum = 0.0;
    for (int k = l; k <= u; k++) {
      sum += m_data[m_ind[k]][c];
    }
    const double average = sum / static_cast<double>(u - l + 1);
    int m = select_on_coordinate_value(c, average, l, u);

    // Move the indices around to cut on dim 'c'.
    node->cut_dim = c;
    node->m_l = l;
    node->m_u = u;
    node->left = build_tree_for_range(l, m, node);
    node->right = build_tree_for_range(m + 1, u, node);

    if (!node->right) {
      node->box = node->left->box;
      node->cut_val = node->left->box[c][1];
      node->cut_val_left = node->cut_val_right = node->cut_val;
    } else if (!node->left) {
      node->box = node->right->box;
      node->cut_val = node->right->box[c][1];
      node->cut_val_left = node->cut_val_right = node->cut_val;
    } else {
      node->cut_val_right = node->right->box[c][0];
      node->cut_val_left  = node->left->box[c][1];
      node->cut_val = 0.5 * (node->cut_val_left + node->cut_val_right); 
      
      // Now recompute true bounding box as union of subtree boxes.
      // This is now faster having built the tree, being logarithmic in
      // N, not linear as would be from naive method.
      for (size_t i = 0; i < m_dim; i++) {
        node->box[i][1] = std::max(node->left->box[i][1],
                                   node->right->box[i][1]);
        node->box[i][0] = std::min(node->left->box[i][0],
                                   node->right->box[i][0]);
      }
    }
  }
  return node;
}

std::array<double, 2> KDTree::spread_in_coordinate(const int c, const int l,
                                                   const int u) const {
  // Return the minimum and maximum of the indexed data between l and u.

  double smin = m_data[m_ind[l]][c];
  double smax = smin;

  // Process two at a time.
  int i; 
  for (i = l + 2; i <= u; i += 2) {
    double lmin = m_data[m_ind[i - 1]][c];
    double lmax = m_data[m_ind[i]][c];
    if (lmin > lmax) std::swap(lmin, lmax); 
    if (smin > lmin) smin = lmin;
    if (smax < lmax) smax = lmax;
  }
  // Is there one more element? 
  if (i == u + 1) {
    double last = m_data[m_ind[u]][c];
    if (smin > last) smin = last;
    if (smax < last) smax = last;
  }
  return {smin, smax};
}

int KDTree::select_on_coordinate_value(int c, double alpha, int l, int u) {
  //  Move indices in ind[l..u] so that the elements in [l .. return]
  //  are <= alpha, and hence are less than the [return + 1 .. u]
  //  elements, viewed across dimension 'c'.
  int lb = l, ub = u;
  while (lb < ub) {
    if (m_data[m_ind[lb]][c] <= alpha) {
      lb++; // good where it is.
    } else {
      std::swap(m_ind[lb], m_ind[ub]);
      ub--;
    }
  }
  // Here ub = lb.
  return m_data[m_ind[lb]][c] <= alpha ? lb : lb - 1;
}

void KDTree::n_nearest(const std::vector<double>& qv, 
                       const unsigned int nn, 
                       std::vector<KDTreeResult>& result) const {
  // Search for n nearest to a given query vector 'qv'.
  std::priority_queue<KDTreeResult> res; 
  double r2 = std::numeric_limits<double>::max();
  m_root->search_n(-1, 0, nn, r2, qv, *this, res);
  result.clear();
  while (!res.empty()) {
    result.push_back(res.top());
    res.pop();
  }
  if (sort_results) sort(result.begin(), result.end());
}

void KDTree::n_nearest_around_point(const unsigned int idx, 
                                    const unsigned int ndecorrel, 
                                    const unsigned int nn,
                                    std::vector<KDTreeResult>& result) const {

  std::priority_queue<KDTreeResult> res;
  double r2 = std::numeric_limits<double>::max();
  m_root->search_n(idx, ndecorrel, nn, r2, m_data[idx], *this, res);
  result.clear(); 
  while (!res.empty()) {
    result.push_back(res.top());
    res.pop();
  }
  if (sort_results) sort(result.begin(), result.end());
}

void KDTree::r_nearest(const std::vector<double>& qv, const double r2, 
                       std::vector<KDTreeResult>& result) const {
  // Search for all within a ball of a certain radius.
  result.clear(); 
  m_root->search_r(-1, 0, r2, qv, *this, result);
  if (sort_results) sort(result.begin(), result.end());
} 

void KDTree::r_nearest_around_point(const unsigned int idx, 
                                    const unsigned int ndecorrel, 
                                    const double r2,
                                    std::vector<KDTreeResult>& result) const {

  result.clear(); 
  m_root->search_r(idx, ndecorrel, r2, m_data[idx], *this, result);
  if (sort_results) sort(result.begin(), result.end());
}

// Constructor
KDTreeNode::KDTreeNode(int dim) : box(dim) {} 

// Destructor
KDTreeNode::~KDTreeNode() {
  if (left) delete left; 
  if (right) delete right; 
}

void KDTreeNode::search_n(const int idx0, const int nd,
                          const unsigned int nn, double& r2, 
                          const std::vector<double>& qv, const KDTree& tree,
                          std::priority_queue<KDTreeResult>& res) const {

  if (!left && !right) {
    // We are on a terminal node.
    process_terminal_node_n(idx0, nd, nn, r2, qv, tree, res);
    return;
  }
  KDTreeNode *ncloser = nullptr;
  KDTreeNode *nfarther = nullptr;

  double extra;
  double qval = qv[cut_dim]; 
  // value of the wall boundary on the cut dimension. 
  if (qval < cut_val) {
    ncloser = left;
    nfarther = right;
    extra = cut_val_right - qval;
  } else {
    ncloser = right;
    nfarther = left;
    extra = qval - cut_val_left; 
  }

  if (ncloser) ncloser->search_n(idx0, nd, nn, r2, qv, tree, res);
  if ((nfarther) && (squared(extra) < r2)) {
    // first cut
    if (nfarther->box_in_search_range(r2, qv)) {
      nfarther->search_n(idx0, nd, nn, r2, qv, tree, res); 
    }      
  }
}

void KDTreeNode::search_r(const int idx0, const int nd, const double r2,
                          const std::vector<double>& qv, const KDTree& tree,
                          std::vector<KDTreeResult>& res) const {

  if (!left && !right) {
    // We are on a terminal node.
    process_terminal_node_r(idx0, nd, r2, qv, tree, res);
    return;
  }
  KDTreeNode *ncloser = nullptr;
  KDTreeNode *nfarther = nullptr;

  double extra;
  double qval = qv[cut_dim];
  // value of the wall boundary on the cut dimension. 
  if (qval < cut_val) {
    ncloser = left;
    nfarther = right;
    extra = cut_val_right - qval;
  } else {
    ncloser = right;
    nfarther = left;
    extra = qval - cut_val_left; 
  }

  if (ncloser) ncloser->search_r(idx0, nd, r2, qv, tree, res);
  if ((nfarther) && (squared(extra) < r2)) {
    // first cut
    if (nfarther->box_in_search_range(r2, qv)) {
      nfarther->search_r(idx0, nd, r2, qv, tree, res); 
    }      
  }
}

inline bool KDTreeNode::box_in_search_range(const double r2,
                                            const std::vector<double>& qv) const {

  // Does the bounding box have any point which is within 'r2' to 'qv'??
 
  const size_t dim = qv.size();
  double dis2 = 0.0; 
  for (size_t i = 0; i < dim; i++) {
    dis2 += squared(dis_from_bnd(qv[i], box[i][0], box[i][1]));
    if (dis2 > r2) return false;
  }
  return true;
}

void KDTreeNode::process_terminal_node_n(const int idx0, const int nd,
    const unsigned int nn, double& r2, const std::vector<double>& qv, 
    const KDTree& tree, std::priority_queue<KDTreeResult>& res) const {

  const size_t dim = tree.m_dim;
  const auto& data = tree.m_data;

  for (int i = m_l; i <= m_u; i++) {
    const int idx = tree.m_ind[i];
    bool early_exit = false;
    double dis = 0.0;
    for (size_t k = 0; k < dim; k++) {
      dis += squared(data[idx][k] - qv[k]);
      if (dis > r2) {
        early_exit = true; 
        break;
      }
    }
    if (early_exit) continue; // next iteration of mainloop

    // Skip points within the decorrelation interval. 
    if (idx0 >= 0 && (abs(idx - idx0) < nd)) continue;

    // Add the point to the list.
    if (res.size() < nn) {
      // The list so far is undersized. 
      KDTreeResult e;
      e.idx = idx;
      e.dis = dis;
      res.push(e); 
      // Set the ball radius to the largest on the list (maximum priority).
      if (res.size() == nn) r2 = res.top().dis;
    } else {
      // if we get here then the current node, has a squared 
      // distance smaller
      // than the last on the list, and belongs on the list.
      KDTreeResult e;
      e.idx = idx;
      e.dis = dis;
      res.pop();
      res.push(e); 
      r2 = res.top().dis;
    }
  } // main loop
}

void KDTreeNode::process_terminal_node_r(const int idx0, const int nd,
    const double r2, const std::vector<double>& qv, const KDTree& tree,
    std::vector<KDTreeResult>& res) const {

  const size_t dim = tree.m_dim;
  const auto& data = tree.m_data;

  for (int i = m_l; i <= m_u; i++) {
    const int idx = tree.m_ind[i]; 
    bool early_exit = false;
    double dis = 0.0;
    for (size_t k = 0; k < dim; k++) {
      dis += squared(data[idx][k] - qv[k]);
      if (dis > r2) {
        early_exit = true; 
        break;
      }
    }
    if (early_exit) continue; // next iteration of mainloop

    // Skip points within the decorrelation interval.   
    if (idx0 >= 0 && (abs(idx - idx0) < nd)) continue;

    KDTreeResult e;
    e.idx = idx;
    e.dis = dis;
    res.push_back(std::move(e));
  }
}

}
