#ifndef G_KDTREE2_H
#define G_KDTREE2_H

// (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
//
// Licensed under the Academic Free License version 1.1 found in file LICENSE
// with additional provisions in that same file.


// Implement a kd tree for fast searching of points in a fixed data base
// in k-dimensional Euclidean space.

#include <vector>
#include <array>
#include <queue>
#include <algorithm>

namespace Garfield {

typedef std::vector<std::vector<double> > KDTreeArray;

class KDTreeNode; 

/// Search result

struct KDTreeResult {
  double dis; //< square Euclidean distance
  size_t idx; //< index
}; 

/// Main k-d tree class.
/// Fast search of points in k-dimensional Euclidean space.

class KDTree {
public: 
  // Reference to the underlying data to be included in the tree.
  const KDTreeArray& m_data;   

  size_t m_dim;
  bool sort_results = false;

public:
  KDTree() = delete;
  /// Constructor.
  KDTree(KDTreeArray& data_in);
  /// Destructor.
  ~KDTree();

  /** Search for nn nearest neighbours around a point.
    * \param qv input point 
    * \param nn number of nearest neighbours
    * \param result indices and distances of the nearest neighbours
    */
  void n_nearest(const std::vector<double>& qv, const unsigned int nn, 
                 std::vector<KDTreeResult>& result) const;

  /** Search for nn nearest neighbours around a node of the input data, 
    * excluding neighbors within a decorrelation interval.
    * \param idx index of the input point
    * \param ndecorrel decorrelation interval
    * \param nn number of nearest neighbours
    * \param result indices and distances of the nearest neighbours
    */
  void n_nearest_around_point(const unsigned int idx, 
                              const unsigned int ndecorrel, 
                              const unsigned int nn,
                              std::vector<KDTreeResult>& result) const;
  
  /** Search for all neighbors in a ball of size r2 
    * \param qv input point
    * \param r2 ball size (square Euclidean distance)
    * \param result indices and distances of the nearest neighbours
    */ 
  void r_nearest(const std::vector<double>& qv, const double r2,
                 std::vector<KDTreeResult>& result) const;

  /// Like r_nearest, but around an existing point, 
  /// with decorrelation interval. 
  void r_nearest_around_point(const unsigned int idx, 
                              const unsigned int ndecorrel, const double r2,
                              std::vector<KDTreeResult>& result) const;

  friend class KDTreeNode;
private:
  KDTreeNode* m_root = nullptr;

  // Index for the tree leaves. Data in a leaf with bounds [l,u] are
  // in 'data[ind[l],*] to data[ind[u],*]
  std::vector<size_t> m_ind; 

  static constexpr int bucketsize = 12; // global constant. 

private:
  KDTreeNode* build_tree_for_range(int l, int u, KDTreeNode* parent);
  int select_on_coordinate_value(int c, double alpha, int l, int u); 
  std::array<double, 2> spread_in_coordinate(const int c, const int l, const int u) const;
};

/// A node in the k-d tree.

class KDTreeNode {
public:
  /// Constructor
  KDTreeNode(int dim);
  /// Destructor
  ~KDTreeNode();

private:
  friend class KDTree;

  // Dimension to cut.
  size_t cut_dim = 0; 
  // Cut value.
  double cut_val = 0.;
  double cut_val_left = 0.;
  double cut_val_right = 0.;
  // Extents in index array for searching
  int m_l = 0;
  int m_u = 0;
  // [min,max] of the box enclosing all points
  std::vector<std::array<double, 2> > box; 

  // Pointers to left and right nodes.
  KDTreeNode *left = nullptr;
  KDTreeNode *right = nullptr;  

  // Recursive innermost core routines for searching.
  void search_n(const int idx0, const int nd,
                const unsigned int nn, double& r2,
                const std::vector<double>& qv, const KDTree& tree, 
                std::priority_queue<KDTreeResult>& res) const; 
  void search_r(const int idx0, const int nd, const double r2,
                const std::vector<double>& qv, const KDTree& tree, 
                std::vector<KDTreeResult>& res) const;
  
  // Return true if the bounding box for this node is within the
  // search range around a point.
  bool box_in_search_range(const double r2, 
                           const std::vector<double>& qv) const;

  // For processing final buckets. 
  void process_terminal_node_n(const int idx0, const int nd,
                               const unsigned int nn, double& r2, 
                               const std::vector<double>& qv, 
                               const KDTree& tree,
                               std::priority_queue<KDTreeResult>& res) const;
  void process_terminal_node_r(const int idx0, const int nd,
                               const double r2, 
                               const std::vector<double>& qv, 
                               const KDTree& tree,
                               std::vector<KDTreeResult>& res) const;

};

}

#endif
