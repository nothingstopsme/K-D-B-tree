# K-D-B-tree
This is a c++ implementation of K-D-B-tree introduced in this [paper](https://dl.acm.org/doi/pdf/10.1145/582318.582321)

## Assumptions
1. All dimensions are represented in `double`
2. The distance (squared) is defined as the Euclidean distance in that multi-dimensional space

## Features
1. Insertion/Deletion with self-balancing.
2. Searching for $k$ nearest neighbours, with $k \ge 1, k \in N$.
3. Searching for near neighbours within a given radius (squared).
4. Searching for all the points inside a given region (defined as a hyperrectangle in the search space), including those right on the region boundaries.

## Notes
1. While one of the advantages of the likes of B-tree is the potential of better memory management via page utilisation, this implementation is not optimised for that.
2. The dimension-choosing strategy for splitting/joining adopted in this implementation is as follows:
   * Splitting
     - When a node collection (corresponding to "point page" in the paper) is split, the splitting dimension/line is picked from the coordinates of those nodes; and among all possible choices of coordinates, the one which can produce the most balanced node distribution, i.e. the number of nodes in the resulting left region and the one in the resulting right region after splitting are the closest, is selected.
     - When a region collection (corresponding to "region page" in the paper) is split, the splitting dimension/line is picked from the boundary lines of those regions; and among all possible choices of boundary lines, the one which borders the most regions, with the tie break being whether any of other regions would be cut through by that splitting (preferring not), is selected.
   * Joining
     
     If a node/region collection is going to be joined with others, the joining dimension/line is picked from the boundary lines of the immediate region where that collection belongs; and among all possible choices of boundary lines, the one which borders the fewest regions is selected, with the corresponding neighbouring regions being those to be joined with it.
     
3. I believe that the minimum size of collections guaranteed for all dimensionalities under any splitting/joining strategies (including the one used in this implementation) would be 1, as there would always be some extreme cases in a high-dimensional space where (nearly) balanced separations are not possible. For example, when it comes to a collection of 5 points $[(1, 1, 1, 1),\ (1, 1, 1, 2),\ (1, 1, 2, 1),\ (1, 2, 1, 1),\ (2, 1, 1, 1)]$, the best result one can achieve, with only one cut along a single dimension, is a 1-4 split.

## Usage
First, inherit `NodeInterface` and provide your node definition
```c++
#include "kdb_tree.h"

class TestNode : public kdb_tree::NodeInterface {
  public:
    TestNode(const double x, const double y, const double z)
      : x_(x), y_(y), z_(z) {
      }

    virtual double operator [](const unsigned int dimension) const override {
      switch(dimension) {
        case 0:
          return x_;
        case 1:
          return y_;
        case 2:
          return z_;
        default:
          throw std::logic_error("A wrong dimension is given");
      }
    }

    // optional
    virtual void Print(std::ostream& os) const override {
      os << "(" << x_ << "," << y_ << "," << z_ << ")";
    }

  private:
    double x_;
    double y_;
    double z_;
};
```
then a tree for this particular node definition can be created/manipulated/queried like this:
```c++
using TestTree = kdb_tree::Tree<3, TestNode>;
TestTree tree;

auto node0 = std::make_shared<TestTree::Node>(0, 0, 0);
auto node1 = std::make_shared<TestTree::Node>(1, 1, 1);
tree.Insert(node0);
tree.Insert(node1);
std::cout << "The number of nodes after insertion = " << tree.NodeCount() << std::endl;

// TestTree::NearList is a list containing all nodes found,
// in ascending order with respect to the distance (squared) to each node
TestTree::NearList nearests = tree.FindNearests(*node0, 2);
std::cout << "The nearest node: " << *nearests.front().first
          << ", dist squared = " << nearests.front().second << std::endl;

TestTree::NearList nears = tree.FindNears(*node1, 3);
std::cout << "The farthest node with a distance squared <= 3 : "
          << *nears.back().first << ", dist squared = " << nears.back().second << std::endl;

tree.Delete(node1);
tree.Delete(node0);
std::cout << "The number of nodes after deletion = " << tree.NodeCount() << std::endl;
```
Please check functions in the [unit_test](kdb_tree.cpp#L15) namespace for more examples
