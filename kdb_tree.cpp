#include "kdb_tree.h"
#include <vector>
#include <algorithm>
#include <random>
#include <sstream>

#define OK_COLOUR "\033[32m"
#define NOTICE_COLOUR "\033[36m"
#define ERROR_COLOUR "\033[31m"
#define WARNING_COLOUR "\033[33m"
#define COLOUR_END "\033[0m"

namespace kdb_tree {

namespace unit_test {
class TestNode : public NodeInterface {
  public:
    TestNode(const double x, const double y, const double z);
    virtual double operator [](const unsigned int dimension) const override;
    virtual void Print(std::ostream& os) const override;

  private:
    double x_;
    double y_;
    double z_;
};

TestNode::TestNode(const double x, const double y, const double z)
: x_(x), y_(y), z_(z)
{}


double TestNode::operator[](const unsigned int dimension) const {
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

void TestNode::Print(std::ostream& os) const {
  os << "(" << x_ << "," << y_ << "," << z_ << ")";
}

using TestTree = Tree<3, TestNode>;

void Test_InsertThenDelete(const unsigned int num_of_nodes, const unsigned int seed) {
  std::minstd_rand engine(seed);
  std::uniform_int_distribution<int> distribution(-10, 10);
  std::cout << NOTICE_COLOUR << __func__ << COLOUR_END << std::endl;
  std::cout << "Inserting " << num_of_nodes << " nodes with integer coordinates, random seed = " << seed << std::endl;
  
  TestTree tree;
  TestTree::NodePtrList all_nodes;
  try {
    unsigned int node_count = 0;
    for (unsigned int n = 0; n < num_of_nodes; ++n) {
      auto node = std::make_shared<TestTree::Node>(distribution(engine), 
          distribution(engine),
          distribution(engine));
      if (tree.Insert(node)) {
        node_count += 1;

        std::cout << NOTICE_COLOUR << "Inserting: " << COLOUR_END << *node << std::endl
                  << tree << std::endl;
        all_nodes.emplace_back(std::move(node));
      }
      else
        std::cout << WARNING_COLOUR << "A duplicated node is detected: " << *node << "; skipping it." << COLOUR_END << std::endl;
    }
    std::cout << "Number of valid nodes = " << node_count << std::endl;

    if (!all_nodes.empty()) {
      std::cout << std::endl << "Deleting inserted nodes, in reverse order" << std::endl;
      for (std::list<TestTree::NodePtr>::iterator iter = --all_nodes.end();; --iter) {
      
        if (tree.Delete(*iter))
          std::cout << NOTICE_COLOUR << "Deleting: " << COLOUR_END << **iter << std::endl
                  << tree << std::endl;
        else {
          std::cout << WARNING_COLOUR << "Failed to delete this node: " << **iter << COLOUR_END << std::endl;
        }

        if (iter == all_nodes.begin())
          break;
      }
    }
  }
  catch (const std::exception& e) {
    std::cout << ERROR_COLOUR << e.what() << COLOUR_END << std::endl;
  }
  

}

template <typename COORD_DICE>
void PrepareRandomTestTree(const unsigned int expected_num_of_nodes, std::minstd_rand& random_engine, COORD_DICE& coord_dice, TestTree& tree, TestTree::NodePtrList& inserted_nodes) {

  inserted_nodes.clear();

  std::vector<TestTree::NodePtrList::iterator> shuffle_vector;
  shuffle_vector.reserve(expected_num_of_nodes*2);

  for (unsigned int n = 0; n < expected_num_of_nodes*2; ++n) {
    bool result = false;
    double coords[] = {coord_dice(random_engine), coord_dice(random_engine), coord_dice(random_engine)};
    auto node = std::make_shared<TestTree::Node>(coords[0], coords[1], coords[2]); 
    if (n % 2 == 0) {
      result = tree.Insert(node);
    }
    else {
      result = tree.Insert(std::move(*node));
      if (result)
        node = tree.FindNearests({coords[0], coords[1], coords[2]}).front().first;
    }
    
    if (result) {
      inserted_nodes.emplace_back(std::move(node));
      shuffle_vector.emplace_back(--inserted_nodes.end());
    }
    else
      std::cout << WARNING_COLOUR << "A duplicated node is detected: " << *node << "; skipping it." << COLOUR_END << std::endl;
  }

  
  std::uniform_int_distribution<unsigned int> shuffle_dice (0, shuffle_vector.size()-1);
  for (unsigned int index = 0; index < shuffle_vector.size(); ++index) {
    unsigned int random_index = shuffle_dice(random_engine);
    if (random_index != index) {
      auto temp = shuffle_vector[index];
      shuffle_vector[index] = shuffle_vector[random_index];
      shuffle_vector[random_index] = temp;
    }
  }

  if (inserted_nodes.size() <= expected_num_of_nodes) {
    
    std::cout << WARNING_COLOUR << "The number of inserted nodes does not exceed the expected one (" << expected_num_of_nodes << ")"
              << "; no deletion is performed" << COLOUR_END << std::endl;
    return;
  }

  unsigned int deletion_index = 0;
  while (inserted_nodes.size() > expected_num_of_nodes) {
    TestTree::NodePtrList::iterator deletion_iter = shuffle_vector[deletion_index];
    if (deletion_index % 2 == 0) {
      tree.Delete(*deletion_iter);
    }
    else {
      tree.Delete({(**deletion_iter)[0], (**deletion_iter)[1], (**deletion_iter)[2]});
    }
    
    inserted_nodes.erase(deletion_iter);
    deletion_index++;
  }  

}

void Test_FindNearestsImpl(const TestTree& tree, const TestTree::NodePtrList& all_nodes, const unsigned int k, const unsigned int num_of_tests, std::minstd_rand& random_engine, std::uniform_real_distribution<double>& coord_dice) {

  std::cout << "Testing with " << num_of_tests << " random centres" << std::endl;

  unsigned int mismatch_count = 0;
  for (unsigned int n = 0; n < num_of_tests; ++n) {
    auto centre = TestNode(coord_dice(random_engine), 
                            coord_dice(random_engine),
                            coord_dice(random_engine));

    //brute-force search
    std::multimap<double, NodeInterfacePtr> bf_nearests;
    for (const auto& candidate : all_nodes) {
      double dist_squared = 0.0;
      for (unsigned int d = 0; d < tree.K; ++d) {
        double diff = centre[d] - (*candidate)[d];
        dist_squared += diff * diff;
      }

      if (bf_nearests.size() < k || dist_squared < (--(bf_nearests.end()))->first) {
        bf_nearests.emplace(dist_squared, candidate);
        if (bf_nearests.size() > k)
          bf_nearests.erase(--(bf_nearests.end()));
      }
    }

    TestTree::NearList kd_nearests; 
    const bool using_coords(n % 2 == 0);
    if (using_coords)
      kd_nearests = tree.FindNearests({centre[0], centre[1], centre[2]}, k);
    else
      kd_nearests = tree.FindNearests(centre, k);
    bool mis_matched = false;
    decltype(kd_nearests)::iterator kd_iter;
    decltype(bf_nearests)::iterator bf_iter;
    if (kd_nearests.size() != bf_nearests.size())
      mis_matched = true;
    else {
      for (kd_iter = kd_nearests.begin(), bf_iter = bf_nearests.begin(); 
            kd_iter != kd_nearests.end() && bf_iter != bf_nearests.end();
            ++kd_iter, ++bf_iter) {        
        
        if (kd_iter->first.get() != bf_iter->second.get() || kd_iter->second != bf_iter->first
          ) {
          std::cout << WARNING_COLOUR << "kd: " << *kd_iter->first << "(" << kd_iter->second  << "), bf: " << *bf_iter->second << "(" << bf_iter->first << ")" << COLOUR_END << std::endl;
          mis_matched = true;
          break;
        }

      }
    }

    
    if (mis_matched) {
      std::cout << WARNING_COLOUR 
      << "Mismatch detected: centre = " << centre << ", using_coords = " << using_coords << ", brute-force and kd produce different lists of " << k << "-nearest neighbours" << COLOUR_END << std::endl;
      mismatch_count += 1;
    }

  }

  if (mismatch_count)
    std::cout << ERROR_COLOUR << "mismatch_count = " << mismatch_count << COLOUR_END << std::endl;
  else
    std::cout << OK_COLOUR << "Pass" << COLOUR_END << std::endl;


}

void Test_FindNearests(const unsigned int k, const unsigned int seed) {
  std::cout << NOTICE_COLOUR << __func__ << COLOUR_END << ": random seed = " << seed  << std::endl;
  std::minstd_rand random_engine(seed);
  std::uniform_real_distribution<double> coord_dice (-20.0, 20.0);
  constexpr unsigned int NUM_OF_NODES = 100;

  TestTree tree;
  TestTree::NodePtrList all_nodes;


  try {
    PrepareRandomTestTree(NUM_OF_NODES, random_engine, coord_dice, tree, all_nodes);
    std::cout << "Number of valid nodes = " << all_nodes.size() << std::endl;  
  }
  catch (const std::exception& e) {
    std::cout << ERROR_COLOUR << e.what() << COLOUR_END << std::endl;    
    return;
  }

  Test_FindNearestsImpl(tree, all_nodes, k, NUM_OF_NODES, random_engine, coord_dice);


}


void Test_FindNearsImpl(const TestTree& tree, const TestTree::NodePtrList& all_nodes, const unsigned int num_of_tests, std::minstd_rand& random_engine, std::uniform_real_distribution<double>& coord_dice) {

  std::cout << "Testing with " << num_of_tests << " random (centre, radius) pairs" << std::endl;
  std::uniform_real_distribution<double> radius_dice (0.0, 200.0);

  unsigned int mismatch_count = 0;
  for (unsigned int n = 0; n < num_of_tests; ++n) {
    auto centre = TestNode(coord_dice(random_engine), 
                            coord_dice(random_engine),
                            coord_dice(random_engine));
    double radius_squared = radius_dice(random_engine); 
    unsigned int dimensionality = tree.K;

    auto comparator = [dimensionality](const TestTree::NearList::value_type& a, const TestTree::NearList::value_type& b) -> bool {
      if (a.second != b.second)
        return a.second < b.second;
      else {
        for (unsigned int d = 0; d < dimensionality; ++d) {
          if ((*a.first)[d] != (*b.first)[d])
            return (*a.first)[d] < (*b.first)[d];
        }

        // for the case where a is completely equal to b
        return false;
      }
    };

    //brute-force search
    TestTree::NearList bf_near_list, kd_near_list;
    for (const auto& candidate : all_nodes) {
      double dist_squared = 0.0;
      for (unsigned int d = 0; d < dimensionality; ++d) {
        double diff = centre[d] - (*candidate)[d];
        dist_squared += diff * diff;
      }

      if (dist_squared <= radius_squared) {
        bf_near_list.emplace_back(candidate, dist_squared);
      }
    }

    const bool using_coords(n % 2 == 0);
    if (using_coords)
      kd_near_list = tree.FindNears({centre[0], centre[1], centre[2]}, radius_squared);
    else
      kd_near_list = tree.FindNears(centre, radius_squared);
    bool mis_matched = false;
    if (bf_near_list.size() != kd_near_list.size()) {
      std::cout << WARNING_COLOUR << "Mismatch detected: (centre, radius^2) = (" << centre << "," << radius_squared << "), using_coords = " << using_coords << "; bf found " << bf_near_list.size() << " nodes; kd found " << kd_near_list.size() << " nodes" << COLOUR_END << std::endl;
      mis_matched = true;
    }
    else {

      bf_near_list.sort(comparator);
      kd_near_list.sort(comparator);

      for (auto bf_iter = bf_near_list.begin(), kd_iter = kd_near_list.begin(); 
            bf_iter != bf_near_list.end() && kd_iter != kd_near_list.end(); ++bf_iter, ++kd_iter) {
        for (unsigned int d = 0; d < dimensionality; ++d) {
          if ((*bf_iter->first)[d] != (*kd_iter->first)[d]) {
            std::cout << WARNING_COLOUR << "Mismatch detected: (centre, radius^2) = (" << centre << "," << radius_squared << "), using_coords =" << using_coords << "; bf and kd did not find the same list of nodes"<< COLOUR_END << std::endl;

            mis_matched = true;
            break;
          }
        }
        if (mis_matched)
          break;
      }
    }
    
    if (mis_matched) {
      mismatch_count += 1;
    }

  }

  if (mismatch_count)
    std::cout << ERROR_COLOUR << "mismatch_count = " << mismatch_count << COLOUR_END << std::endl;
  else
    std::cout << OK_COLOUR << "Pass" << COLOUR_END << std::endl;

}

void Test_FindNears(const unsigned int seed) {
  std::cout << NOTICE_COLOUR << __func__ << COLOUR_END << ": random seed = " << seed  << std::endl;
  std::minstd_rand random_engine(seed);
  std::uniform_real_distribution<double> coord_dice (-20.0, 20.0);
  constexpr unsigned int NUM_OF_NODES = 100;

  TestTree tree;
  TestTree::NodePtrList all_nodes;


  try {
    PrepareRandomTestTree(NUM_OF_NODES, random_engine, coord_dice, tree, all_nodes);
    std::cout << "Number of valid nodes = " << all_nodes.size() << std::endl;  
  }
  catch (const std::exception& e) {
    std::cout << ERROR_COLOUR << e.what() << COLOUR_END << std::endl;    
    return;
  }

  Test_FindNearsImpl(tree, all_nodes, NUM_OF_NODES, random_engine, coord_dice);


}


} // namespace unit_test

TreeBase::Region::Region(const unsigned int dimensionality)
: boundaries_(dimensionality) {
}

TreeBase::Region::Region(const Region& other)
: boundaries_(other.boundaries_.size()) {

  for (unsigned int d = 0; d < other.boundaries_.size(); ++d) {    
    const Boundary& other_boundary = other.boundaries_.at(d);
    if (other_boundary.first)
      boundaries_[d].first = std::make_shared<Boundary::first_type::element_type>(*other_boundary.first);
    if (other_boundary.second)
      boundaries_[d].second = std::make_shared<Boundary::second_type::element_type>(*other_boundary.second);
  }
}

TreeBase::Region::Region(Region&& other)
: boundaries_(std::move(other.boundaries_)) {
}

const TreeBase::Region::Boundary& TreeBase::Region::operator[](const unsigned int dimension) const {
  return boundaries_.at(dimension);
}

bool TreeBase::Region::Contain(const NodeInterface& node) const {

  for (unsigned int d = 0; d < boundaries_.size(); ++d) {
    const Boundary& boundary = boundaries_.at(d);
     
    if ((boundary.first && node[d] < *boundary.first) || (boundary.second && node[d] >= *boundary.second))
      return false;
  }

  return true;
}

bool TreeBase::Region::MightIntersect(const NodeInterface& centre, const double radius_squared, const bool inclusive) const {

  for(unsigned int d = 0; d < boundaries_.size(); ++d) {
    const Boundary& boundary = boundaries_.at(d);
    double dist = 0.0;
    
    if (boundary.first && centre[d] < *boundary.first)
      dist = *boundary.first - centre[d];
    else if (boundary.second && centre[d] >= *boundary.second)
      dist = centre[d] - *boundary.second;
    
    double dist_squared = dist * dist;
    if (!inclusive && dist_squared == radius_squared)
      return false;    
    else if (dist_squared > radius_squared)
      return false;

  }

  return true;

}


bool TreeBase::Region::BoundAboveBy(const TreeBase::SplitDesc& split) const {

  const Boundary& boundary = boundaries_.at(split.dimension);

  if (boundary.second && split.at >= *boundary.second)
    return true;
  else
    return false;

}

bool TreeBase::Region::BoundBelowBy(const TreeBase::SplitDesc& split) const {

  const Boundary& boundary = boundaries_.at(split.dimension);

  if (boundary.first && split.at <= *boundary.first)
    return true;
  else
    return false;

}


TreeBase::Region TreeBase::Region::Split(const TreeBase::SplitDesc& split) {
  Region left(*this);
  Boundary& right_boundary = boundaries_.at(split.dimension);
  Boundary& left_boundary = left.boundaries_.at(split.dimension);
  if (right_boundary.first)
    *right_boundary.first = split.at;
  else
    right_boundary.first = std::make_shared<double>(split.at);
  
  if (left_boundary.second)
    *left_boundary.second = split.at;
  else
    left_boundary.second = std::make_shared<double>(split.at);

  return left;
}

void TreeBase::Region::ExpandToCover(const unsigned int dimension, const Region::Boundary& boundary) {
  Boundary& current_boundary = boundaries_.at(dimension);

  if (!boundary.first)
    current_boundary.first.reset();
  else if (current_boundary.first && *current_boundary.first > *boundary.first)
    *current_boundary.first = *boundary.first;

  if (!boundary.second)
    current_boundary.second.reset();
  else if (current_boundary.second && *current_boundary.second < *boundary.second)
    *current_boundary.second = *boundary.second;

}

void TreeBase::Region::Print(std::ostream& os) const {
  os << "[";
  for (unsigned int d = 0; d < boundaries_.size(); ++d) {
    const Boundary& boundary = boundaries_[d];
    if (boundary.first)
      os << *boundary.first;
    else
      os << "MIN";

    os << ",";

    if (boundary.second)
      os << *boundary.second;
    else
      os << "MAX";

    if (d + 1 < boundaries_.size())
      os << ";";
  }
  os << "]";
}

std::ostream& operator<<(std::ostream& os, const TreeBase::Region& region) {
   region.Print(os);
   return os;
}


void TreeBase::Print(std::ostream& os) const {

  std::list<RegionRecord> root_list{{Region(dimensionality_), root_}};
  std::list<const std::list<RegionRecord>*> curr_list{&root_list}, next_list;
  os << NOTICE_COLOUR << "Tree: " << COLOUR_END << "K = " << dimensionality_ << ", B = " << collection_size_max_ << std::endl;
  os << NOTICE_COLOUR << "Node count = " << COLOUR_END << node_count_ << std::endl;
  while (!curr_list.empty()) {
    for (const std::list<RegionRecord>* sub_list : curr_list) {
      for (const RegionRecord& rr : *sub_list) {
        os << NOTICE_COLOUR << "  Region: " << rr.first << COLOUR_END;
        os << "<";
        if (rr.second->is_leaf) {        
          auto node_collection = std::static_pointer_cast<NodeCollection>(rr.second);
          for (const NodeInterfacePtr& node : node_collection->nodes) {
            os << *node;
          }
        }
        else {
          auto region_collection = std::static_pointer_cast<RegionCollection>(rr.second);
          next_list.emplace_back(&region_collection->regions);
        }
        os << ">";
      }
      os << NOTICE_COLOUR << " |" << COLOUR_END;
    }
    os << std::endl;
    curr_list.clear();
    curr_list.swap(next_list);
  }
  
}



bool TreeBase::Insert(const NodeInterfacePtr& new_node) {
  bool result = Insert(root_, new_node);
  if (result) {
    node_count_++;
    Region full(dimensionality_);
    RegionRecordPtr new_sub_region = SplitIfNeeded(full, root_);
    if (new_sub_region) {
      auto new_root = std::make_shared<RegionCollection>(false);
      
      new_root->regions.emplace_back(std::move(*new_sub_region));
      new_root->regions.emplace_back(std::move(full), std::move(root_));
      root_ = std::static_pointer_cast<Record>(new_root);
    }
  }

  return result;
}

TreeBase::RegionRecordPtr TreeBase::Split(Region& region, const RecordPtr& record, const SplitDesc& split) {


  if (record->is_leaf) {
    auto right_nc = std::static_pointer_cast<NodeCollection>(record);
    

    Region left_region = region.Split(split);

    auto left_nc = std::make_shared<NodeCollection>(true);
    for(NodeCollection::List::iterator node_iter = right_nc->nodes.begin(); node_iter != right_nc->nodes.end();) {
      if ((**node_iter)[split.dimension] < split.at) {
        left_nc->nodes.emplace_back(std::move(*node_iter));
        node_iter = right_nc->nodes.erase(node_iter);
      }
      else {
        ++node_iter;
      }
    }

    return std::make_shared<RegionRecord>(std::move(left_region), std::static_pointer_cast<Record>(left_nc));
  }
  else {
    auto right_rc = std::static_pointer_cast<RegionCollection>(record);
    
    Region left_region = region.Split(split);
    auto left_rc = std::make_shared<RegionCollection>(false);



    for (auto sub_region_iter = right_rc->regions.begin(); sub_region_iter != right_rc->regions.end();) {
      if (sub_region_iter->first.BoundAboveBy(split)) {
        left_rc->regions.emplace_back(std::move(*sub_region_iter));
        sub_region_iter = right_rc->regions.erase(sub_region_iter);
      }
      else {
        if (!sub_region_iter->first.BoundBelowBy(split)) {
          RegionRecordPtr sub_region_to_left = Split(sub_region_iter->first, sub_region_iter->second, split);
          if (sub_region_to_left) {
            left_rc->regions.emplace_back(std::move(*sub_region_to_left));
          }
        }

        ++sub_region_iter;
      }

    }

    return std::make_shared<RegionRecord>(std::move(left_region), std::static_pointer_cast<Record>(left_rc));
  
  }
}

TreeBase::RegionRecordPtr TreeBase::SplitIfNeeded(Region& region, const RecordPtr& record) {
  TreeBase::RegionRecordPtr new_rr;

  if (record->is_leaf) {
    auto node_collection = std::static_pointer_cast<NodeCollection>(record);
    if (node_collection->nodes.size() > collection_size_max_) {
      // Finding a coordinate at which a split can produce the most balanced node distribution,
      // i.e. the number of nodes in the left region is as close as possible to 
      // the number of nodes in the right region after the split

      SplitDesc split;
      unsigned int best_lr_diff = node_collection->nodes.size();
      unsigned int theoretical_best_lr_diff = node_collection->nodes.size() % 2;

      for (unsigned int d = 0; d < dimensionality_ && best_lr_diff > theoretical_best_lr_diff; ++d) {
        std::unordered_map<double, unsigned int> coord_counts;
        
        const Region::Boundary& boundary = region[d];
        for (const NodeInterfacePtr& existing_node : node_collection->nodes) {
          double coord = (*existing_node)[d];
          // Since no point sits at the upper boundary of regions, we can skip that part of checking
          if (!boundary.first || *boundary.first != coord) {

            std::pair<decltype(coord_counts)::iterator, bool> emplacing_result (coord_counts.emplace(coord, 1));
            if (!emplacing_result.second)
              emplacing_result.first->second += 1;
          }
        }

        if (coord_counts.size() > 1) {
          std::vector<std::pair<double, unsigned int>> candidate_points(coord_counts.begin(), coord_counts.end());
          unsigned int splitter_at = candidate_points.size() / 2;

          std::nth_element(candidate_points.begin(), candidate_points.begin()+splitter_at, candidate_points.end(), [](const decltype(candidate_points)::value_type& a, const decltype(candidate_points)::value_type& b) -> bool {
            return a.first < b.first;
          });

          unsigned int left = 0, lr_diff = 0;
          for (unsigned int p = 0; p < splitter_at; ++p) {
            left += candidate_points[p].second;
          }
          left *= 2;          
          if (left > node_collection->nodes.size())
            lr_diff = left - node_collection->nodes.size();
          else
            lr_diff = node_collection->nodes.size() - left;

          if (lr_diff < best_lr_diff) {
            best_lr_diff = lr_diff;
            split.dimension = d;
            split.at = candidate_points[splitter_at].first;
          }
        }
      }

      // If best_lr_diff >= node_collection->nodes.size(), 
      // it means a split point that produces a new region containing at least one node can not be found,
      // and this should not happen as there is no duplication
      if (best_lr_diff >= node_collection->nodes.size()) {
        throw std::runtime_error("Can not find a proper dimension to split; this might well mean tree data is corrupted!");
      }


      new_rr = Split(region, record, split);

    }

  }
  else {
    auto region_collection = std::static_pointer_cast<RegionCollection>(record);
    if (region_collection->regions.size() > collection_size_max_) {
      // Finding a boundary which borders the most sub-regions;
      // also boundaries that does not cut through any sub-regions when the containing region
      // is split at it are preferred to those does

      SplitDesc split;
      unsigned int max_shared_count = 0;
      char cutting_through = 0;

      for (unsigned int d = 0; d < dimensionality_ && max_shared_count < region_collection->regions.size(); ++d) {
        std::unordered_map<double, unsigned int> border_counts;
        decltype(border_counts)::iterator max_iter = border_counts.end();

        const Region::Boundary& boundary = region[d];
        for (const RegionRecord& rr : region_collection->regions) {

          if (cutting_through < 0) {
            // when cutting_through < 0, it means split has been set in the previous run,
            // and we check if that split could cut through any sub-regions in this run
            if (!rr.first.BoundAboveBy(split) && !rr.first.BoundBelowBy(split))
              cutting_through = 1;
          }

          const Region::Boundary& sub_boundary = rr.first[d];

          Region::Boundary checking_array[] = {
            {boundary.first, sub_boundary.first}, {boundary.second, sub_boundary.second}
          };


          for (const Region::Boundary& where : checking_array) {
            if (where.second && (!where.first || *where.first != *where.second)) {
              std::pair<decltype(border_counts)::iterator, bool> emplacing_result(border_counts.emplace(*where.second, 1));
              if (!emplacing_result.second)
                emplacing_result.first->second += 1;
              
              if (max_iter == border_counts.end() || max_iter->second < emplacing_result.first->second)
                max_iter = emplacing_result.first;
            }
          }
        }

        // split has been checked in the loop above, and no region is found being cut through by it;
        // thus setting cutting_through to 0
        if (cutting_through < 0)
          cutting_through = 0;
          
        if (max_iter != border_counts.end()) {
          if (max_shared_count < max_iter->second || (max_shared_count == max_iter->second && cutting_through == 1)) {
            max_shared_count = max_iter->second;
            split.dimension = d;
            split.at = max_iter->first;

            // Since split has been updated, setting cutting_through to -1 to trigger the cutting-through check in the next run
            cutting_through = -1;
          }
        }

      }

      if (max_shared_count < 2) {
        throw std::runtime_error("Can not find a proper dimension to split; this might well mean tree data is corrupted!");
      }

      new_rr = Split(region, record, split);
  
    }

  }

  return new_rr;
}


bool TreeBase::Insert(const RecordPtr& record, const NodeInterfacePtr& new_node) {

  if (record->is_leaf) {
    auto node_collection = std::static_pointer_cast<NodeCollection>(record);
    // node duplication (in terms of coordinates) is not allowed
    
    for (const NodeInterfacePtr& existing_node : node_collection->nodes) {
      bool duplicated = true;
      for (unsigned int d = 0; d < dimensionality_; ++d) {
        if ((*existing_node)[d] != (*new_node)[d]) {
          duplicated = false;
          break;
        }
      }
      if (duplicated)
        return false;
          
    }
    node_collection->nodes.emplace_back(new_node);
    
  }
  else {
    auto region_collection = std::static_pointer_cast<RegionCollection>(record);
    RegionCollection::List::iterator sub_region_iter;
    for (sub_region_iter = region_collection->regions.begin(); sub_region_iter != region_collection->regions.end(); ++sub_region_iter) {
      if (sub_region_iter->first.Contain(*new_node)) {

        if (!Insert(sub_region_iter->second, new_node))
          return false;
        
        RegionRecordPtr new_sub_region = SplitIfNeeded(sub_region_iter->first, sub_region_iter->second);
        if (new_sub_region) {        
          region_collection->regions.emplace_back(std::move(*new_sub_region));
        }

        break;
      }
    }

    if (sub_region_iter == region_collection->regions.end()) {
      std::ostringstream string_builder;
      string_builder << "Can not find a region to accommodate this node " << *new_node << "; this might well mean tree data is corrupted!"; 
      throw std::runtime_error(string_builder.str());
    }  
  }

  return true;

}

bool TreeBase::Delete(const NodeInterfacePtr& target_node) {
  Region full(dimensionality_);
  bool result = Delete(full, root_, target_node);
  if (result) {
    node_count_--;
    if (root_->Size() == 1 && !root_->is_leaf) {      
      auto region_collection = std::static_pointer_cast<RegionCollection>(root_);
      root_ = region_collection->regions.front().second;
    }
  }

  return result;
}


bool TreeBase::Delete(const Region& region, const RecordPtr& record, const NodeInterfacePtr& target_node) {
  if (record->is_leaf) {
    auto node_collection = std::static_pointer_cast<NodeCollection>(record);
    
    for (NodeCollection::List::iterator node_iter = node_collection->nodes.begin();
          node_iter != node_collection->nodes.end(); ++node_iter) {
      bool matched = true;
      for (unsigned int d = 0; d < dimensionality_; ++d) {
        if ((**node_iter)[d] != (*target_node)[d]) {
          matched = false;
          break;
        }
      }
      if (matched) {
        node_collection->nodes.erase(node_iter);
        return true;
      }
    }
  }
  else {
    auto region_collection = std::static_pointer_cast<RegionCollection>(record);
    for (RegionCollection::List::iterator sub_region_iter = region_collection->regions.begin(); sub_region_iter != region_collection->regions.end(); ++sub_region_iter) {
      if (sub_region_iter->first.Contain(*target_node)) {
        bool result = Delete(sub_region_iter->first, sub_region_iter->second, target_node);
        if (result) {
          if (sub_region_iter->second->Size() < collection_size_low_) {          
            if (region_collection->regions.size() > 1) {
              // Finding a boundary of the sub-region (sub_region_iter->first) which is shared
              // with the fewest neighbouring sub-regions, and joining them together
              
              SplitDesc split;
              std::list<RegionCollection::List::iterator> neighbour_list;
              // initialising min_shared_count slightly bigger than region_collection->regions.size(),
              // so that the following loop runs even if there are only 2 regions
              unsigned int min_shared_count = region_collection->regions.size()+1;

              for (unsigned int d = 0; d < dimensionality_ && min_shared_count > 2; ++d) {
                std::unordered_map<double, decltype(neighbour_list)> border_counts;

                const Region::Boundary& containing_boundary = region[d];
                const Region::Boundary& boundary = sub_region_iter->first[d];
                // Only boundaries inside the containing region need to be checked
                if (boundary.first && (!containing_boundary.first || *boundary.first != *containing_boundary.first)) {
                  border_counts[*boundary.first].emplace_back(sub_region_iter);
                }

                if (boundary.second && (!containing_boundary.second || *boundary.second != *containing_boundary.second)) {
                  border_counts[*boundary.second].emplace_back(sub_region_iter);
                }


                for (RegionCollection::List::iterator neighbour_iter = region_collection->regions.begin(); neighbour_iter != region_collection->regions.end(); ++neighbour_iter) {
                  if (neighbour_iter == sub_region_iter)
                    continue;

                  const Region::Boundary& neighbour_boundary = neighbour_iter->first[d];

                  if (neighbour_boundary.first) {
                    decltype(border_counts)::iterator search_iter = border_counts.find(*neighbour_boundary.first);
                    if (search_iter != border_counts.end())
                      search_iter->second.emplace_back(neighbour_iter);
                  }
                  if (neighbour_boundary.second) {
                    decltype(border_counts)::iterator search_iter = border_counts.find(*neighbour_boundary.second);
                    if (search_iter != border_counts.end())
                      search_iter->second.emplace_back(neighbour_iter);
                  }
                }

                for (decltype(border_counts)::value_type& border_pair : border_counts) {
                  if (border_pair.second.size() < min_shared_count) {
                    min_shared_count = border_pair.second.size();
                    neighbour_list.swap(border_pair.second);
                    split.dimension = d;
                    split.at = border_pair.first;
                  }
                }
              }

              if (neighbour_list.empty()) {
                throw std::runtime_error("Can not find a neighbour to join; this might well mean tree data is corrupted!");
              }

              RegionCollection::List::iterator joined_iter = neighbour_list.front();
              for (decltype(neighbour_list)::iterator to_be_joined_iter = ++neighbour_list.begin(); to_be_joined_iter != neighbour_list.end(); ++to_be_joined_iter) {
                joined_iter->first.ExpandToCover(split.dimension, (*to_be_joined_iter)->first[split.dimension]);
                joined_iter->second->Takeover((*to_be_joined_iter)->second);

                region_collection->regions.erase(*to_be_joined_iter);
              }

            }
          }

          // As the joining operation triggered by deletion might produce regions containing nodes/sub-regions of which the number exceeds collection_size_max_, 
          // calling SplitIfNeeded() to make sure they still conform to the requirement
          std::list<RegionCollection::List::iterator> might_need_split_list;
          might_need_split_list.emplace_back(sub_region_iter);
          while(!might_need_split_list.empty()) {
            const decltype(might_need_split_list)::value_type& rr_iter = might_need_split_list.front();
            
            RegionRecordPtr new_rr = SplitIfNeeded(rr_iter->first, rr_iter->second);
            if (new_rr) {
              region_collection->regions.emplace_back(std::move(*new_rr));
              might_need_split_list.emplace_back(--region_collection->regions.end());
            }
            else {
              might_need_split_list.pop_front();
            }
          }

        }

        return result;
      }
    }
  }

  return false;

}

void TreeBase::FindNearests(const NodeInterface& centre, const unsigned int k, std::multimap<double, NodeInterfacePtr>& nearests) const {
  nearests.clear();
  FindNearests(root_, centre, k, nearests);
}

void TreeBase::FindNearests(const RecordPtr& record, const NodeInterface& centre, const unsigned int k, std::multimap<double, NodeInterfacePtr>& nearests) const {
  if (record->is_leaf) {
    auto node_collection = std::static_pointer_cast<NodeCollection>(record);
    for (const NodeInterfacePtr& node : node_collection->nodes) {
      double dist_squared = 0.0;
      for (unsigned int d = 0; d < dimensionality_; ++d) {
        double diff = (*node)[d] - centre[d];
        dist_squared += diff * diff;
      }

      if (nearests.size() < k || dist_squared < (--(nearests.end()))->first) {
        nearests.emplace(dist_squared, node);
        if (nearests.size() > k) {
          nearests.erase(--(nearests.end()));
        }
      }
    }
  }
  else {
    auto region_collection = std::static_pointer_cast<RegionCollection>(record);
    RegionCollection::List::iterator sub_region_iter;
    for (sub_region_iter = region_collection->regions.begin(); sub_region_iter != region_collection->regions.end(); ++sub_region_iter) {          
      if (sub_region_iter->first.Contain(centre)) {        
        FindNearests(sub_region_iter->second, centre, k, nearests);
        break;
      }
    }

    for (RegionCollection::List::iterator other_iter = region_collection->regions.begin(); other_iter != region_collection->regions.end(); ++other_iter) {
      if (other_iter == sub_region_iter)
        continue;

      double checking_dist_squared = std::numeric_limits<double>::max();
      if (nearests.size() == k)
        checking_dist_squared = (--(nearests.end()))->first;

      if (other_iter->first.MightIntersect(centre, checking_dist_squared, false))
        FindNearests(other_iter->second, centre, k, nearests);          

    }

  }
}

void TreeBase::FindNears(const NodeInterface& centre, const double radius_squared, NearReportInterface& report) const {
  
  FindNears(root_, centre, radius_squared, report);
}


void TreeBase::FindNears(const RecordPtr& record, const NodeInterface& centre, const double radius_squared, NearReportInterface& report) const {
  if (record->is_leaf) {
    auto node_collection = std::static_pointer_cast<NodeCollection>(record);
    for (const auto& node : node_collection->nodes) {
      double dist_squared = 0.0;
      for (unsigned int d = 0; d < dimensionality_; ++d) {
        double diff = (*node)[d] - centre[d];
        dist_squared += diff * diff;
      }
      if (dist_squared <= radius_squared) {
        report(node, dist_squared);
      }
    }
  }
  else {
    auto region_collection = std::static_pointer_cast<RegionCollection>(record);
    RegionCollection::List::iterator sub_region_iter;
    for (sub_region_iter = region_collection->regions.begin(); sub_region_iter != region_collection->regions.end(); ++sub_region_iter) {          
      if (sub_region_iter->first.Contain(centre)) {        
        FindNears(sub_region_iter->second, centre, radius_squared, report);
        break;
      }
    }

    for (RegionCollection::List::iterator other_iter = region_collection->regions.begin(); other_iter != region_collection->regions.end(); ++other_iter) {
      if (other_iter == sub_region_iter)
        continue;

      if (other_iter->first.MightIntersect(centre, radius_squared, true))
        FindNears(other_iter->second, centre, radius_squared, report);

    }

  }
}


} // namespace kdb_tree
