#ifndef KDB_TREE_H
#define KDB_TREE_H

#include <iostream>
#include <memory>
#include <list>
#include <vector>
#include <unordered_map>
#include <map>
#include <chrono>

namespace kdb_tree {

class TreeBase;
class NodeInterface;
using NodeInterfacePtr = std::shared_ptr<NodeInterface>;

class NodeInterface {
  public:
    virtual double operator [](const unsigned int dimension) const = 0;
    virtual void Print(std::ostream& os) const {
    }


};

inline std::ostream& operator<<(std::ostream& os, const NodeInterface& node) {
  node.Print(os);
  return os;
}

class TreeBase {
  private:
    // forward declaration
    class NodeCollection;

  public:

    virtual void Print(std::ostream& os) const;

    TreeBase(const unsigned int dimensionality, const unsigned int collection_size)
    : dimensionality_(dimensionality),
      collection_size_max_(collection_size),
      collection_size_low_((collection_size+1)/2),
      root_(std::make_shared<NodeCollection>(true)),
      node_count_(0) {

    }

    virtual ~TreeBase() {
    };

    inline unsigned int NodeCount() const {
      return node_count_;
    }

    inline void Clear() {
      root_ = std::make_shared<NodeCollection>(true);
      node_count_ = 0;
    }



  protected:
    bool Insert(const NodeInterfacePtr& new_node);
    bool Delete(const NodeInterfacePtr& target_node);
    void FindNearests(const NodeInterface& centre, const unsigned int k, std::multimap<double, NodeInterfacePtr>& nearests) const;
    void FindNears(const NodeInterface& centre, const double radius_squared, std::multimap<double, NodeInterfacePtr>& nears) const;

    class SplitDesc {
      public:
        unsigned int dimension;
        double at;
    };


    class Region {
      public:
        using Boundary = std::pair<std::shared_ptr<double>, std::shared_ptr<double>>;

        Region(const unsigned int dimensionality);
        Region(const Region& other);
        Region(Region&& other);

        const Boundary& operator[](const unsigned int dimension) const;
        Region Split(const SplitDesc& split);
        void ExpandToCover(const unsigned int dimension, const Region::Boundary& boundary);

        bool Contain(const NodeInterface& node, const bool including_upper_bound = false) const;
        bool Intersect(const Region& other, const bool inclusive) const;
        bool Intersect(const NodeInterface& centre, const double radius_squared, const bool inclusive) const;
        bool BoundAboveBy(const SplitDesc& split) const;
        bool BoundBelowBy(const SplitDesc& split) const;

        void Print(std::ostream& os) const;


      protected:
        Region() = default;

        void SetBoundaries(std::vector<Boundary>&& boundaries);


      private:
        std::vector<Boundary> boundaries_;


    };
    // making the global function "operator<<" be the frend of TreeBase, so that it can see "Region" class
    friend std::ostream& operator<<(std::ostream& os, const Region& region);

    void FindContained(const Region& bounding_region, std::list<NodeInterfacePtr>& inside_list) const;

  private:


    class Record;
    using RecordPtr = std::shared_ptr<Record>;

    using RegionRecord = std::pair<Region, RecordPtr>;
    using RegionRecordPtr = std::shared_ptr<RegionRecord>;


    class Record {
      public:
        Record(const bool is_leaf): is_leaf(is_leaf) {
        }
        virtual void Takeover(const RecordPtr& other) = 0;
        virtual std::size_t Size() const = 0;

        bool is_leaf;
    };


    class NodeCollection : public Record {
      public:
        using List = std::list<NodeInterfacePtr>;

        NodeCollection(const bool is_leaf)
        : Record(is_leaf) {
        }

        inline virtual void Takeover(const RecordPtr& other) override {
          auto other_nc = std::static_pointer_cast<NodeCollection>(other);
          for (List::value_type& value : other_nc->nodes) {
            nodes.emplace_back(std::move(value));
          }
          other_nc->nodes.clear();
        }

        inline virtual std::size_t Size() const override {
          return nodes.size();
        }

        List nodes;

    };

    using NodeCollectionPtr = std::shared_ptr<NodeCollection>;

    class RegionCollection : public Record {
      public:
        using List = std::list<RegionRecord>;

        RegionCollection(const bool is_leaf)
        : Record(is_leaf) {
        }

        inline virtual void Takeover(const RecordPtr& other) override {
          auto other_rc = std::static_pointer_cast<RegionCollection>(other);
          for (List::value_type& value : other_rc->regions) {
            regions.emplace_back(std::move(value));
          }

          other_rc->regions.clear();
        }

        inline virtual std::size_t Size() const override {
          return regions.size();
        }

        List regions;
    };

    using RegionCollectionPtr = std::shared_ptr<RegionCollection>;


    RegionRecordPtr Split(Region& region, const RecordPtr& record, const SplitDesc& split);
    RegionRecordPtr SplitIfNeeded(Region& region, const RecordPtr& record);

    bool Insert(const RecordPtr& record, const NodeInterfacePtr& new_node);

    bool Delete(const Region& region, const RecordPtr& record, const NodeInterfacePtr& target_node);

    void FindNearests(const RecordPtr& record, const NodeInterface& centre, const unsigned int k, std::multimap<double, NodeInterfacePtr>& nearests) const;
    void FindNears(const RecordPtr& record, const NodeInterface& centre, const double radius_squared, std::multimap<double, NodeInterfacePtr>& nears) const;
    void FindContained(const RecordPtr& record, const Region& bounding_region, std::list<NodeInterfacePtr>& inside_list) const;

    const unsigned int dimensionality_;
    const unsigned int collection_size_max_;
    const unsigned int collection_size_low_;
    RecordPtr root_;
    unsigned int node_count_;
};


template <unsigned int CONSTANT_K, typename NODE_TYPE, unsigned int CONSTANT_B=3>
class Tree : public TreeBase {
  public:
    using Node = NODE_TYPE;
    using NodePtr = std::shared_ptr<Node>;
    using NodePtrList = std::list<NodePtr>;

    using NearInfo = std::pair<NodePtr, double>;
    using NearList = std::list<NearInfo>;

    static constexpr unsigned int K = CONSTANT_K;
    static constexpr unsigned int B = CONSTANT_B;

  private:
    class LocaterNode : public NodeInterface {
      public:
        LocaterNode(const std::initializer_list<double> coords)
        : coords_(coords) {
          if (coords_.size() != CONSTANT_K)
            throw std::logic_error("The number of coordinates is wrong");
        }
        virtual double operator [](const unsigned int dimension) const override {
          return coords_[dimension];
        }
        virtual void Print(std::ostream& os) const override {
          os << "(";
          for (unsigned int d = 0; d < coords_.size(); ++d) {
            os << coords_[d];
            if (d+1 < coords_.size())
              os << ",";
          }
          os << ")";
        }
      private:
        std::vector<double> coords_;
    };

    static inline NearList ConvertToNearList(const std::multimap<double, NodeInterfacePtr>& map) {
      NearList near_list;
      for (const auto& pair : map) {
        near_list.emplace_back(std::static_pointer_cast<Node>(pair.second), pair.first);
      }

      return near_list;

    }

    static inline NodePtrList ConvertToNodeList(const std::list<NodeInterfacePtr>& list) {
      NodePtrList node_ptr_list;
      for (const auto& each : list) {
        node_ptr_list.push_back(std::static_pointer_cast<Node>(each));
      }

      return node_ptr_list;

    }


  public:
    class BoundingRegion : public TreeBase::Region {
    public:
      BoundingRegion(const std::initializer_list<std::pair<double, double>> region_definition) {

        if (region_definition.size() != CONSTANT_K)
          throw std::logic_error("The number of boundaries is wrong");

        std::vector<Boundary> boundaries;
        boundaries.reserve(region_definition.size());
        for (const std::pair<double, double>& def : region_definition) {
          boundaries.emplace_back(std::make_shared<double>(def.first), std::make_shared<double>(def.second));
        }

        SetBoundaries(std::move(boundaries));

      }

      // Actually making the global function "operator<<" be the friend of BoundingRegion is redundant,
      // but here this is used as a trick to trigger the "name injection" rule of c++,
      // to solve the failure of template argument deduction that happens when defining "operator<<" with template arguments in the global namespace
      friend inline std::ostream &operator<<(std::ostream& os, const BoundingRegion& bounding_region) {
        bounding_region.Print(os);
        return os;
      }
    };

    Tree() : TreeBase(CONSTANT_K, CONSTANT_B) {
    }


    bool Insert(Node&& node) {
      return Insert(std::make_shared<Node>(std::move(node)));
    }

    bool Insert(const NodePtr& node) {
      return TreeBase::Insert(std::static_pointer_cast<NodeInterface>(node));
    }


    bool Delete(const std::initializer_list<double> coords) {
      return TreeBase::Delete(std::static_pointer_cast<NodeInterface>(std::make_shared<LocaterNode>(coords)));
    }

    bool Delete(const NodePtr& node) {
      return TreeBase::Delete(std::static_pointer_cast<NodeInterface>(node));
    }



    NearList FindNearests(const std::initializer_list<double> coords, const unsigned int k=1) const {
      std::multimap<double, NodeInterfacePtr> nearests;
      TreeBase::FindNearests(LocaterNode(coords), k, nearests);
      return ConvertToNearList(nearests);
    }


    NearList FindNearests(const Node& centre, const unsigned int k=1) const {
      std::multimap<double, NodeInterfacePtr> nearests;
      TreeBase::FindNearests(centre, k, nearests);
      return ConvertToNearList(nearests);
    }


    NearList FindNears(const std::initializer_list<double> coords, const double radius_squared) const {
      std::multimap<double, NodeInterfacePtr> nears;
      TreeBase::FindNears(LocaterNode(coords), radius_squared, nears);
      return ConvertToNearList(nears);
    }


    NearList FindNears(const Node& centre, const double radius_squared) const {
      std::multimap<double, NodeInterfacePtr> nears;
      TreeBase::FindNears(centre, radius_squared, nears);
      return ConvertToNearList(nears);
    }

    NodePtrList FindContained(const BoundingRegion& bounding_region) const {
      std::list<NodeInterfacePtr> inside_list;
      TreeBase::FindContained(bounding_region, inside_list);
      return ConvertToNodeList(inside_list);
    }


};

template <unsigned int CONSTANT_K, typename NODE_TYPE, unsigned int CONSTANT_B>
inline std::ostream &operator<<(std::ostream& os, const Tree<CONSTANT_K, NODE_TYPE, CONSTANT_B>& tree) {
     tree.Print(os);
     return os;
}


namespace unit_test {
  void Test_InsertThenDelete(const unsigned int num_of_nodes=10, const unsigned int seed=std::chrono::system_clock::now().time_since_epoch().count());
  void Test_FindNearests(const unsigned int k=1, const unsigned int seed=std::chrono::system_clock::now().time_since_epoch().count());
  void Test_FindNears(const unsigned int seed=std::chrono::system_clock::now().time_since_epoch().count());
  void Test_FindContained(const unsigned int seed=std::chrono::system_clock::now().time_since_epoch().count());

} // namespace unit_test

} // namespace kdb_tree

#endif
