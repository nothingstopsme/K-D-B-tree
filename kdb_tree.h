#ifndef KDB_TREE_H
#define KDB_TREE_H

#include <iostream>
#include <memory>
#include <list>
#include <vector>
#include <limits>
#include <unordered_set>
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

    class NearReportInterface {
      public:
        virtual void operator()(const NodeInterfacePtr& near, const double dist_squared) = 0;
    };
    void FindNears(const NodeInterface& centre, const double radius_squared, NearReportInterface& report) const;

  private:
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
        
        bool Contain(const NodeInterface& node) const;
        bool MightIntersect(const NodeInterface& centre, const double radius_squared, const bool inclusive) const;
        bool BoundAboveBy(const SplitDesc& split) const;
        bool BoundBelowBy(const SplitDesc& split) const;

        void Print(std::ostream& os) const;

      private:
        std::vector<Boundary> boundaries_;
        
    };
    friend std::ostream& operator<<(std::ostream& os, const Region& region);

    

    
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
    void FindNears(const RecordPtr& record, const NodeInterface& centre, const double radius_squared, NearReportInterface& report) const;

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

    class NearReport : public NearReportInterface {
      public:
        inline virtual void operator()(const NodeInterfacePtr& near, const double dist_squared) override {
          list.emplace_back(std::static_pointer_cast<Node>(near), dist_squared);
        }

        NearList list;
    };


  public:

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
      NearList near_list;
      for (const auto& pair : nearests) {
        near_list.emplace_back(std::static_pointer_cast<Node>(pair.second), pair.first);
      }
      return near_list;

    }
    
    
    NearList FindNearests(const Node& centre, const unsigned int k=1) const {
      std::multimap<double, NodeInterfacePtr> nearests;
      TreeBase::FindNearests(centre, k, nearests);
      NearList near_list;
      for (const auto& pair : nearests) {
        near_list.emplace_back(std::static_pointer_cast<Node>(pair.second), pair.first);
      }
      return near_list;

    }

   
    NearList FindNears(const std::initializer_list<double> coords, const double radius_squared) const {      
      NearReport report;
      TreeBase::FindNears(LocaterNode(coords), radius_squared, report);
      return report.list;
    }

    
    NearList FindNears(const Node& centre, const double radius_squared) const {
      NearReport report;
      TreeBase::FindNears(centre, radius_squared, report);
      return report.list;
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

} // namespace unit_test

} // namespace kdb_tree

#endif
