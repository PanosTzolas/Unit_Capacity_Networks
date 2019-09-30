#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

template<typename nTy = unsigned, typename wTy = unsigned>
struct Graph {

	typedef nTy NodeType;
	typedef wTy WeightType;

	struct NodeData {
		NodeData(){}
		NodeData(const NodeType& n)
			: nodeId(n){}
		
		bool operator < (const NodeData& other) const{
			return this->nodeId < other.nodeId ? true : false;
		}

		bool operator == (const NodeData& other) const{
			return this->nodeId == other.nodeId ? true : false;
		}		


		NodeType nodeId;
		
	};

	struct EdgeProperties{
		EdgeProperties(){}
		EdgeProperties(const WeightType& w)
			: constWeight(w){}
		WeightType constWeight;

	};

	struct EdgeData : public EdgeProperties {
		EdgeData(){}
		EdgeData(const NodeData& n)
			: node(n){}
		
		EdgeData(const NodeData& n, const WeightType& w)
			: EdgeProperties(w), node(n){} 


		bool operator < (const EdgeData& other) const{
			return this->node < other.node ? true : false ; 
		}

		bool operator == (const EdgeData& other) const{
			return this->node == other.node ? true : false;
		}

		friend ostream& operator << (ostream& os, const EdgeData& edge){
			return os << edge.node.nodeId;
		}
		
		NodeData node;
	};


    Graph()
        : num_nodes(0), num_edges(0)  {}

    Graph(const unsigned& n, const unsigned& m)
        : num_nodes(n), num_edges(m) {
        Edges.resize(n);
    }

    void generateRandomGraph() {
        unsigned counter = 0;
        do {
            for(unsigned i = 0; i < num_nodes; ++i) {
                for(unsigned j = 0; j < getRandomNumber(); ++j) {
                    EdgeData item(getRandomNode(), getRandomWeight());
                   Edges[i].push_back(item);
                    counter++;
                    if(counter == num_edges)break;
                }
                std::sort(Edges[i].begin(), Edges[i].end());
                unsigned previousSize = Edges[i].size();
                Edges[i].erase( std::unique(Edges[i].begin(), Edges[i].end()),
                        Edges[i].end());
                counter -= (previousSize - Edges[i].size());
                if(counter == num_edges)break;
            }
        } while (counter != num_edges);
    }

    void printEdges() {
        unsigned counter = 0;
        for(unsigned i = 0; i < num_nodes; ++i) {
            std::vector<EdgeData>& curSet = Edges[i];
            for(auto itr = curSet.begin(); itr != curSet.end(); ++itr) {
                const EdgeData& item = *itr;
                std::cout << counter++ << ". ( " << i << ", " << item << " )" << 
                    "\n"; 
            }
        }
    }
    
    unsigned getRandomNumber() { return std::rand() % 500; }
    NodeType getRandomNode(){ return std::rand() % num_nodes; }
    WeightType getRandomWeight(){return std::rand()%100; }

    unsigned num_nodes;
    unsigned num_edges;
    std::vector<std::vector<EdgeData>> Edges;
};





#endif
