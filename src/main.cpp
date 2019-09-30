#include <iostream>

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include <boost/property_map/property_map.hpp>
#include <boost/array.hpp>
#include <boost/graph/graphviz.hpp>

#include <boost/heap/priority_queue.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/lockfree/queue.hpp>

#include <limits> 
#include <stdlib.h> 
#include <cstdlib>
#include <list> 
#include <array>
#include <vector>

#include <math.h>
#include <time.h> 
#include <ctime>

#include <boost/config.hpp>
#include <utility>

#include "LEDA/graph/graph.h"
#include "LEDA/core/random_source.h"
#include "LEDA/core/p_queue.h"
#include "LEDA/graph/node_partition.h"
#include "LEDA/graph/node_list.h"
#include <LEDA/core/dictionary.h>
#include <LEDA/graph/node_pq.h>
#include <LEDA/graph/node_array.h>
#include <LEDA/graph/max_flow.h>
#include <string>


#define UNIT_CAPACITY 1
#define GAP 20
#define GAP_2 1000

using namespace std;
using namespace boost;
using namespace boost::heap;



typedef property<edge_weight_t, int,
		property<edge_index_t, int,
		property<edge_capacity_t, int,
		property<edge_residual_capacity_t, int,
		property<edge_reverse_t, int> > > > > EdgeProperties;

typedef adjacency_list_traits<vecS, vecS, directedS> Traits;

typedef adjacency_list_traits<vecS, vecS,		// PROSOXH ISWS listS
	directedS>::vertex_descriptor vertex_descriptor;

typedef adjacency_list_traits<vecS, vecS,		// isws listS
	directedS>::edge_descriptor edge_descriptor;

typedef pair<vertex_descriptor, vertex_descriptor> my_e;
typedef std::list <my_e> klist;



typedef property <vertex_name_t, klist,	// here int
	property<vertex_index2_t, int, 
	property<vertex_distance_t, int,
	property<vertex_predecessor_t, vertex_descriptor,		/*  Vertex or swich to vertex_descriptor  */
	property<vertex_color_t, int> > > > >VertexProperties;

typedef adjacency_list< vecS, vecS,
	bidirectionalS,							/* PROSOXH GIA NA XRHSIMOPOIOUME TO in_edges() alliws, directedS */
	VertexProperties,
	EdgeProperties
> Graph;


typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::edge_iterator edge_iterator;



typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertex_iterator vertex_iter;



typedef pair<int, Vertex> pi;
typedef pair<Vertex, Edge> pve;

/*
	This function creates a basic graph to start.
*/


void create_graph(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);


	std::pair<vertex_iter, vertex_iter> vp;

	vp = vertices(g);

	s = *vp.first;  // node start s
	
	add_edge(*vp.first, *vp.first + 1, 1, g);
	add_edge(*vp.first, *vp.first + 2, 1, g);
	color[*vp.first] = 0; // set color of all nodes to 0
	vp.first++;
	color[*vp.first] = 0;
	add_edge(*vp.first, *vp.first + 2, 1, g);
	add_edge(*vp.first, *vp.first + 3, 1, g);
	vp.first++;
	color[*vp.first] = 0;
	add_edge(*vp.first, *vp.first + 1, 1, g);
	add_edge(*vp.first, *vp.first + 2, 1, g);
	vp.first++;
	color[*vp.first] = 0;
	add_edge(*vp.first, *vp.first + 2, 1, g);
	vp.first++;
	color[*vp.first] = 0;
	add_edge(*vp.first, *vp.first + 1, 1, g);

	vp.first++;
	color[*vp.first] = 0;
	t = *vp.first; // node sink t

	std::ofstream gout;
	gout.open("test.txt");
	boost::write_graphviz(gout, g);


}

/*
	This function create the Graph of Figure 7.4 , page 231 
*/

void create_figure_graph(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(g);
	s = *vp.first;
	
	int i = 0;
	vp.first++;
	for (i = 1; i <= GAP; i++) {
		add_edge( s, *vp.first, 1, g);
		add_edge(*vp.first, *vp.first + GAP, 1, g);
		vp.first++;
	}

	for (i = 1; i <= GAP; i++) {
		add_edge(*vp.first, *vp.first + GAP + 1 - i, 1, g);
		vp.first++;
	}

	t = *vp.first;
	//cout << "t:\t" << t << endl;

	std::ofstream gout;
	gout.open("test.txt");
	boost::write_graphviz(gout, g);
}

void create_figure_2(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);


	int scaling= 5;

	std::pair<vertex_iter, vertex_iter> vp;

	vp = vertices(g);

	t = *vp.first;
	vp.first++;;
	s = *vp.first;


	int i = 0;
	vp.first++;
	for (i = 0; i <= GAP; i++) {
		add_edge(s, *vp.first, 1, g);
		add_edge(*vp.first, *vp.first + 1, 1, g);
		vp.first++;
		for (int k = 0; k <= 5; k++) {
			add_edge(*vp.first, *vp.first + 1, 1, g);
			vp.first++;
			add_edge(*vp.first, *vp.first + 1, 1, g);
			vp.first++;
		}

		if (i == 4) {
			add_edge(*vp.first, *vp.first + 1, 1, g);
			vp.first++;
		}
		add_edge(*vp.first, t, 1, g);
		vp.first++;
	}

	


	std::ofstream gout;
	gout.open("test2.txt");
	boost::write_graphviz(gout, g);

}

void create_figure_3(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);


	int scaling = 5;
	

	std::pair<vertex_iter, vertex_iter> vp;

	vp = vertices(g);

	t = *vp.first;
	vp.first++;;
	s = *vp.first;

	

	int i = 0;
	vp.first++;
	for (i = 0; i <= GAP; i++) {
		add_edge(s, *vp.first, 1, g);
		add_edge(*vp.first, *vp.first + 1, 1, g);
		vp.first++;
		for (int k = 0; k <= GAP_2; k++) {
			add_edge(*vp.first, *vp.first + 1, 1, g);
			vp.first++;
			add_edge(*vp.first, *vp.first + 1, 1, g);
			vp.first++;
		}
		add_edge(*vp.first, t, 1, g);
		vp.first++;
	}


	std::ofstream gout;
	gout.open("test3.txt");
	boost::write_graphviz(gout, g);

}
void generate_random_graph(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);
	srand(time(NULL));
	int gen_num = rand() % 100;
	std::list<Vertex> myL;
	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(g);

	for (vp = vertices(g); vp.first != vp.second; ++vp.first) {

		gen_num = rand() % 100;

		if (gen_num % 2 == 0) {
			myL.push_back(*vp.first);
		}
	}
	Vertex temp;
	int count = 1;
	
	for (auto it = myL.begin(); it != myL.end(); ++it) {
		if (count % 2 == 0) {
			if (temp != t && *it != s) {
				add_edge(temp, *it, 1, g);
			}
		}
		temp = *it;
		count++;
	}

	cout<<"Num of edges =\t"<< num_edges(g)<<endl;

	std::ofstream gout;
	gout.open("test4.txt");
	boost::write_graphviz(gout, g);

}

/*
	This function set capacities to 0 and flows to 0 of all edges
*/

 void set_vertices(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(g);
	for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
		color[*vp.first] = 0;
	}

}

void set_edges(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	std::pair<edge_iterator, edge_iterator> ei = edges(g);

	//std::cout << "Number of edges= " << num_edges(g) << endl;
	//std::cout << "Edge list:\n";

	int starting_flow = 0;
	int starting_capacity = UNIT_CAPACITY;

	/*  print all edges  */
	for (edge_iterator it = ei.first; it != ei.second; ++it) {

		c[*it] = starting_capacity;
		f[*it] = starting_flow;
		r[*it] = c[*it] - f[*it];
		//cout << *it;
		//cout << "\t" << edge_w[*it] << endl;
		//cout << "capacity: \t" << c[*it] << "\t flow: \t" << f[*it] << endl;
		//cout << " Target: \t" << boost::target(*it, g) << endl;
	}



}
void print_vertices(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(g);
	for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
		//color[*vp.first] = 0;
		cout << "color[" << *vp.first << "] =\t" << color[*vp.first] << endl;
	}


}


void my_leda_alg(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);


	leda::graph G1;

	int num_of_nodes = (int)num_vertices(g);
	
	
	leda::node v1[num_of_nodes+2];


	for(int i = 0; i < num_of_nodes; i++){
		v1[i] = G1.new_node();
	}


	int kappa = index[s];

	leda::node start_s = v1[index[s]];
	leda::node sink_t = v1[index[t]];

	//leda::edge_array<int> cap(G1);
	
	std::pair<edge_iterator, edge_iterator> ei = edges(g);
	for (edge_iterator it = ei.first; it != ei.second; ++it) {
			G1.new_edge(v1[index[source(*it,g)]], v1[index[target(*it,g)]]);
	}
	leda::edge_array<int> flow(G1);
	leda::edge_array<int> cap(G1);


	int klap = 1;	
	leda::edge my_edge;
	forall_edges(my_edge,G1){
		cap[my_edge] = (int) UNIT_CAPACITY;
	}

	clock_t start_leda, stop_leda;	
	double duration;

	start_leda = clock();
	int flow_value = leda::MAX_FLOW(G1, start_s, sink_t, cap, flow);
	stop_leda = clock();
	duration = (stop_leda - start_leda)/(double)CLOCKS_PER_SEC;

	cout<< "flow =\t" << flow_value << endl;
	cout<< "duration leda =\t"<< duration << endl;

		
}




void my_BFS(int &flow, Graph &g, Vertex &s, Vertex &t, int dbound, std::map<int, int> &numb,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);


	std::queue<Vertex> mQ;

	std::pair<vertex_iter, vertex_iter> kp;
	kp = vertices(g);

	color[s] = 1;  // mark s
	
	int lvl = 0;
	Vertex u;
	mQ.push(s);
	layer[s] = lvl;
	//u = mQ.back();

	for (int k = 0; k <= dbound; k++) {
		numb[k] = 0;
	}

	int k = 1;
	while (!mQ.empty()) {
		k++;

		u = mQ.front();
		mQ.pop();

		graph_traits<Graph>::out_edge_iterator ei, ei_end;

		for (boost::tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {

			if (color[target(*ei, g)] == 1) {
				color[target(*ei, g)] = 0;
				mQ.push(target(*ei, g));
				//my_map[d[target(*ei, g)]].push_back(target(*ei, g));
				numb[d[target(*ei, g)]] = numb[d[target(*ei, g)]] + 1;
			}
		}
	}

	
	
}




void reverse_BFS(Graph &g, Vertex &s, Vertex &t, int dbound,  std::map<int, int> &numb,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	Vertex v;
	std::queue<Vertex> Q;
	Q.push(t);
	color[t] = 1;
	d[t] = 0;

	
	while (!Q.empty()) {
		v = Q.front();
		//cout << "v=\t" << v << endl;
		
		Q.pop();
		boost::graph_traits<Graph>::in_edge_iterator ei, edge_end;
		for (boost::tie(ei, edge_end) = in_edges(v, g); ei != edge_end; ++ei) {
			
			
			if (color[source(*ei, g)] == 0) {
				color[source(*ei, g)] = 1;
				Q.push(source(*ei, g));
				d[source(*ei, g)] = d[v] + 1;
				numb[d[source(*ei, g)]] = numb[d[source(*ei, g)]] + 1;
				//cout << "\tcolor =\t"<< color[source(*ei, g)] << endl;
				//cout << "\td[" << source(*ei, g) << "] =\t" << d[source(*ei, g)] << endl;
			}
			else if (color[source(*ei, g)] == 1) {
				//cout << "\tcolor =\t" << color[source(*ei, g)] << endl;
			}
			//cout << "\t v = \t" << v  << "\t source:\t"<<source(*ei, g) << endl;
			//cout << "\t Q last:\t" << Q.back() << endl;

		}
		//cout << endl;
		//if (color[s] == 1) break;
	}

	//cout << "\t End of reverse_BFS" << endl;
}


void advance(Vertex &v, Vertex &u, my_e &pe, Edge &e,Graph &g, Vertex &s, Vertex &t, std::map<int, std::list<Edge>> &my_map,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	/*	let (i, j) be an admissible arc in A(i);  */

	//cout << "\tadvance(" << v << ")" << endl;

	my_map[index[v]].push_back(e);

	x[v].push_back(pe);

	//x[v] = e;
	//A[v] = e;

	p[u] = v;			//	pred(j) : = i
	//cout << "p[" << u << "] = \t" << v << "\t u =\t" << u << endl;
	//	tha xreiastei isws kai h akmh
	v = u;		//	i: = j;


}

void retreat(Vertex &v, Vertex &u, Edge &e, Graph &g, Vertex &s, Vertex &t, std::map<int, std::list<Edge>> &my_map,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	//cout << "\tretreat(" << v << ")" << endl;

	/*       d(i) : = min{d(j) + 1 : (i, j) E A(/) and rij> O};             */
	
	bool flag_end_a = false;
	std::list<Edge> helper;

	std::pair<Edge, bool> my_edge_pair;
	int i = 0;
	bool flag = true;
	list<int> mylist{ 1, 2, 3, 4, 5 };

	// using begin() to print list 
	Vertex min_v;
	for (auto it = my_map[index[v]].begin(); it != my_map[index[v]].end(); ++it) {
		if (source(*it, g) == v) {
			if (r[*it] > 0 && flag==true) {
				flag = false;
				min_v = target(*it, g);
			}
			if (flag == false) {
				if (r[*it] > 0 && d[target(*it, g)] < d[min_v]) {
					min_v = target(*it, g);
				}	
			}
		}
		else {
			if (r[*it] > 0 && flag == true) {
				flag = false;
				min_v = source(*it, g);
			}
			if (flag == false) {
				if (r[*it] > 0 && d[source(*it, g)] < d[min_v]) {
					min_v = source(*it, g);
				}
			}
		}
	}
	
	if (flag == false) {
		d[v] = d[min_v] + 1;
	}
	else {
		d[v] = d[v] + 1;
	}
	
	//d[v] = d[u] + 1;

	/* if(i!=s) then i := pred(i)	*/	

	if (v != s) {
		v = p[v];
		//cout << "if on retreat()" << endl;
	}

	
}

void retreat_b(Vertex &v, Vertex &u, Edge &e, Graph &g, Vertex &s, Vertex &t, int &vary, int dbound, std::map<int, int> &numb, std::map<int, std::list<Edge>> &my_map,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	//cout << "\tretreat(" << v << ")" << endl;

	/*       d(i) : = min{d(j) + 1 : (i, j) E A(/) and rij> O};             */

	bool flag = true;

	std::pair<Edge, bool> my_edge_pair;
	int i = 0;
	Vertex min_v, temp;
	for (auto it = my_map[index[v]].begin(); it != my_map[index[v]].end(); ++it) {
		if (source(*it, g) == v) {
			if (r[*it] > 0 && flag == true) {
				flag = false;
				min_v = target(*it, g);
			}
			if (flag == false) {
				if (r[*it] > 0 && d[target(*it, g)] < d[min_v]) {
					min_v = target(*it, g);
				}
			}
		}
		else {
			if (r[*it] > 0 && flag == true) {
				flag = false;
				min_v = source(*it, g);
			}
			if (flag == false) {
				if (r[*it] > 0 && d[source(*it, g)] < d[min_v]) {
					min_v = source(*it, g);
				}
			}
		}

	}
	int temp1, temp2;
	if (flag == false) {
		vary = d[v];
		numb[d[v]] = numb[d[v]] - 1;
		numb[d[min_v]] = numb[d[min_v]] + 1;
		d[v] = d[min_v] + 1;
	}
	else {
		vary = d[v];
		temp1 = d[v];
		numb[d[v]] = numb[d[v]] - 1;
		numb[temp1] = numb[temp1] + 1;
		d[v] = d[v] + 1;
	}

	

	//d[v] = d[u] + 1;


	/* if(i!=s) then i := pred(i)	*/

	if (v != s) {
		v = p[v];
		//cout << "if on retreat()" << endl;
	}


}


void augment2(int &flow,Vertex &v, Vertex &u, Graph &g, Vertex &s, Vertex &t, 
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	Vertex h1;

	h1 = t;

	while (h1 != s) {
		boost::graph_traits<Graph>::in_edge_iterator ei, edge_end;
		for (boost::tie(ei, edge_end) = in_edges(h1, g); ei != edge_end; ++ei) {
			if (source(*ei, g) == p[target(*ei, g)]) {
				f[*ei] = f[*ei] + c[*ei];
				r[*ei] = c[*ei] - f[*ei];
				h1 = source(*ei, g);
				break;
			}
		
		}

	}
	flow = flow + 1;

}

void phase_b(int &flow, Graph &g, Vertex &v, Vertex &s, Vertex &t, int dbound, std::map<int, int> &numb, std::map<int, std::list<Edge>> &my_map,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
) {
	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	bool my_flag = true;
	int vary = 1;
	//my_BFS(flow, g, s, t, dbound, numb, x , layer, d, p, color, edge_w, f, c, r, reverse);

	int my_d = (int) GAP * (int) GAP_2;	

	Vertex u;
	while (my_flag == true && d[s] <= my_d) {
		

		graph_traits<Graph>::out_edge_iterator ei, ei_end;
		Edge e;
		my_e pe;
		bool FLAG = true;
		for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {		//	if i has an admissible arc then

			u = target(*ei, g);
			FLAG = true;
			if ((d[v] == d[u] + 1) && (r[*ei] > 0)) {	//d[v] == d[u] + 1 &&  f[*ei] == 0 && r[*ei] > 0
				e = *ei;
				pe = make_pair(v, u);
				advance(v, u, pe, e, g, s, t, my_map, x, layer, d, p, color, edge_w, f, c, r, reverse);
				if (v == t) {
					augment2(flow, v, u, g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);
					v = s;
				}
				FLAG = false;
				break;
			}

		}
		if (FLAG == true) {
			retreat_b(v, u, e, g, s, t, vary, dbound, numb, my_map, x, layer, d, p, color, edge_w, f, c, r, reverse);
			FLAG = false;
		}

		if (numb[vary] < 0) {
			my_flag = false;
			break;
		}
		
	}
	
}


void shortest_augmenting_path_alg(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &layer,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w,
	property_map<Graph, edge_index_t>::type &f,
	property_map<Graph, edge_capacity_t>::type &c,
	property_map<Graph, edge_residual_capacity_t>::type &r,
	property_map<Graph, edge_reverse_t>::type &reverse
	) {
	
	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	std::pair<vertex_iter, vertex_iter> kp;
	kp = vertices(g);

	std::map<int, std::list<Edge>> my_map;

	std::queue <pve> A;
	
	int flow = 0;	// x:= 0

	int n = (int)num_vertices(g);
	int m = (int)num_edges(g);
	double  my_val = (double)2 / 3;
	double my_val2 = (double)1 / 2;
	double f1 = (double)pow(n, my_val) * 2;
	double f2 = (double)pow(m, my_val2);
	double f1_c = ceil(f1);
	double f2_c = ceil(f2);
	int d_bound = (int)min(f1_c, f2_c);
	int vary = 1;

	std::map<int, int> numb;

	reverse_BFS(g, s, t, d_bound, numb, x, layer, d, p, color, edge_w, f, c, r, reverse);	//  obtain the exact distance labels

	Vertex v,u;	
	v = s;			// i:= s
	
	int debug_counter = 0;

	//cout << "my_val = \t"<< my_val <<"\tn =\t"<< n <<"\tf1_c =\t"<< f1_c<<"\tf1 =\t" << f1 << endl;
	//cout << "mu_val2 =\t" << my_val2 << "\t m =\t" << m <<"\tf2_c =\t"<< f2_c<< " \tf2 =\t" << f2 << endl;

	bool my_flag = true;
	
	
	

	//cout << "d_bound = \t" << d_bound << endl;

		/*		d[s] <= d_bound && my_flag==true		*/		
	while (d[s] <= d_bound ) {		//	while d(s) < n do
		
		graph_traits<Graph>::out_edge_iterator ei, ei_end;
		Edge e;
		my_e pe;
		bool FLAG = true;
		for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {		//	if i has an admissible arc then
			
			u = target(*ei, g);
			FLAG = true;
			if ( (d[v] == d[u] + 1) && (r[*ei] > 0)) {	//d[v] == d[u] + 1 &&  f[*ei] == 0 && r[*ei] > 0
				e = *ei;
				
				pe = make_pair(v, u);
				advance(v, u, pe, e, g, s, t, my_map, x, layer, d, p, color, edge_w, f, c, r, reverse);
				if (v == t) {
					augment2(flow, v, u, g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);
					v = s;
					
				}
				FLAG = false;
				break;
			}



		}
		if (FLAG == true) {
			retreat_b(v, u, e, g, s, t, vary, d_bound, numb, my_map, x, layer, d, p, color, edge_w, f, c, r, reverse);
			//retreat(v, u, e, g, s, t, my_map, x, layer, d, p, color, edge_w, f, c, r, reverse);
			FLAG = false;
		}
		
	}
	
	phase_b(flow, g, v, s, t, d_bound, numb, my_map, x, layer, d, p, color, edge_w, f, c, r, reverse);


	
	cout << "Flow =\t" << flow << endl;

}


int main(int argc, char** argv) {

	Graph g;

	property_map<Graph, vertex_name_t>::type
		x = get(vertex_name, g);
	property_map<Graph, vertex_index2_t>::type
		layer = get(vertex_index2, g);
	property_map<Graph, vertex_distance_t>::type
		d = get(vertex_distance, g);
	property_map<Graph, vertex_predecessor_t>::type
		p = get(vertex_predecessor, g);
	property_map<Graph, vertex_color_t>::type
		color = get(vertex_color, g);
	property_map<Graph, edge_weight_t>::type
		edge_w = get(edge_weight, g);
	property_map<Graph, edge_index_t>::type
		f = get(edge_index, g);
	property_map<Graph, edge_capacity_t>::type
		c = get(edge_capacity, g);
	property_map<Graph, edge_residual_capacity_t>::type
		r = get(edge_residual_capacity, g);
	property_map<Graph, edge_reverse_t>::type
		reverse = get(edge_reverse, g);


	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	/*   Create Graph  -- Create basic graph to start    */

	Vertex s, t;
	//create_graph(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);

	//create_figure_graph(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);

	//create_figure_2(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);

	create_figure_3(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);

	generate_random_graph(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);

	//graph_comp(g, s, t);

	/*  Create Graph */

	/* SET VERTICES AND EDGES  */

	set_vertices(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);
	set_edges(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);

	/* SET EDGES  */



	//my_BFS(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);
	//reverse_BFS(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);

	std::clock_t start, stop;
	double duration;

	start = std::clock();
	shortest_augmenting_path_alg(g, s, t, x, layer, d, p, color, edge_w, f, c, r, reverse);
	stop = std::clock();
	duration = (stop - start) / (double)CLOCKS_PER_SEC;

	cout << "TIME: \t" << duration << endl;

	my_leda_alg(g, s, t, x, layer, d, p, color, edge_w, f, c, r);


	return 0;
}
