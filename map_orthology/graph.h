#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include "gen_dna_funcs.h"
#include <queue>
#include <map>

using namespace::std;

#ifndef ___GRAPH_H___
#define ___GRAPH_H___

#ifndef GCC_COMPILE

#ifdef TRUE
#undef TRUE
#endif

#ifdef FALSE
#undef FALSE
#endif
#endif
//enum BOOL {FALSE, TRUE};

enum GENE_TYPE {TRANS_FACT, REG_GENE, DUMMY};

struct Node_dist {
	int dist;
	BOOL connected;
};

template <class LIST_UNIT>
class Linked_list {
	public:
		Linked_list *next, *last;

		Linked_list()   {};
		Linked_list(Linked_list *new_last, LIST_UNIT *new_element);
		LIST_UNIT* get_element()      {return(this_element);};
	~Linked_list();
	
	protected:
		LIST_UNIT *this_element;
};

template <class NODE_UNIT>
class Edge;

template <class NODE_UNIT>
class Node {
public:
	Edge<NODE_UNIT> *best_edge;
	double best_edge_score;
	
	Node();
	Node(BOOL undir);
	Node(int num, NODE_UNIT *node_element);
	Node(int num, NODE_UNIT *node_element, BOOL undir);
	int get_num_edges()         {return(num_edges);};
	void reset_edges()                  {trans_edges=start_edges;};
	void initialize_path_lens(int num_nodes);
	Node_dist* get_dist_to_node(int node_num);
	void set_dist_to_node(int node_num, int dist, BOOL connected)    {path_lengths[node_num].dist=dist; path_lengths[node_num].connected=connected;};
	Node * get_next_edge();	
	Edge<NODE_UNIT> * get_next_edge_obj();
	int get_node_num()                  {return(node_num);};
	void set_node_num(int num)          {node_num=num;};
	BOOL edge_exists(Node *link_node);
	Edge<NODE_UNIT> * get_edge(Node *link_node);
	int delete_edge(Node *old_edge_node);
	int delete_edge(Node *old_edge_node, int level);
	void add_edge(Node *new_edge_node);
	void add_edge(Node *new_edge_node, Edge<NODE_UNIT> *new_edge);
	void add_edge_internal(Node *new_edge_node, Edge<NODE_UNIT> *new_edge);
	void add_edge_allow_redund(Node *new_edge_node, Edge<NODE_UNIT> *new_edge);
	void add_edge(Node *new_edge_node, int weight);
	BOOL is_covered()               {return(node_covered);};
	BOOL is_undirected()			{return(undirected);};
	void set_covered()			    {node_covered=TRUE;};
	void set_uncovered()			    {node_covered=FALSE;};
	int get_component_num()         {return(component_num);};
	int get_score(Node<NODE_UNIT> *other) {if (edge_exists(other)==TRUE) return(1); else return(0);};
	double get_clustering_coeff();
	void calc_clustering_coeff();
	void set_component_num(int c)        {component_num=c;};
	NODE_UNIT* get_element()      {return(this_element);};
	void set_element(NODE_UNIT *ele)     {(*this_element)=(*ele);};
	void check_edges();

	
	void reset_out_edges()  {reset_edges();};
	void reset_in_edges()	{trans_in_edges=start_in_edges;};
	int get_num_out_edges() {return(get_num_edges());};
	int get_num_in_edges()  {return(num_in_edges);};
	Node * get_next_out_edge() {return(get_next_edge());};	
	Edge<NODE_UNIT> * get_next_out_edge_obj() {return(get_next_edge_obj());};
	Node * get_next_in_edge();	
	Edge<NODE_UNIT> * get_next_int_edge_obj();
	BOOL out_edge_exists(Node *link_node) {return(edge_exists(link_node));};
	BOOL in_edge_exists(Node *link_node);
	Edge<NODE_UNIT> * get_out_edge(Node *link_node) {return(get_edge(link_node));};
	Edge<NODE_UNIT> * get_in_edge(Node *link_node);
	int delete_out_edge(Node *old_edge_node) {return(delete_edge(old_edge_node));};
	int delete_in_edge(Node *in_edge_node);
	void add_out_edge(Node *new_edge_node) {add_edge(new_edge_node);};
	void add_in_edge(Node *new_edge_node);
	void add_in_edge(Node *new_edge_node, Edge<NODE_UNIT> *new_edge);
    int get_num_min_paths ()      {return(num_min_paths);};
    void set_num_min_paths(int num)   {num_min_paths=num;};
    void set_BFS_parent(Node * new_par)  {BFS_parent=new_par;};
    Node * get_BFS_parent()              {return(BFS_parent);};

	~Node();
	
protected:
	int num_edges, node_num, component_num, num_in_edges, num_min_paths;
	double clustering_coeff;
	BOOL node_covered, undirected, clustering_coeff_set, i_own_element;
	NODE_UNIT *this_element;
	Node_dist *path_lengths;
    Node *BFS_parent;
	//Linked_list<Node> *the_edges, *trans_edges, *start_edges;
	Linked_list< Edge<NODE_UNIT> > *the_edges, *trans_edges, *start_edges,
				*the_in_edges, *trans_in_edges, *start_in_edges;
};


template <class NODE_UNIT>
class Edge
{
	public:
		Edge();
		Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode);
		Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode, int weight);
		Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode, double weight);
		Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode, char *new_name);
        Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode,  string new_name);
		Edge & operator=(Edge &assign_from);

		int get_int_weight()					{return(int_weight);};
		double get_double_weight()					{return(double_weight);};
        char* get_c_name()						{return(name.c_str);};
        string get_name()                       {return(name);};
		void set_int_weight(int newweight)   {int_weight=newweight;};
		void set_double_weight(double newweight)		{double_weight=newweight;};
		Node<NODE_UNIT>* get_first_node()	{return(node1);};
		Node<NODE_UNIT>* get_second_node()	{return(node2);};
		Node<NODE_UNIT>* get_other_node(Node<NODE_UNIT> *me);

	protected:
		int int_weight;
		double double_weight;
		string name;
		Node<NODE_UNIT> *node1, *node2;

};


template <class NODE_UNIT>
class Graph {
	public:
		Linked_list< Linked_list< Node<NODE_UNIT> > > *fully_connected_comps, *start_comps;


		Graph();
		Graph(int n_nodes);
		Graph(int n_nodes, BOOL undir);
		Graph(int n_nodes, Node<NODE_UNIT> **the_nodes);
		Graph & operator=(Graph &assign_from);
		int get_num_nodes()             {return(num_nodes);};
		Node<NODE_UNIT>* get_node(int n);
		Node<NODE_UNIT>* find_node(NODE_UNIT* element);
		void delete_node (int node_num);
		void delete_node (Node<NODE_UNIT> *node);
		void add_edge(Node<NODE_UNIT> *node1, Node<NODE_UNIT> *node2);
		void add_edge(int node_num1, int node_num2);
		void delete_edge(Node<NODE_UNIT> *node1, Node<NODE_UNIT> *node2);
		void delete_edge(int node_num1, int node_num2);
		void reset_clust_coeff_calc()		{have_avg_clust_coeff=FALSE;};
		void reset_component_nums();
		void remove_self_interactions();
		int get_num_components()         {return(total_components);};
		void number_components();
        void number_components_iter();
        void process_component_node (Node<NODE_UNIT> *current, std::queue< Node<NODE_UNIT> *> *to_process, int component_num);
		double get_avg_clustering_coeff();
		int get_num_shared_edges(Node<NODE_UNIT> *node1, Node<NODE_UNIT> *node2);
		int get_num_shared_edges(int node1_num, int node2_num);
		void recurse_component(Node<NODE_UNIT> *current, int component);
		double get_min_path_lengths();
        void assign_num_min_paths();
		int get_total_edges();
		BOOL is_undirected();
		void find_fully_connected_components();
		int num_fully_connected_components()           {return(num_fully_connected_comps);};
		Linked_list< Node<NODE_UNIT> > * get_fully_connected_component_n(int n);
		int get_fully_connected_component_size(int n);
		void check_edges();
		void set_created_empty()   {created_empty=TRUE;};

		~Graph();

	protected:
		int num_nodes, total_components, total_edges, num_fully_connected_comps;
		double avg_clust_coeff;
		Node<NODE_UNIT> **nodes;
		BOOL counted_edges, created_empty, have_avg_clust_coeff, have_node_dists;
		Linked_list<int> *comp_sizes, *start_comp_sizes;
		
		void count_edges();
		void recurse_fully_connected_comps(Node<NODE_UNIT> *base_node, Linked_list< Node<NODE_UNIT> > *comp_so_far, int min_node_num);
		void prune_redundant_fully_connected_comps();
		BOOL is_list_subset(Linked_list <Node<NODE_UNIT> > *partial_list, Linked_list<Node<NODE_UNIT> > *complete_list);
		void delete_list(Linked_list <Node<NODE_UNIT> >*the_list);
};



template <class NODE_UNIT>
class Read_Graph {
	public:
	Read_Graph();
	Read_Graph(BOOL undir);

	virtual Graph<NODE_UNIT>* get_graph(char *nodefile, char *edgefile);
    virtual Graph<NODE_UNIT>* get_graph(string nodefile,  string edgefile);
	~Read_Graph();

	protected:
	int num_nodes, num_edges, *edge_weights;
    double *float_edge_weights;
    string node_filename, edge_filename, *edge_names;
    BOOL undirected, have_edge_names;
    ifstream nodein, edgein;
    NODE_UNIT *units, **node_edges;
    Node<NODE_UNIT> **the_nodes;
    Graph<NODE_UNIT> *the_graph;
    

	virtual int get_node_list()=0;
	virtual void get_edge_list()=0;
    virtual Node<NODE_UNIT>* find_node(NODE_UNIT *match);
		
};

template <class NODE_UNIT>
class EdgePair {
public:
    EdgePair() {empty=TRUE;}
    EdgePair<NODE_UNIT>& operator= (EdgePair<NODE_UNIT> &assign_from);
    
    BOOL empty;
    NODE_UNIT id1, id2;
};

template <class NODE_UNIT>
class Read_EdgeFile_Graph {
public:
    Read_EdgeFile_Graph();
    Read_EdgeFile_Graph(BOOL undir);
    Read_EdgeFile_Graph(BOOL undir, BOOL head);
    
    virtual Graph<NODE_UNIT>* get_graph(string edgefile);
    ~Read_EdgeFile_Graph();
    
protected:
    int num_nodes, num_edges;
    string edge_filename, *edgenames;
    BOOL undirected, have_edge_names, header;
    ifstream edgein;
    NODE_UNIT *units, **node_edges;
    Node<NODE_UNIT> **the_nodes;
    Graph<NODE_UNIT> *the_graph;
   
    
    virtual void read_edge_pair(EdgePair<NODE_UNIT>&)=0;
};


#endif
