#include <iostream>
#include "graph.h"
#include "linked_list.h"
#include <queue>


using namespace::std;

template <class LIST_UNIT>
Linked_list<LIST_UNIT>::Linked_list(Linked_list *new_last, LIST_UNIT *new_element)
{
	last=new_last;
	this_element=new_element;
	next=0;
}


template<class LIST_UNIT>
Linked_list<LIST_UNIT>::~Linked_list()
{
	//if (this_element != 0) {
	//	delete this_element;
	//	this_element=0;
	//}
		
}

template <class NODE_UNIT>
Node<NODE_UNIT>::Node()
{
	path_lengths=0;
	component_num =-1;
	num_edges=num_in_edges=0;
	undirected=TRUE;
	i_own_element=TRUE;
	
	the_edges=start_edges=trans_edges=the_in_edges=start_in_edges=trans_in_edges=0;
	node_covered=FALSE;
	clustering_coeff_set=FALSE;
    BFS_parent=0;

	//Call to default constructor of NODE_UNIT class--if bad idea, create element 
	//first and call Node<NODE_UNIT>(int num, NODE_UNIT *node_element)
	this_element=new NODE_UNIT;
}

template <class NODE_UNIT>
Node<NODE_UNIT>::Node(BOOL undir)
{
	path_lengths=0;
	component_num =-1;
	num_edges=num_in_edges=0;
	undirected=undir;
	i_own_element=TRUE;
    BFS_parent=0;
	
	the_edges=start_edges=trans_edges=the_in_edges=start_in_edges=trans_in_edges=0;
	node_covered=FALSE;
	clustering_coeff_set=FALSE;

	//Call to default constructor of NODE_UNIT class--if bad idea, create element 
	//first and call Node<NODE_UNIT>(int num, NODE_UNIT *node_element)
	this_element=new NODE_UNIT;
}



template <class NODE_UNIT>
Node<NODE_UNIT>::Node(int num, NODE_UNIT *node_element)
{
	path_lengths=0;
	node_num=num;
	this_element = node_element;
	component_num =-1;
	num_edges=num_in_edges=0;
	undirected=TRUE;
	i_own_element=FALSE;
    BFS_parent=0;
		
	the_edges=start_edges=trans_edges=the_in_edges=start_in_edges=trans_in_edges=0;
	node_covered=FALSE;
	clustering_coeff_set=FALSE;
}

template <class NODE_UNIT>
Node<NODE_UNIT>::Node(int num, NODE_UNIT *node_element, BOOL undir)
{
	path_lengths=0;
	node_num=num;
	this_element = node_element;
	component_num =-1;
	num_edges=num_in_edges=0;
	undirected=undir;
	i_own_element=FALSE;
    BFS_parent=0;
		
	the_edges=start_edges=trans_edges=the_in_edges=start_in_edges=trans_in_edges=0;
	node_covered=FALSE;
	clustering_coeff_set=FALSE;

}


template <class NODE_UNIT>
Node<NODE_UNIT> * Node<NODE_UNIT>::get_next_edge() 
{
	Node *retval;
	if (trans_edges != 0) {
		retval=trans_edges->get_element()->get_other_node(this);
		trans_edges=trans_edges->next;
	}
	else
		retval=0;

	
	return(retval);
}


template <class NODE_UNIT>
Node<NODE_UNIT> * Node<NODE_UNIT>::get_next_in_edge()
{
	Node *retval;
	
	if (trans_in_edges != 0) {
		retval=trans_in_edges->get_element()->get_other_node(this);
		trans_in_edges=trans_in_edges->next;
	}
	else
		retval=0;

	
	return(retval);

}


template <class NODE_UNIT>
Edge<NODE_UNIT> * Node<NODE_UNIT>::get_next_edge_obj()
{
	Edge<NODE_UNIT> *retval;

	if (trans_edges != 0) {
		retval=trans_edges->get_element();
		trans_edges=trans_edges->next;
	}
	else
		retval=0;

	
	return(retval);

}


template <class NODE_UNIT>
Edge<NODE_UNIT> * Node<NODE_UNIT>::get_next_int_edge_obj()
{
	Edge<NODE_UNIT> *retval;

	if (trans_in_edges != 0) {
		retval=trans_in_edges->get_element();
		trans_in_edges=trans_in_edges->next;
	}
	else
		retval=0;

	
	return(retval);


}


template <class NODE_UNIT>
BOOL Node<NODE_UNIT>::edge_exists(Node *link_node)
{
	BOOL retval=FALSE;
	Linked_list< Edge<NODE_UNIT> > *old_trans;

	old_trans=trans_edges;

	trans_edges=start_edges;

	while ((trans_edges != 0) && (retval == FALSE))
	{
		if (trans_edges->get_element()->get_other_node(this) == link_node)
			retval=TRUE;
		trans_edges=trans_edges->next;
	}
	trans_edges=old_trans;
	
	return(retval);
}



template <class NODE_UNIT>
BOOL Node<NODE_UNIT>::in_edge_exists(Node *link_node)
{
	BOOL retval=FALSE;
	Linked_list< Edge<NODE_UNIT> > *old_trans;

	old_trans=trans_in_edges;

	trans_in_edges=start_in_edges;

	while ((trans_in_edges != 0) && (retval == FALSE))
	{
		if (trans_in_edges->get_element()->get_other_node(this) == link_node)
			retval=TRUE;
		trans_in_edges=trans_in_edges->next;
	}
	trans_in_edges=old_trans;
	
	return(retval);



}


template <class NODE_UNIT>
Edge<NODE_UNIT> * Node<NODE_UNIT>::get_edge(Node *link_node)
{
	Edge<NODE_UNIT> *retval=0;
	Linked_list< Edge<NODE_UNIT> > *old_trans;

	old_trans=trans_edges;

	trans_edges=start_edges;

	while ((trans_edges != 0) && (retval == 0))
	{
		if (trans_edges->get_element()->get_other_node(this) == link_node)
			retval=trans_edges->get_element();
		trans_edges=trans_edges->next;
	}
	trans_edges=old_trans;
	
	return(retval);
}


template <class NODE_UNIT>
Edge<NODE_UNIT> * Node<NODE_UNIT>::get_in_edge(Node *link_node)
{
	Edge<NODE_UNIT> *retval=0;
	Linked_list< Edge<NODE_UNIT> > *old_trans;

	old_trans=trans_in_edges;

	trans_in_edges=start_in_edges;

	while ((trans_in_edges != 0) && (retval == 0))
	{
		if (trans_in_edges->get_element()->get_other_node(this) == link_node)
			retval=trans_in_edges->get_element();
		trans_in_edges=trans_in_edges->next;
	}
	trans_in_edges=old_trans;
	
	return(retval);

}



template <class NODE_UNIT>
void Node<NODE_UNIT>::add_edge(Node *new_edge_node)
{
	Edge<NODE_UNIT> *new_edge;

	if (edge_exists(new_edge_node) == FALSE) {
		new_edge=new Edge<NODE_UNIT>(this, new_edge_node);

		if (the_edges == 0) {
			the_edges=new Linked_list< Edge<NODE_UNIT> >(0, new_edge);
			start_edges=the_edges;

			if (undirected == TRUE)
				start_in_edges=the_in_edges=the_edges;
		}
		else {
			the_edges->next=new Linked_list< Edge<NODE_UNIT> >(the_edges, new_edge);
			the_edges=the_edges->next;

			if (undirected == TRUE) {
				the_in_edges->next=the_edges;
				the_in_edges=the_in_edges->next;
			}
		}
		num_edges++;
		if (undirected == TRUE) {
			new_edge_node->add_edge_internal(this, new_edge);
			num_in_edges++;
		}
		else
			new_edge_node->add_in_edge(this);
	}
}


template <class NODE_UNIT>
void Node<NODE_UNIT>::add_edge(Node *new_edge_node, Edge<NODE_UNIT> *new_edge)
{
	Edge<NODE_UNIT> *first_edge;
	
	if (edge_exists(new_edge_node) == FALSE) {
		first_edge=new Edge<NODE_UNIT>(this, new_edge_node);
		(*first_edge)=(*new_edge);
		
		if (the_edges == 0) {
			the_edges=new Linked_list< Edge<NODE_UNIT> >(0, first_edge);
			start_edges=the_edges;

			if (undirected == TRUE) 
				the_in_edges=start_in_edges=the_edges;
		}
		else {
			the_edges->next=new Linked_list< Edge<NODE_UNIT> >(the_edges, first_edge);
			the_edges=the_edges->next;

			if (undirected == TRUE) {
				the_in_edges->next=the_edges;
				the_in_edges=the_in_edges->next;
			}
			
		}
		num_edges++;
		if (undirected == TRUE) {
			if (! new_edge_node->edge_exists(this)) {
				new_edge_node->add_edge_internal(this, first_edge);
			}
			num_in_edges++;
		}
		else {
			if (! new_edge_node->in_edge_exists(this)) {
					new_edge_node->add_in_edge(this, first_edge);
			}
		}
	}
}


template <class NODE_UNIT>
void Node<NODE_UNIT>::add_edge_internal(Node *new_edge_node, Edge<NODE_UNIT> *first_edge)
{
	if (edge_exists(new_edge_node) == FALSE) {
		if (the_edges == 0) {
			the_edges=new Linked_list< Edge<NODE_UNIT> >(0, first_edge);
			start_edges=the_edges;
			
			if (undirected == TRUE) 
				the_in_edges=start_in_edges=the_edges;
		}
		else {
			the_edges->next=new Linked_list< Edge<NODE_UNIT> >(the_edges, first_edge);
			the_edges=the_edges->next;
			
			if (undirected == TRUE) {
				the_in_edges->next=the_edges;
				the_in_edges=the_in_edges->next;
			}
			
		}
		num_edges++;
		if (undirected == TRUE) {
			num_in_edges++;
		}
	}
}




template <class NODE_UNIT>
void Node<NODE_UNIT>::add_edge_allow_redund(Node *new_edge_node, Edge<NODE_UNIT> *new_edge)
{
	Edge<NODE_UNIT> *first_edge;
	
	first_edge=new Edge<NODE_UNIT>(this, new_edge_node);
	(*first_edge)=(*new_edge);
	
		if (the_edges == 0) {
		the_edges=new Linked_list< Edge<NODE_UNIT> >(0, first_edge);
		start_edges=the_edges;
		
		if (undirected == TRUE) 
			the_in_edges=start_in_edges=the_edges;
	}
	else {
		the_edges->next=new Linked_list< Edge<NODE_UNIT> >(the_edges, first_edge);
		the_edges=the_edges->next;
		
		if (undirected == TRUE) {
			the_in_edges->next=the_edges;
			the_in_edges=the_in_edges->next;
		}
		
	}
	num_edges++;
	if (undirected == TRUE) {
		new_edge_node->add_edge(this);
		num_in_edges++;
	}
	else
		new_edge_node->add_in_edge(this);

}



template <class NODE_UNIT>
void Node<NODE_UNIT>::add_edge(Node *new_edge_node, int weight)
{
	Edge<NODE_UNIT> *new_edge;

	if (edge_exists(new_edge_node) == FALSE) {
		new_edge=new Edge<NODE_UNIT>(this, new_edge_node, weight);

		if (the_edges == 0) {
			the_edges=new Linked_list< Edge<NODE_UNIT> >(0, new_edge);
			start_edges=the_edges;

			if (undirected == TRUE)
				the_in_edges=start_in_edges=the_edges;
		}
		else {
			the_edges->next=new Linked_list< Edge<NODE_UNIT> >(the_edges, new_edge);
			the_edges=the_edges->next;

			if (undirected == TRUE) {
				the_in_edges->next=the_edges;
				the_in_edges=the_in_edges->next;
			}
		}
		num_edges++;
		if (undirected == TRUE) {
			new_edge_node->add_edge_internal(this, new_edge);
			num_in_edges++;
		}
		else
			//WARNING--BAD--WEIGHT IGNORED!!!
			new_edge_node->add_in_edge(this);
	}
}


template <class NODE_UNIT>
void Node<NODE_UNIT>::add_in_edge(Node *new_edge_node)
{
	Edge<NODE_UNIT> *new_edge;

	if (in_edge_exists(new_edge_node) == FALSE) {
		new_edge=new Edge<NODE_UNIT>(this, new_edge_node);

		if (the_in_edges == 0) {
			the_in_edges=new Linked_list< Edge<NODE_UNIT> >(0, new_edge);
			start_in_edges=the_in_edges;

		}
		else {
			the_in_edges->next=new Linked_list< Edge<NODE_UNIT> >(the_in_edges, new_edge);
			the_in_edges=the_in_edges->next;	
		}
		num_in_edges++;
	}
}


template <class NODE_UNIT>
void Node<NODE_UNIT>::add_in_edge(Node *new_edge_node, Edge<NODE_UNIT> *new_edge)
{
	Edge<NODE_UNIT> *first_edge;
	
	if (in_edge_exists(new_edge_node) == FALSE) {
		first_edge=new Edge<NODE_UNIT>(this, new_edge_node);
		(*first_edge)=(*new_edge);
		if (the_in_edges == 0) {
			the_in_edges=new Linked_list< Edge<NODE_UNIT> >(0, first_edge);
			start_in_edges=the_in_edges;
		}
		else {
			the_in_edges->next=new Linked_list< Edge<NODE_UNIT> >(the_in_edges, first_edge);
			the_in_edges=the_in_edges->next;	
		}
		num_in_edges++;
	}
}



template <class NODE_UNIT>
int Node<NODE_UNIT>::delete_edge(Node *old_edge_node)
{
	return(delete_edge(old_edge_node, 0));
}

template <class NODE_UNIT>
int Node<NODE_UNIT>::delete_edge(Node *old_edge_node, int level)
{
	int ret_val=0;
	Edge<NODE_UNIT> *the_edge;
	Linked_list< Edge<NODE_UNIT> > *trans_list, *save_last, *save_next;
	
	if (start_edges != 0) {
		trans_list=start_edges;
		while((trans_list->next != 0) &&
				( trans_list->get_element()->get_other_node(this) != old_edge_node )) 
					trans_list=trans_list->next;
		if (trans_list->get_element()->get_other_node(this) == old_edge_node )
		{
			the_edge=trans_list->get_element();

			if (trans_list == start_edges) {
				start_edges=trans_list->next;
				if (undirected == TRUE)
					start_in_edges=trans_list->next;
			}
			if (trans_list == the_edges) 
					the_edges = trans_list->last;
			num_edges--;
			if (undirected == TRUE) {
				num_in_edges--;
				if (trans_list == the_in_edges)
					the_in_edges = trans_list->last;
			}
			ret_val=1;
			save_last=trans_list->last;
			save_next=trans_list->next;
			
			//Only delete the edge object once both sides of the undirected relationship are 
			//correct
			if ((undirected == FALSE) || (level != 0))
				delete the_edge;
			delete trans_list;
			if (save_last != 0)
				save_last->next=save_next;
			if(save_next != 0)
				save_next->last=save_last;
		}

		//In an undirected graph, delete the reciprocal edge
		if ((undirected == TRUE) && (level==0))
				old_edge_node->delete_edge(this, level+1);
		else if (undirected == FALSE)
			old_edge_node->delete_in_edge(this);
	}

	return(ret_val);
}



template <class NODE_UNIT>
int Node<NODE_UNIT>::delete_in_edge(Node *old_edge_node)
{
	int ret_val=0;
	Edge<NODE_UNIT> *the_edge;
	Linked_list< Edge<NODE_UNIT> > *trans_list, *save_last, *save_next;
	
	if (start_in_edges != 0) {
		trans_list=start_in_edges;
		while((trans_list->next != 0) &&
				( trans_list->get_element()->get_other_node(this) != old_edge_node )) 
					trans_list=trans_list->next;
		if (trans_list->get_element()->get_other_node(this) == old_edge_node )
		{
			the_edge=trans_list->get_element();

			if (trans_list == start_in_edges)
				start_in_edges=trans_list->next;
			if (trans_list == the_in_edges)
					the_in_edges = trans_list->last;
			num_in_edges--;
			ret_val=1;
			save_last=trans_list->last;
			save_next=trans_list->next;
			
		
			delete the_edge;
			delete trans_list;
			if (save_last != 0)
				save_last->next=save_next;
			if(save_next != 0)
				save_next->last=save_last;
		}
	
	}

	return(ret_val);
}


template <class NODE_UNIT>
void Node<NODE_UNIT>::check_edges() 
{
	Linked_list< Edge<NODE_UNIT> > *trans_list;

	if (num_edges > 1) {
		trans_list=start_edges;

		while(trans_list != 0) {
			if ((trans_list->next ==0) && (trans_list->last == 0))
				cerr<<"Invalid edge list\n";
			trans_list=trans_list->next;
		}
	}

		if (num_in_edges > 1) {
		trans_list=start_in_edges;

		while(trans_list != 0) {
			if ((trans_list->next ==0) && (trans_list->last == 0))
				cerr<<"Invalid edge list\n";
            
            if ((trans_list->get_element()->get_first_node() == 0) || (trans_list->get_element()->get_first_node() == 0)) {cerr<<"ERROR edge from node "<<node_num<<" has missing target\n";}
			trans_list=trans_list->next;
		}
	}
}


template <class NODE_UNIT>
void Node<NODE_UNIT>::initialize_path_lens(int num_nodes)
{
	int i;
	Node<NODE_UNIT> *curr_edge;
	
	if (path_lengths == 0)
		path_lengths = new Node_dist [num_nodes];

	for(i=0; i<num_nodes; i++) {
		path_lengths[i].dist=0;
		path_lengths[i].connected=FALSE;
	}
	
	reset_edges();

	path_lengths[node_num].dist=0;
	path_lengths[node_num].connected=TRUE;

	for(i=0; i<num_edges; i++) {
		curr_edge=get_next_edge();
		if (curr_edge->get_node_num() != node_num) {
			path_lengths[curr_edge->get_node_num()].dist=1;
			path_lengths[curr_edge->get_node_num()].connected=TRUE;
		}
	}
}


template <class NODE_UNIT>
Node_dist * Node<NODE_UNIT>::get_dist_to_node(int node_num)
{
	if (path_lengths!=0)
		return(&path_lengths[node_num]);
	else
		return(0);
}

template <class NODE_UNIT>
double Node<NODE_UNIT>::get_clustering_coeff()
{
	if (clustering_coeff_set == FALSE)
		calc_clustering_coeff();
	return(clustering_coeff);
}


template <class NODE_UNIT>
void Node<NODE_UNIT>::calc_clustering_coeff()
{
	int i,j, count_edges=0;
	Node<NODE_UNIT> *first_node, *next_node;

	reset_edges();

	for(i=0; i<num_edges; i++) {
		first_node=get_next_edge();

		for(j=i+1; j<num_edges; j++) {
			next_node=get_next_edge();
			if (first_node->edge_exists(next_node) == TRUE)
				count_edges++;

		}
		reset_edges();
		for(j=0; j<=i; j++)
			next_node=get_next_edge();

	}

	//cout<<node_num<<": count: "<<count_edges<<endl;
	if (num_edges > 1)
		clustering_coeff=(double)(2*count_edges)/(double)(num_edges*(num_edges-1));
	else
		clustering_coeff=1.0;
}



template <class NODE_UNIT>
Node<NODE_UNIT>::~Node()
{
	Linked_list< Edge<NODE_UNIT> > *temp;

	if (i_own_element == TRUE) delete this_element;
	
	trans_edges = start_edges;
	
	while (trans_edges != 0) {
		temp=trans_edges->next;

		delete trans_edges;
		trans_edges=temp;
	}

	if (undirected == FALSE) {
		trans_in_edges = start_in_edges;
	
		while (trans_in_edges != 0) {
			temp=trans_in_edges->next;
			delete trans_in_edges;
			trans_in_edges=temp;

		}
	}

	if(path_lengths != 0)
		delete[] path_lengths;
}

template <class NODE_UNIT>
Edge<NODE_UNIT>::Edge()
{
	cerr<<"Error: call to default constructor of class Edge\n";
	node1=node2=0;
}

template <class NODE_UNIT>
Edge<NODE_UNIT>::Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode)
{
	node1=fnode;
	node2=snode;
	int_weight=1;
	double_weight=1.0;
}
		
template <class NODE_UNIT>
Edge<NODE_UNIT>::Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode, int weight)
{
	node1=fnode;
	node2=snode;
	int_weight=weight;
}
		
template <class NODE_UNIT>
Edge<NODE_UNIT>::Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode, double weight)
{
	node1=fnode;
	node2=snode;
	double_weight=weight;
}


template <class NODE_UNIT>
Edge<NODE_UNIT>::Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode, char * new_name)
{
    string *temp;
    
	node1=fnode;
	node2=snode;
	temp= new string (new_name);
    name=temp;
    delete temp;
}


template <class NODE_UNIT>
Edge<NODE_UNIT>::Edge(Node<NODE_UNIT> *fnode, Node<NODE_UNIT> *snode, string new_name)
{
	node1=fnode;
	node2=snode;
    name=new_name;
}

template <class NODE_UNIT>
Edge<NODE_UNIT> & Edge<NODE_UNIT>::operator=(Edge &assign_from)
{
	BOOL save_nodes=FALSE;
	if (node1 == assign_from.get_first_node()) {
			if (node2 == assign_from.get_second_node())
				save_nodes=TRUE;
	}
	else if (node1 == assign_from.get_second_node()) {
			if (node2 == assign_from.get_first_node())
				save_nodes=TRUE;
	}
	
	if (save_nodes == FALSE) {
		node1=assign_from.get_first_node();
		node2=assign_from.get_second_node();
	}
	
	name = assign_from.get_name();
	double_weight=assign_from.get_double_weight();	
	int_weight=assign_from.get_int_weight();
    return(*this);
}


template <class NODE_UNIT>
Node<NODE_UNIT> * Edge<NODE_UNIT>::get_other_node(Node<NODE_UNIT> *me)
{
	if (node1==me) 
		return(node2);
	else
		return(node1);
}

template <class NODE_UNIT>
Graph<NODE_UNIT>::Graph()                         
{
	total_edges=0;
	counted_edges=FALSE;
	have_avg_clust_coeff=FALSE;
	fully_connected_comps=start_comps=0;
	comp_sizes=start_comp_sizes=0;
	nodes=0;
	num_nodes=0;
    have_node_dists=FALSE;

}


//Creates an "empty" graph where the nodes have no data in them
template <class NODE_UNIT>
Graph<NODE_UNIT>::Graph(int n_nodes)
{
	int i;
	
	created_empty=TRUE;
	have_avg_clust_coeff=FALSE;
	counted_edges=FALSE;
	num_nodes=n_nodes;
	nodes = new Node<NODE_UNIT>* [num_nodes];		
	fully_connected_comps=start_comps=0;
	comp_sizes=start_comp_sizes=0;
    have_node_dists=FALSE;

	for (i=0; i<num_nodes; i++) {
		nodes[i]=new Node<NODE_UNIT>;
		nodes[i]->set_node_num(i);	
	}

}

//Creates an "empty" graph where the nodes have no data in them with directed/undirected value
//taken from undir
template <class NODE_UNIT>
Graph<NODE_UNIT>::Graph(int n_nodes, BOOL undir)
{
	int i;

	created_empty=TRUE;
	have_avg_clust_coeff=FALSE;
	counted_edges=FALSE;
	num_nodes=n_nodes;
	have_node_dists=FALSE;
    
	nodes = new Node<NODE_UNIT>* [num_nodes];	

	fully_connected_comps=start_comps=0;
	comp_sizes=start_comp_sizes=0;

	for (i=0; i<num_nodes; i++) {
	
		nodes[i]=new Node<NODE_UNIT>(undir);
	

		nodes[i]->set_node_num(i);	
	}

}



template <class NODE_UNIT>
Graph<NODE_UNIT>::Graph(int n_nodes, Node<NODE_UNIT> **the_nodes)
{
	int i;
	
	created_empty=FALSE;
	have_avg_clust_coeff=FALSE;
	counted_edges=FALSE;
	num_nodes=n_nodes;
	nodes = new Node<NODE_UNIT>* [num_nodes];		
	fully_connected_comps=start_comps=0;
	comp_sizes=start_comp_sizes=0;
    have_node_dists=FALSE;
    
	for (i=0; i<num_nodes; i++) {
		nodes[i]=the_nodes[i];
		the_nodes[i]->set_node_num(i);	
	}

}



template <class NODE_UNIT>
Graph<NODE_UNIT> & Graph<NODE_UNIT>::operator=(Graph<NODE_UNIT> &assign_from)
{
	int i, j, k;
	NODE_UNIT new_element;
	Node<NODE_UNIT> *edge_dest, *curr_edge;

	if (num_nodes != assign_from.get_num_nodes()) {
		cerr<<"Error: assignment from graph of unequal size\n";
	}
	else if (created_empty == FALSE)
	{
		cerr<<"Error: Can't perform assignment on graph with non-empty members\n";
	}
	else
	{
		for (i=0; i<num_nodes; i++)
		{
			nodes[i]->reset_edges();
			curr_edge=nodes[i]->get_next_edge();
			while(curr_edge !=0)
			{
				nodes[i]->delete_edge(curr_edge);
				nodes[i]->reset_edges();
				curr_edge=nodes[i]->get_next_edge();
			}
		}

		total_components=assign_from.get_num_components();
		for (i=0; i<num_nodes; i++)
		{
			new_element = (*assign_from.get_node(i)->get_element());
			nodes[i]->set_element(&new_element);
			nodes[i]->set_node_num(assign_from.get_node(i)->get_node_num());
            nodes[i]->set_component_num(assign_from.get_node(i)->get_component_num());
		}
		for (i=0; i<num_nodes; i++) {
			assign_from.get_node(i)->reset_edges();
			
			for(j=0; j<assign_from.get_node(i)->get_num_edges(); j++)
			{
				edge_dest=assign_from.get_node(i)->get_next_edge();
				k=0;
				while(nodes[k]->get_node_num() != edge_dest->get_node_num()) {k++;}

				nodes[i]->add_edge(nodes[k]);
				nodes[i]->get_edge(nodes[k])->set_int_weight(assign_from.get_node(i)->get_edge(edge_dest)->get_int_weight());
				nodes[i]->get_edge(nodes[k])->set_double_weight(assign_from.get_node(i)->get_edge(edge_dest)->get_double_weight());

			}


		}

	}
	return(*this);
}


template <class NODE_UNIT>
Node<NODE_UNIT>* Graph<NODE_UNIT>::get_node(int n)
{
	if (n<num_nodes)
		return(nodes[n]);
	else
		return(0);
}


template <class NODE_UNIT>
Node<NODE_UNIT>* Graph<NODE_UNIT>::find_node(NODE_UNIT *element) 
{
	int i=0;

	while ((i<num_nodes) && (!((*element) == (* nodes[i]->get_element())))) i++;
	
	if (i<num_nodes) return(nodes[i]);
	else return(0);
}

template <class NODE_UNIT>
void Graph<NODE_UNIT>::reset_component_nums()
{
	int i;

	for(i=0; i<num_nodes; i++) {
		nodes[i]->set_uncovered();
		nodes[i]->set_component_num(-1);
	}
}

template <class NODE_UNIT>
void Graph<NODE_UNIT>::number_components()
{
	int i;

	total_components=0;
	for(i=0; i<num_nodes; i++) 
		if (nodes[i]->get_component_num() == -1)
			recurse_component(nodes[i], total_components++);

	
}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::number_components_iter()
{
    int i, j;
    BOOL all_unvisited;
    Node<NODE_UNIT> *next_node;
    std::queue< Node<NODE_UNIT> *> to_process;
    
    
    total_components=0;
    for(i=0; i<num_nodes; i++) {
        if (nodes[i]->get_component_num() == -1) {
            process_component_node(nodes[i], &to_process, total_components);
            
            while (!to_process.empty()) {
                next_node=to_process.front();
                to_process.pop();
                cout<<"Processing: "<<next_node->get_node_num()<<endl;
                
                process_component_node(next_node, &to_process, total_components);
            }
            
            total_components++;
        }
    }
    
}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::process_component_node (Node<NODE_UNIT> *current, std::queue< Node<NODE_UNIT> *> *to_process, int component_num)
{
    Node<NODE_UNIT> *next_node;
    
    if (current->is_covered() == FALSE) {
        current->set_covered();
        current->set_component_num(component_num);
        current->reset_edges();
        
        next_node=current->get_next_edge();
        while (next_node != 0) {
            if (next_node->is_covered() == FALSE)
                to_process->push(next_node);
            next_node=current->get_next_edge();
        }
    }
}

template <class NODE_UNIT>
void Graph<NODE_UNIT>::delete_node (int node_num)
{
	delete_node(nodes[node_num]);

}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::delete_node (Node<NODE_UNIT> *node)
{
	int i, node_cnt=0;
	Node<NODE_UNIT> **new_node_array;

	new_node_array=new Node<NODE_UNIT>* [num_nodes-1];

	node->reset_edges();
	for(i=0; i<node->get_num_edges(); i++) 
		node->delete_edge(node->get_next_edge());
	

	for(i=0; i<num_nodes; i++) {
		if (nodes[i] != node) {
			new_node_array[node_cnt]=nodes[i];
			new_node_array[node_cnt]->set_node_num(node_cnt);
			node_cnt++;
		}

	}
	delete[] nodes;
	nodes=new_node_array;
	delete node;
	num_nodes--;
	counted_edges=FALSE;
}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::add_edge(Node<NODE_UNIT> *node1, Node<NODE_UNIT> *node2)
{
	counted_edges=FALSE;
	node1->add_edge(node2);
}
	

template <class NODE_UNIT>
void Graph<NODE_UNIT>::add_edge(int node_num1, int node_num2)
{
	add_edge(nodes[node_num1], nodes[node_num2]);
}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::delete_edge(Node<NODE_UNIT> *node1, Node<NODE_UNIT> *node2)
{
	if (node1->edge_exists(node2) == TRUE)  {
		counted_edges=FALSE;
		node1->delete_edge(node2);
	}
}
		
template <class NODE_UNIT>
void Graph<NODE_UNIT>::delete_edge(int node_num1, int node_num2)
{
	delete_edge(nodes[node_num1], nodes[node_num2]);
}



template <class NODE_UNIT>
void Graph<NODE_UNIT>::find_fully_connected_components()
{
	int i, j;
	Linked_list< Node<NODE_UNIT> > *new_list;

	num_fully_connected_comps=0;

	for (i=0; i<num_nodes; i++) {
		//Stupid--but makes sure the nodes are search in their numbered order
		//which may not correspond to their order in the nodes array
		j=0;
		while(nodes[j]->get_node_num() != i)
			j++;


		new_list=new Linked_list< Node<NODE_UNIT> > (0, nodes[j]);
		recurse_fully_connected_comps(nodes[j], new_list, nodes[j]->get_node_num());
		delete new_list;
	}

	//Some of our connected components may be subsets of each other--fix this
	prune_redundant_fully_connected_comps();
}

		

template <class NODE_UNIT>
Linked_list< Node<NODE_UNIT> > * Graph<NODE_UNIT>::get_fully_connected_component_n(int n)
{
	int i=0;

	if (n<num_fully_connected_comps) {
		fully_connected_comps=start_comps;
		while(i<n) {
			fully_connected_comps=fully_connected_comps->next;
			i++;
		}
		return(fully_connected_comps->get_element());
	}
	else {
		cerr<<"ERROR: request for fully-connected component "<<n<<" when there are only "
			<<num_fully_connected_comps<<" total fully connected components\n";
		return(0);
	}

}

template <class NODE_UNIT>
int Graph<NODE_UNIT>::get_fully_connected_component_size(int n)
{
		int i=0;

	if (n<num_fully_connected_comps) {
		comp_sizes=start_comp_sizes;
		while(i<n) {
			comp_sizes=comp_sizes->next;
			i++;
		}
		return((*comp_sizes->get_element()));
	}
	else {
		cerr<<"ERROR: request for size of fully-connected component "<<n<<" when there are only "
			<<num_fully_connected_comps<<" total fully connected components\n";
		return(0);
	}

}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::recurse_component(Node<NODE_UNIT> *current, int component)
{
	Node<NODE_UNIT> *next_node;


	
	current->set_covered();
	current->set_component_num(component);
	if (current->get_num_edges() != 0) {
		current->reset_edges();
		next_node=current->get_next_edge();

	
		while (next_node != 0) {
			if (next_node->is_covered() == FALSE) 
				recurse_component(next_node, component);
				
			next_node=current->get_next_edge();
		}
	}
}


template <class NODE_UNIT>
double Graph<NODE_UNIT>::get_min_path_lengths()
{
	int i, j ,k, l, cnt;
	double retval;

	for(i=0; i<num_nodes; i++)
		nodes[i]->initialize_path_lens(num_nodes);

	for(i=0; i<num_nodes; i++) {
		if (nodes[i]->get_num_edges() != 0) {
			for(j=0; j<num_nodes; j++) {
				for(k=0; k<num_nodes; k++) {
					if ((nodes[j]->get_dist_to_node(i)->connected == TRUE) && (nodes[i]->get_dist_to_node(k)->connected == TRUE)) {
						if ((nodes[j]->get_dist_to_node(k)->connected == FALSE) || 
							((nodes[j]->get_dist_to_node(i)->dist + nodes[i]->get_dist_to_node(k)->dist) 
							< nodes[j]->get_dist_to_node(k)->dist))
							nodes[j]->set_dist_to_node(k, nodes[j]->get_dist_to_node(i)->dist + 
							nodes[i]->get_dist_to_node(k)->dist, TRUE);
					}

				}
			}
		}

	}

	retval=0.0;
	cnt=0;

	for(i=0; i<num_nodes; i++) {
	//	cout<<i<<"\t";
		for(j=0; j<num_nodes; j++) {
			if (nodes[i]->get_dist_to_node(j)->connected == TRUE) {
	//			cout<<nodes[i]->get_dist_to_node(j)->dist<<"\t";
				retval+=nodes[i]->get_dist_to_node(j)->dist;
				cnt++;
			}
		//	else
	///			cout<<"INF\t";
		}
	//	cout<<endl;
	}
    have_node_dists=TRUE;

	retval=retval/(1.0*cnt);
	return(retval);

}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::assign_num_min_paths()
{
    int i,j,k;
    if (have_node_dists==FALSE) get_min_path_lengths();
    
    for(i=0; i<num_nodes; i++) {
        nodes[i]->set_num_min_paths(0);
        
        for (j=0; j<num_nodes; j++) {
            if (i!=j) {
                if (nodes[i]->get_dist_to_node(j)->connected == TRUE) {
                    for(k=j+1; k<num_nodes; k++) {
                        if ((k!=i) &&(nodes[k]->get_dist_to_node(j)->connected==TRUE)) {
                            if ((nodes[i]->get_dist_to_node(j)->dist + nodes[i]->get_dist_to_node(k)->dist) == nodes[j]->get_dist_to_node(k)->dist)
                                nodes[i]->set_num_min_paths(nodes[i]->get_num_min_paths()+1);
                        }
                    }
                }
            }
        }
    }
    
}


template <class NODE_UNIT>
int Graph<NODE_UNIT>::get_total_edges()
{
	if (counted_edges == FALSE)
		count_edges();
	return(total_edges);
}


template <class NODE_UNIT>
int Graph<NODE_UNIT>::get_num_shared_edges(Node<NODE_UNIT> *node1, Node<NODE_UNIT> *node2)
{
	int i, num_shared=0;
	Node<NODE_UNIT> *edge_connect_node;

	
	node1->reset_edges();
	for(i=0; i<node1->get_num_edges(); i++) {
		edge_connect_node=node1->get_next_edge();

		if (node2->edge_exists(edge_connect_node)==TRUE) 
			num_shared++;
	}

	return(num_shared);
}
		
template <class NODE_UNIT>
int Graph<NODE_UNIT>::get_num_shared_edges(int node1_num, int node2_num)
{
	return(get_num_shared_edges(nodes[node1_num], nodes[node2_num]));
}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::check_edges() 
{
	int i;

	for (i=0; i<num_nodes; i++)
		nodes[i]->check_edges();
}



template <class NODE_UNIT>
void Graph<NODE_UNIT>::remove_self_interactions()
{
	int i;
	Node<NODE_UNIT> *delete_edge;

	for(i=0; i<num_nodes; i++) {
		nodes[i]->reset_edges();
		delete_edge=nodes[i]->get_next_edge();
		while(delete_edge != 0) {
			if (nodes[i]->get_node_num() == delete_edge->get_node_num()) {
				nodes[i]->delete_edge(delete_edge);
				nodes[i]->reset_edges();
			}
			delete_edge=nodes[i]->get_next_edge();
		}

	}	
}


template <class NODE_UNIT>
Graph<NODE_UNIT>::~Graph()
{
	int i, j;
	Node<NODE_UNIT> *curr_node, *other_node;
	Linked_list<int> *tempsize;
	Linked_list<Linked_list< Node<NODE_UNIT> > > *temp;

	for(i=0; i<num_nodes; i++) {
		curr_node=nodes[i];
		curr_node->reset_edges();
		other_node=curr_node->get_next_edge();
		while (other_node != 0) {
			curr_node->delete_edge(other_node);
			curr_node->reset_edges();
			other_node=curr_node->get_next_edge();
		}
	}
	
	
	
	if (nodes != 0) {
		if (created_empty == TRUE) {
			for(i=0; i<num_nodes; i++)	delete nodes[i];
		}
		delete[] nodes;

	}
	fully_connected_comps=start_comps;
	comp_sizes=start_comp_sizes;
	while (fully_connected_comps != 0) {		
		temp=fully_connected_comps->next;
		delete_list(fully_connected_comps->get_element());
		delete fully_connected_comps;
		fully_connected_comps=temp;

		tempsize=comp_sizes->next;
		delete comp_sizes;
		comp_sizes=tempsize;
	}
}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::count_edges()
{
	int i;
	Node<NODE_UNIT> *curr_node, *other_node;
	
	total_edges=0;
	counted_edges=TRUE;

	if (nodes[0]->is_undirected() == FALSE) {
	for(i=0; i<num_nodes; i++) 
			total_edges+=nodes[i]->get_num_edges();
	}
	else {
		for(i=0; i<num_nodes; i++)  {
			curr_node=nodes[i];
			curr_node->reset_edges();
			other_node=curr_node->get_next_edge();
			while(other_node != 0) {
				if (curr_node->get_node_num() <= other_node->get_node_num()) total_edges++;
				other_node=curr_node->get_next_edge();
			}
		}
	}

	//if (nodes[0]->is_undirected() == TRUE)
	//	total_edges/=2;
}


template <class NODE_UNIT>
double Graph<NODE_UNIT>::get_avg_clustering_coeff()
{
	int i;

	avg_clust_coeff=0.0;


	if (have_avg_clust_coeff==FALSE) {
		for(i=0; i<num_nodes; i++) 
			avg_clust_coeff+=nodes[i]->get_clustering_coeff();
		avg_clust_coeff/=num_nodes;
		have_avg_clust_coeff=TRUE;
	}
	return(avg_clust_coeff);

}

template <class NODE_UNIT> 
BOOL Graph<NODE_UNIT>::is_undirected()
{
	return(nodes[0]->is_undirected());
}

template <class NODE_UNIT>
void Graph<NODE_UNIT>::recurse_fully_connected_comps(Node<NODE_UNIT> *base_node, 
													Linked_list< Node<NODE_UNIT> > *comp_so_far, int min_node_num)
{
	int i, j, *size;
	BOOL connect_to_all, at_least_one_added;
	Node<NODE_UNIT> *next_node;
	Linked_list<Node<NODE_UNIT> > *start_list, *list, *new_list, *start_new_list;
	
	list=start_list=comp_so_far;

	//Will indicate if we've added anything to the component we received
	at_least_one_added=FALSE;

	for (i=0; i<base_node->get_num_edges(); i++) {
		//We allow each of our "base_nodes" edges to be part of a possible component
		//At each step we re-walk through the edge list because our position in 
		//the edge list may change due to recursion
		base_node->reset_edges();
		for(j=0; j<=i; j++)
			next_node=base_node->get_next_edge();
		
		if (next_node->get_node_num() > min_node_num) {
			//This node could be added to the fully connected component
			connect_to_all=TRUE;
			
			//Check to see if this node is connected to all nodes already in the component
			while((list != 0) && (connect_to_all == TRUE)) {
				if (next_node->edge_exists(list->get_element()) == FALSE)
					connect_to_all=FALSE;
				list=list->next;
			}

			if (connect_to_all==TRUE) {
			//Recurse on this new partial connected component

				//This component is not complete
				at_least_one_added=TRUE;

				//First we copy the current list
				list=start_list;
				new_list=new Linked_list< Node<NODE_UNIT> > (0, list->get_element());
				start_new_list=new_list;
				list=list->next;

				while(list != 0) {
					new_list->next=new Linked_list< Node<NODE_UNIT> > (new_list, list->get_element());
					new_list=new_list->next;
					list=list->next;
				}

				//Now add the new node
				new_list->next=new Linked_list< Node<NODE_UNIT> > (new_list, next_node);
				new_list=start_new_list;

				//Recurse
				recurse_fully_connected_comps(base_node, new_list, next_node->get_node_num());

				//Clean up
				new_list=start_new_list;
				while(new_list != 0) {
					start_new_list=new_list->next;
					delete new_list;
					new_list=start_new_list;
				}
			}
		}

	}

	if (at_least_one_added == FALSE) {
		//We received a complete component--store it
		//First copy it to make deallocation of intermediate lists easy
		size=new int;
		(*size)=1;
		
		list=new Linked_list <Node<NODE_UNIT> >(0, comp_so_far->get_element());
		start_list=list;
		comp_so_far=comp_so_far->next;

		while(comp_so_far != 0) {
			list->next=new Linked_list< Node<NODE_UNIT> > (list, comp_so_far->get_element());
			list=list->next;
			comp_so_far=comp_so_far->next;
			(*size)++;
		}
		
		
		if (num_fully_connected_comps == 0) {
				fully_connected_comps=new Linked_list< Linked_list< Node<NODE_UNIT> > >(0, start_list);
				start_comps=fully_connected_comps;
				comp_sizes=new Linked_list<int> (0, size);
				start_comp_sizes=comp_sizes;
		}
		else {
			fully_connected_comps->next=new Linked_list< Linked_list< Node<NODE_UNIT> > >(fully_connected_comps, start_list);
			fully_connected_comps=fully_connected_comps->next;
			comp_sizes->next=new Linked_list<int> (comp_sizes, size);
			comp_sizes=comp_sizes->next;
		}
		num_fully_connected_comps++;
	}



}


template <class NODE_UNIT>
void Graph<NODE_UNIT>::prune_redundant_fully_connected_comps()
{
	int i, max_size=0;
	BOOL removed;
	Linked_list<int> *other_sizes, *tempsize;
	Linked_list<Linked_list< Node<NODE_UNIT> > > *other_list_copy, *templist;

	//Find the largest component
	comp_sizes=start_comp_sizes;

	while(comp_sizes !=0 ) {
		if ((*comp_sizes->get_element()) > max_size)
			max_size=(*comp_sizes->get_element());
		comp_sizes=comp_sizes->next;
	}

	//Progressively try to merge small components into larger ones (which are completely known
	for(i=1; i<max_size; i++) {
		fully_connected_comps=start_comps;
		comp_sizes=start_comp_sizes;


		while(fully_connected_comps != 0) {
			if ((*comp_sizes->get_element()) == i) {
				other_list_copy=start_comps;
				other_sizes=start_comp_sizes;

				removed=FALSE;

				while((other_list_copy != 0) &&( removed == FALSE)) {
					if ((*other_sizes->get_element()) > i) {
					//Maybe this component (in "fully_connected_comps") is a subset of 
					//the one in other_list_copy
						if (is_list_subset(fully_connected_comps->get_element(), other_list_copy->get_element()) == TRUE) {
							removed=TRUE;
							if (fully_connected_comps->last != 0)
								fully_connected_comps->last->next=fully_connected_comps->next;
							if (fully_connected_comps->next != 0)
								fully_connected_comps->next->last=fully_connected_comps->last;
							templist=fully_connected_comps->last;
							delete fully_connected_comps;
							fully_connected_comps=templist;
							
							if (comp_sizes->last !=0)
								comp_sizes->last->next=comp_sizes->next;
							if (comp_sizes->next != 0)
								comp_sizes->next->last=comp_sizes->last;
							tempsize=comp_sizes->last;
							delete comp_sizes;
							comp_sizes=tempsize;

							num_fully_connected_comps--;		

						}
					}
					other_list_copy=other_list_copy->next;
					other_sizes=other_sizes->next;

				}

			}
			comp_sizes=comp_sizes->next;
			fully_connected_comps=fully_connected_comps->next;
		}

	}


}

template <class NODE_UNIT>
BOOL Graph<NODE_UNIT>::is_list_subset(Linked_list <Node<NODE_UNIT> > *partial_list, Linked_list<Node<NODE_UNIT> > *complete_list)
{
	Linked_list<Node<NODE_UNIT> >  *start_complete;
	BOOL is_subset=TRUE, element_present;


	start_complete=complete_list;

	while( (partial_list != 0) && (is_subset == TRUE) ) {
		element_present = FALSE;

		complete_list=start_complete;

		while((complete_list !=0) && (element_present == FALSE)) {
			if (complete_list->get_element()->get_node_num() == partial_list->get_element()->get_node_num())
				element_present = TRUE;
			complete_list=complete_list->next;
		}

		if (element_present == FALSE)
			is_subset=FALSE;

		partial_list=partial_list->next;

	}
	return(is_subset);
}

template <class NODE_UNIT>
void Graph<NODE_UNIT>::delete_list(Linked_list <Node<NODE_UNIT> > *the_list)
{
	Linked_list<Node<NODE_UNIT> > *temp;

	while(the_list != 0)
	{
		temp=the_list->next;
		delete the_list;
		the_list=temp;
	}

}


template <class NODE_UNIT>
Read_Graph<NODE_UNIT>::Read_Graph() 
{
	have_edge_names=FALSE;
	undirected=TRUE;
	num_nodes=0;
	num_edges=0;
	edge_weights=0;
    float_edge_weights=0;

}
	

template <class NODE_UNIT>
Read_Graph<NODE_UNIT>::Read_Graph(BOOL undir)  
{
	have_edge_names=FALSE;
	undirected=undir;
	num_nodes=0;
	num_edges=0;
	edge_weights=0;
    float_edge_weights=0;
}



template <class NODE_UNIT>
Graph<NODE_UNIT>* Read_Graph<NODE_UNIT>::get_graph(char *nodefile, char *edgefile)
{
	string *tempnode, *tempedge;
	
	tempnode = new string(nodefile);
    tempedge=new string(edgefile);
    return(get_graph(*tempnode, *tempedge));
    
    delete tempnode;
    delete tempedge;
}


template <class NODE_UNIT>
Graph<NODE_UNIT>* Read_Graph<NODE_UNIT>::get_graph(string nodefile,  string edgefile)
{
	int i, j;
	Node<NODE_UNIT> *node1, *node2;
	Edge<NODE_UNIT> *new_edge;
	
	if (edge_weights !=0) delete[] edge_weights;
    if(float_edge_weights != 0) delete[] float_edge_weights;
    
	node_filename = nodefile;
	edge_filename = edgefile;
    
    
	//Read in the nodes (as NODE_UNIT array) and store in units
	if (get_node_list() != 0) {
		cerr<<"Error: invalid node file\n";
		return(0);
	}
    
	//Create the nodes
	the_nodes=new Node<NODE_UNIT>* [num_nodes];
    
	//cout<<"Creating node array\n"<<flush;
	
	for(i=0; i<num_nodes; i++)
 		the_nodes[i]=new Node<NODE_UNIT>(i, &units[i], undirected);
	//	the_nodes[i]=new Node<NODE_UNIT>(i, &units[i]);
	
	//cout<<"Done\n"<<flush;
	
	the_graph=new Graph<NODE_UNIT>(num_nodes, the_nodes);
    
    
	//If edge file fails, we have a graph w/o edges
	//Otherwise, we return a paired list of node *elements*:
	//We search the graph to find the corresponding node and add the edge
	get_edge_list();
	//cout<<"Done with edge list\n"<<flush;
	
    
	for(i=0; i<num_edges; i++) {
		//cout<<"Looking for nodes "<<node_edges[i][0].get_name()<<": "<<node_edges[i][1].get_name()<<endl;
        //j=0;
		//while (!((*the_nodes[j]->get_element()) == node_edges[i][0])) {j++;};
		//node1=the_nodes[j];
        node1=find_node(&node_edges[i][0]);
		
        //j=0;
		//while (!((*the_nodes[j]->get_element()) == node_edges[i][1])) {j++;};
		//node2=the_nodes[j];
        node2=find_node(&node_edges[i][1]);
        
		if (edge_weights ==0) {
            if (float_edge_weights != 0) {
                new_edge=new Edge<NODE_UNIT>(node1, node2, float_edge_weights[i]);
            }
            else {
                if (have_edge_names == FALSE)
                    new_edge=new Edge<NODE_UNIT>(node1, node2);
                else
                    new_edge=new Edge<NODE_UNIT>(node1, node2, edge_names[i]);
            }
		}
		else
			new_edge=new Edge<NODE_UNIT>(node1, node2, edge_weights[i]);
		//cout<<"Creating edge "<<new_edge->get_int_weight()<<endl;
		node1->add_edge(node2, new_edge);
	}
	//cout<<"Done assigning edges\n";
	return(the_graph);
    
}




template <class NODE_UNIT>
Read_Graph<NODE_UNIT>::~Read_Graph()
{
	int i;

	for(i=0; i<num_edges; i++) {
		delete[] node_edges[i];
	}

	if (num_edges > 0) {
		delete[] node_edges;
		if (have_edge_names == TRUE) 
			delete[] edge_names;
	}
}


template <class NODE_UNIT>
Node<NODE_UNIT>* Read_Graph<NODE_UNIT>::find_node(NODE_UNIT *match)
{
    int i;
    
    i=0;
    while (!((*the_nodes[i]->get_element()) == (*match))) {i++;};
    
    return(the_nodes[i]);
}


template <class NODE_UNIT>
EdgePair<NODE_UNIT>& EdgePair<NODE_UNIT>::operator= (EdgePair<NODE_UNIT> &assign_from)
{
    id1=assign_from.id1;
    id2=assign_from.id2;
    empty=assign_from.empty;
    return(*this);
}

template <class NODE_UNIT>
Read_EdgeFile_Graph<NODE_UNIT>::Read_EdgeFile_Graph()
{
    have_edge_names=FALSE;
    undirected=TRUE;
    num_nodes=0;
    num_edges=0;
    header=FALSE;
  }


template <class NODE_UNIT>
Read_EdgeFile_Graph<NODE_UNIT>::Read_EdgeFile_Graph(BOOL undir)
{
    have_edge_names=FALSE;
    undirected=undir;
    num_nodes=0;
    num_edges=0;
    header=FALSE;
}

template <class NODE_UNIT>
Read_EdgeFile_Graph<NODE_UNIT>::Read_EdgeFile_Graph(BOOL undir, BOOL head)
{
    have_edge_names=FALSE;
    undirected=undir;
    num_nodes=0;
    num_edges=0;
    header=head;
}



template <class NODE_UNIT>
Graph<NODE_UNIT>* Read_EdgeFile_Graph<NODE_UNIT>::get_graph(string edgefile)
{
    int i, node_id1, node_id2;
    string dummy;
    EdgePair<NODE_UNIT> *read_edge;
    NODE_UNIT *my_unit;
    Edge<NODE_UNIT> *new_edge;
    std::map<NODE_UNIT, int> node_hash;
    typename std::map<NODE_UNIT, int>::iterator it;
    List< EdgePair<NODE_UNIT> > edge_list;
    
    
    node_hash.clear();
    
    num_nodes=0;
    num_edges=0;
    
    edgein.open(edgefile.c_str());
    
    if(edgein.fail())
        cerr<<"Error: invalid interaction file\n";
    
    if (header ==TRUE)
        getline(edgein, dummy);
    while(!(edgein.eof())) {
        read_edge=new EdgePair<NODE_UNIT>;
        read_edge_pair(*read_edge);
        if (read_edge->empty == FALSE) {
        
            if(node_hash.find(read_edge->id1) == node_hash.end()){
                node_hash[read_edge->id1] = num_nodes;
                num_nodes++;
            }
            
            if(node_hash.find(read_edge->id2) == node_hash.end()){
                node_hash[read_edge->id2] = num_nodes;
                num_nodes++;
            }
            
            edge_list.add_to_list(*read_edge);
        }
        
    }
    
    edgein.close();
    
    the_nodes=new Node<NODE_UNIT>* [num_nodes];
    
    //cout<<"Creating a graph of "<<num_nodes<<" nodes\n";
    
    for (it=node_hash.begin(); it!=node_hash.end(); ++it) {
        my_unit=new NODE_UNIT();
        (*my_unit)=it->first;
        the_nodes[it->second]=new Node<NODE_UNIT>(it->second, my_unit, undirected);
    }

    
    the_graph=new Graph<NODE_UNIT>(num_nodes, the_nodes);
    
    for(i=0; i<edge_list.get_list_length(); i++) {
        node_id1=node_hash[edge_list.get_nth_element(i)->item()->id1];
        node_id2=node_hash[edge_list.get_nth_element(i)->item()->id2];
        
        //cout<<"Edge count "<<i<<" adding edge for ids "<<node_id1<<" and "<<node_id2<<endl;
       
        if (the_graph->get_node(node_id1)->edge_exists(the_graph->get_node(node_id2)) == FALSE) {
            new_edge=new Edge<NODE_UNIT>(the_graph->get_node(node_id1), the_graph->get_node(node_id2));
            the_graph->get_node(node_id1)->add_edge(the_graph->get_node(node_id2), new_edge);
            num_edges++;
        }
        
    }
    
    //cout<<"Done assigning edges\n";
    return(the_graph);
    
}





template <class NODE_UNIT>
Read_EdgeFile_Graph<NODE_UNIT>::~Read_EdgeFile_Graph()
{
   
}


