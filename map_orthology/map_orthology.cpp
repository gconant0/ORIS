#include <iostream>
#include <fstream>
#include <map>
#include "graph.cpp"
#include "linked_list.cpp"

using namespace::std;

#define MAX_DOUBLE 1e308

//ORPHANS WITH ONLY SELF GENOME MATCHES???

class Ensembl_gene {
public:
	BOOL no_order;
	Ensembl_gene() {no_order=TRUE; chrom=-1;};
	Ensembl_gene& operator=(Ensembl_gene &assign_from);
	int chrom, position;
	char name[50];
	
};

class Orthology_node {
public:	
	int region_size;
	BOOL left_complete, right_complete, in_ortho_list, orphan, no_order;
	
	Orthology_node();
	Orthology_node(Ensembl_gene *new_gene);
	Orthology_node& operator=(Orthology_node &assign_from);
	void add_gene(Ensembl_gene *new_gene);
	int get_chromosome()  {return(chrom);};
	int get_position()	{return(position);};
	int get_num_genes()	{return(num_genes);};
	int get_node_num()	{return(node_num);};
	int get_gene_pos()	{return(gene_pos);};
	int get_organism()	{return(organism);};
	Orthology_node* get_ortholog()	{return(my_ortholog);};
	BOOL omit_obj()		{return(omit);};
	Ensembl_gene * get_gene(int gene_num);
	void set_position(int new_pos) {position=new_pos;};
	void set_chromosome(int new_chrom) {chrom=new_chrom;};
	void set_omit()   {omit=TRUE;};
	void set_node_num(int num)	{node_num=num;};
	void set_gene_pos(int pos)	{gene_pos=pos;};
	void set_organism(int org)  {organism=org;};
	void set_ortholog(Orthology_node *new_ortho)  {my_ortholog=new_ortho;};
	~Orthology_node();
	
protected:
	int num_genes, chrom, position, node_num, gene_pos, organism;
	BOOL omit;
	Orthology_node *my_ortholog;
	Ensembl_gene **genes;	
};

class Ortholog_pair {
public:
	Ortholog_pair();
	Ortholog_pair& operator=(Ortholog_pair &assign_from);
	int operator==(Ortholog_pair test)  {if ((ortho_pair[0] == test.ortho_pair[0]) && (ortho_pair[1] == test.ortho_pair[1])) return(0); else return(-1);};
	Orthology_node *ortho_pair[2];
};


BOOL read_genefile(char *filename, Orthology_node *&initial_genes, int &num_genes);
void assign_positions(int num_genes, Orthology_node **gene_list);
void collapse_tandems(int num_genes, int &nnew_num_objs, Orthology_node *&new_list, Graph<Orthology_node> *the_graph);
BOOL read_parafile(char *filename, Graph<Orthology_node> *the_graph, double ks_thresh);
void merge_nodes (Orthology_node *node1, Orthology_node *node2);
BOOL read_orthofile(char *filename, Graph<Orthology_node> *the_graph, Graph<Orthology_node> *single_gene_graph, BOOL use_ka);

void find_ortholog_self_matches(Graph<Orthology_node> * the_graph, List<Ortholog_pair> *orthologs, double ks_thresh);

void find_more_anchors(Graph<Orthology_node> * the_graph, List<Ortholog_pair> *orthologs, double ks_thresh);


void clear_ortho_edges(Graph<Orthology_node> *the_graph, List<Ortholog_pair> *orthologs);
void remove_orphans(Graph<Orthology_node> *the_graph, Orthology_node **full_list, int full_size, int &new_num_objs, Orthology_node **&new_list, int &num_orphans, Orthology_node **&orphan_list);
void number_regions(Orthology_node **gene_list, int num_objs);
void rank_edges(Graph<Orthology_node> *the_graph, Graph<Orthology_node> *single_gene_graph);



int main (int argc, char **argv) 
{
	int i, j, k, num_genes1, num_genes2, num_reduced1, num_reduced2, chroms1, chroms2, num_completed,
		gene_loc, curr_chrom[2], curr_min[2], curr_max[2], num_in_region, num_curr1, num_curr2, num_orph1, num_orph2, cnt_ortho1, cnt_ortho2;
	char  genefile1[100], genefile2[100], parafile1[100], parafile2[100], orthofile[100];
    string output_header, outfile_tandem, outfile_ortho, outfile_amb1, outfile_amb2, outfile_orph1, outfile_orph2, ks_string;
	double ks_thresh, val1, val2, tandem_thresh;
    BOOL use_ka=FALSE;
	ifstream infile;
	ofstream outfile;
	Orthology_node *initial_genes1, *initial_genes2, *full_initial, *reduced_genes1, *reduced_genes2, **reduced_p1, **reduced_p2, **curr_genes1, **curr_genes2, **orphans1, **orphans2;
	Ortholog_pair *anchor, *a_pair, pair1, pair2;
	Node<Orthology_node> *node1, *node2, *save_node, *ortho_node1, *ortho_node2;
	Graph<Orthology_node> *graph1, *graph2, *the_graph, *single_gene_graph;
	List<Ortholog_pair> *orthologs, *new_orthologs, *completed_orthologs;	

	if (argc>6) {

		//cout<<"Reading file "<<argv[1]<<endl;
		strcpy(genefile1, argv[1]);
		strcpy(genefile2, argv[2]);
		strcpy(parafile1, argv[3]);
		strcpy(parafile2, argv[4]);
		strcpy(orthofile, argv[5]);
        output_header=argv[6];
        if ((argv[7][1] == 'a') || (argv[7][1] == 'A')) {
            use_ka=TRUE;
        }
        ks_string=argv[7];
        ks_thresh=string_to_float(ks_string.substr(3, ks_string.length()).c_str());
    
        //tandem search uses a different thresh that is more strict.
        tandem_thresh = ks_thresh / 2.0;

		
		if (read_genefile(genefile1, initial_genes1, num_genes1) == FALSE) {
			cerr<<"Error reading gene file "<<genefile1<<endl;
			return(-1);
		}
		
		if (read_genefile(genefile2, initial_genes2, num_genes2) == FALSE) {
			cerr<<"Error reading gene file "<<genefile2<<endl;
			return(-1);
		}
		
		
		curr_genes1=new Orthology_node * [num_genes1];
		for(i=0; i<num_genes1; i++)
			curr_genes1[i]=&initial_genes1[i];
		
		
		curr_genes2=new Orthology_node * [num_genes2];
		for(i=0; i<num_genes2; i++)
			curr_genes2[i]=&initial_genes2[i];
		
		
		full_initial=new Orthology_node [num_genes1+num_genes2];
		for(i=0; i<num_genes1; i++)
			full_initial[i]=initial_genes1[i];
		for(i=0; i<num_genes2; i++)
			full_initial[i+num_genes1]=initial_genes2[i];
		single_gene_graph=new Graph<Orthology_node>(num_genes1+num_genes2, TRUE);
		for(i=0; i<num_genes1+num_genes2; i++)
			single_gene_graph->get_node(i)->set_element(&full_initial[i]);


		assign_positions(num_genes1, curr_genes1);
		assign_positions(num_genes2, curr_genes2);
		
		delete[] curr_genes1;
		delete[] curr_genes2;
		
		chroms1=initial_genes1[num_genes1-1].get_chromosome();
		chroms2=initial_genes2[num_genes2-1].get_chromosome();
		
		cout<<"Species 1 has "<<chroms1<<endl;
		cout<<"Species 2 has "<<chroms2<<endl;
		
		graph1=new Graph<Orthology_node>(num_genes1, TRUE);
		graph2=new Graph<Orthology_node>(num_genes2, TRUE);
		
		for(i=0; i<num_genes1; i++)
			graph1->get_node(i)->set_element(&initial_genes1[i]);
		
		for(i=0; i<num_genes2; i++)
			graph2->get_node(i)->set_element(&initial_genes2[i]);
		
		
		//changed to tandem_thresh
		read_parafile(parafile1, graph1, ks_thresh);
		//read_parafile(parafile1, graph1, tandem_thresh);
		cout<<"Read first paralog file\n"<<flush;
		read_parafile(parafile2, graph2, ks_thresh);
		//read_parafile(parafile2, graph2, tandem_thresh);
		
		collapse_tandems(num_genes1, num_reduced1, reduced_genes1, graph1);
		collapse_tandems(num_genes2, num_reduced2, reduced_genes2, graph2);
		
		
		
		cout<<"Genome 1 reduced from "<<num_genes1<<" to "<<num_reduced1<<endl;
		cout<<"Genome 2 reduced from "<<num_genes2<<" to "<<num_reduced2<<endl;
		delete[] initial_genes1;
		delete[] initial_genes2;
		
		curr_genes1=new Orthology_node * [num_reduced1];
		for(i=0; i<num_reduced1; i++) {
			curr_genes1[i]=&reduced_genes1[i];
			curr_genes1[i]->set_organism(0);
		}
		
		curr_genes2=new Orthology_node * [num_reduced2];
		for(i=0; i<num_reduced2; i++) {
			curr_genes2[i]=&reduced_genes2[i];
			curr_genes2[i]->set_organism(1);
		}
		assign_positions(num_reduced1, curr_genes1);
		assign_positions(num_reduced2, curr_genes2);
		
        outfile_tandem=output_header + "_tandem_groups.txt";
		outfile.open(outfile_tandem.c_str());
		for(i=0; i<num_reduced1; i++) {
			for(j=0; j<reduced_genes1[i].get_num_genes(); j++)
				outfile<<reduced_genes1[i].get_gene(j)->name<<"\t"<<i<<endl;
		}
		for(i=0; i<num_reduced2; i++) {
			for(j=0; j<reduced_genes2[i].get_num_genes(); j++)
				outfile<<reduced_genes2[i].get_gene(j)->name<<"\t"<<i<<endl;
		}
		outfile.close();
		
		for(i=0; i<num_reduced1; i++) {
			if (reduced_genes1[i].get_num_genes() >1) {
				cout<<"Tandem: ";
				for (j=0; j<reduced_genes1[i].get_num_genes(); j++)
					cout<<reduced_genes1[i].get_gene(j)->name<<"\t";
				cout<<endl;
			}
			
		}
		
		for(i=0; i<num_reduced2; i++) {
			if (reduced_genes2[i].get_num_genes() >1) {
				cout<<"Tandem: ";
				for (j=0; j<reduced_genes2[i].get_num_genes(); j++)
					cout<<reduced_genes2[i].get_gene(j)->name<<"\t";
				cout<<endl;
			}
			
		}
		
		the_graph = new Graph<Orthology_node>(num_reduced1+num_reduced2, TRUE);
		
		for(i=0; i<num_reduced1; i++) {
			reduced_genes1[i].set_node_num(i);
			the_graph->get_node(i)->set_element(&reduced_genes1[i]);
		}
		
		for(i=0; i<num_reduced2; i++) {
			reduced_genes2[i].set_node_num(i+num_reduced1);
			the_graph->get_node(i+num_reduced1)->set_element(&reduced_genes2[i]);
		}
		
		
		
		reduced_p1=new Orthology_node*[num_reduced1];
		reduced_p2=new Orthology_node*[num_reduced2];
		
		for(i=0; i<the_graph->get_num_nodes(); i++) {
				if (the_graph->get_node(i)->get_element()->get_organism() == 0)
					reduced_p1[the_graph->get_node(i)->get_element()->get_gene_pos()]=the_graph->get_node(i)->get_element();
				else
					reduced_p2[the_graph->get_node(i)->get_element()->get_gene_pos()]=the_graph->get_node(i)->get_element();
		}
		
		delete[] reduced_genes1;
		delete[] reduced_genes2;
		
		for(i=0; i<the_graph->get_num_nodes(); i++) 
			the_graph->get_node(i)->get_element()->set_node_num(i);
		
		if (read_orthofile(orthofile, the_graph, single_gene_graph, use_ka) == FALSE) {
			cerr<<"Error reading orthology file "<<orthofile<<endl;
			return(-1);
		}
		
		orthologs=new List<Ortholog_pair>();
		new_orthologs=new List<Ortholog_pair>();
		completed_orthologs=new List<Ortholog_pair>();
		
		find_ortholog_self_matches(the_graph, orthologs, ks_thresh);
		
		clear_ortho_edges(the_graph, orthologs);
		delete[] curr_genes1;
		delete[] curr_genes2;
		num_curr1=num_reduced1;
		num_curr2=num_reduced2;
		remove_orphans(the_graph, reduced_p1, num_reduced1, num_curr1, curr_genes1, num_orph1, orphans1);
		remove_orphans(the_graph, reduced_p2, num_reduced2, num_curr2, curr_genes2, num_orph2, orphans2);
		assign_positions(num_curr1, curr_genes1);
		assign_positions(num_curr2, curr_genes2);
		
		
		rank_edges(the_graph, single_gene_graph);
		cout<<"Found "<<orthologs->get_list_length()<<" ortholog pairs in first pass\n";
		
		//HANDLE DCB  verses BCD problem!!!!!

		for(k=0; k<3; k++) {
			do {
				num_completed=0;
				for(i=0; i<orthologs->get_list_length(); i++) {
					anchor=orthologs->get_nth_element(i)->item();
					
					//Anchor right
					if ((anchor->ortho_pair[0]->get_gene_pos() != (num_curr1-1)) && (anchor->ortho_pair[0]->right_complete == FALSE)) {
						if ((curr_genes1[anchor->ortho_pair[0]->get_gene_pos()+1]->get_chromosome() == anchor->ortho_pair[0]->get_chromosome()) &&
							(curr_genes1[anchor->ortho_pair[0]->get_gene_pos()+1]->get_position() == (anchor->ortho_pair[0]->get_position()+1))) {
							node1=the_graph->get_node(curr_genes1[anchor->ortho_pair[0]->get_gene_pos()+1]->get_node_num());
							
							if (node1->get_element()->get_ortholog() != 0) {
								if ((node1->get_element()->get_ortholog()->get_chromosome() == anchor->ortho_pair[1]->get_chromosome()) && 
									((node1->get_element()->get_ortholog()->get_position() == (anchor->ortho_pair[1]->get_position()+1)) || 
									 (node1->get_element()->get_ortholog()->get_position() ==(anchor->ortho_pair[1]->get_position()-1)))) {
									anchor->ortho_pair[0]->right_complete=TRUE;
									node1->get_element()->left_complete=TRUE;
								}
							}
							else {
								node1->reset_edges();
								pair1.ortho_pair[0]=0;
								pair2.ortho_pair[0]=0;

								for(j=0; j<node1->get_num_edges(); j++) {
									node2=node1->get_next_edge();
									
									if (node2->get_element()->get_chromosome() == anchor->ortho_pair[1]->get_chromosome()) {
										if (node2->get_element()->get_position() == (anchor->ortho_pair[1]->get_position()+1)) {
											val1=node2->get_edge(node1)->get_double_weight();
											ortho_node1=node2;
											pair1.ortho_pair[0]=node1->get_element();
											pair1.ortho_pair[1]=node2->get_element();
										} 
										if (node2->get_element()->get_position() == (anchor->ortho_pair[1]->get_position()-1)) {
											val2=node2->get_edge(node1)->get_double_weight();
											ortho_node2=node2;
											pair2.ortho_pair[0]=node1->get_element();
											pair2.ortho_pair[1]=node2->get_element();
										}
									}
								}
								
								a_pair=0;
								
								if ((pair1.ortho_pair[0] != 0) && (pair2.ortho_pair[0] !=0)) {
									if (val2 < val1) {
										a_pair = &pair2;
										ortho_node1=ortho_node2;
									}
									else
										a_pair = &pair1;
								}
								else 
								{
									if (pair1.ortho_pair[0] != 0) 
										a_pair = &pair1;
									
									else  if (pair2.ortho_pair[0] != 0) { 
										a_pair = &pair2;
										ortho_node1=ortho_node2;
									}
								}
								if (a_pair != 0) {
									if ((k== 1) ||  (node1->get_edge(ortho_node1) == node1->best_edge)) {
										if (a_pair->ortho_pair[1]->in_ortho_list == FALSE) {
											anchor->ortho_pair[0]->right_complete=TRUE;
											if (node1->get_element()->in_ortho_list == FALSE) {
												a_pair->ortho_pair[0]->left_complete=TRUE;
												a_pair->ortho_pair[1]->left_complete=TRUE;
												a_pair->ortho_pair[0]->right_complete=FALSE;
												a_pair->ortho_pair[1]->right_complete=FALSE;
												a_pair->ortho_pair[0]->in_ortho_list=TRUE;
												a_pair->ortho_pair[1]->in_ortho_list=TRUE;
												new_orthologs->add_to_list(*a_pair);
												a_pair->ortho_pair[0]->set_ortholog(a_pair->ortho_pair[1]);
												a_pair->ortho_pair[1]->set_ortholog(a_pair->ortho_pair[0]);
												num_completed++;
											}
										}
									}
								}
							}
						}
					}
				
					//Anchor left
					if ((anchor->ortho_pair[0]->get_gene_pos() != 0) && (anchor->ortho_pair[0]->left_complete == FALSE)) {
						if ((curr_genes1[anchor->ortho_pair[0]->get_gene_pos()-1]->get_chromosome() == anchor->ortho_pair[0]->get_chromosome()) &&
							(curr_genes1[anchor->ortho_pair[0]->get_gene_pos()-1]->get_position() == (anchor->ortho_pair[0]->get_position()-1))) {
							node1=the_graph->get_node(curr_genes1[anchor->ortho_pair[0]->get_gene_pos()-1]->get_node_num());
							
							if (node1->get_element()->get_ortholog() != 0) {
								if ((node1->get_element()->get_ortholog()->get_chromosome() == anchor->ortho_pair[1]->get_chromosome()) && 
									((node1->get_element()->get_ortholog()->get_position() == (anchor->ortho_pair[1]->get_position()+1)) || 
									 (node1->get_element()->get_ortholog()->get_position() ==(anchor->ortho_pair[1]->get_position()-1)))) {
									anchor->ortho_pair[0]->left_complete=TRUE;
									node1->get_element()->right_complete=TRUE;
								}
							}
							else {
								node1->reset_edges();
								pair1.ortho_pair[0]=0;
								pair2.ortho_pair[0]=0;
								
								for(j=0; j<node1->get_num_edges(); j++) {
									node2=node1->get_next_edge();
									
									if (node2->get_element()->get_chromosome() == anchor->ortho_pair[1]->get_chromosome()) {
										if (node2->get_element()->get_position() == (anchor->ortho_pair[1]->get_position()+1)) {
											val1=node2->get_edge(node1)->get_double_weight();
											ortho_node1=node2;
											pair1.ortho_pair[0]=node1->get_element();
											pair1.ortho_pair[1]=node2->get_element();
										} 
										if (node2->get_element()->get_position() == (anchor->ortho_pair[1]->get_position()-1)) {
											val2=node2->get_edge(node1)->get_double_weight();
											ortho_node2=node2;
											pair2.ortho_pair[0]=node1->get_element();
											pair2.ortho_pair[1]=node2->get_element();
										}
									}
								}
								
								a_pair=0;
								
								if ((pair1.ortho_pair[0] != 0) && (pair2.ortho_pair[0] !=0)) {
									if (val2 < val1)  {
										a_pair = &pair2;
										ortho_node1=ortho_node2;
									}
									else
										a_pair = &pair1;
								}
								else 
								{
									if (pair1.ortho_pair[0] != 0) 
										a_pair = &pair1;
									
									else if (pair2.ortho_pair[0] != 0) {
										a_pair = &pair2;
										ortho_node1=ortho_node2;
									}
								}
								if (a_pair != 0) {
									if ((k== 1) ||  (node1->get_edge(ortho_node1) == node1->best_edge)) {
										if (a_pair->ortho_pair[1]->in_ortho_list == FALSE) {
											anchor->ortho_pair[0]->left_complete=TRUE;
											if (node1->get_element()->in_ortho_list == FALSE) {
												a_pair->ortho_pair[0]->left_complete=FALSE;
												a_pair->ortho_pair[1]->left_complete=FALSE;
												a_pair->ortho_pair[0]->right_complete=TRUE;
												a_pair->ortho_pair[1]->right_complete=TRUE;
												a_pair->ortho_pair[0]->in_ortho_list=TRUE;
												a_pair->ortho_pair[1]->in_ortho_list=TRUE;
												new_orthologs->add_to_list(*a_pair);
												a_pair->ortho_pair[0]->set_ortholog(a_pair->ortho_pair[1]);
												a_pair->ortho_pair[1]->set_ortholog(a_pair->ortho_pair[0]);
												num_completed++;
											}
										}
									}
								}
							}
						}
					}
				}		
				
				clear_ortho_edges(the_graph, new_orthologs);
				
				for(i=0; i<orthologs->get_list_length(); i++) {
					if ((orthologs->get_nth_element(i)->item()->ortho_pair[0]->left_complete == FALSE) ||
						(orthologs->get_nth_element(i)->item()->ortho_pair[0]->right_complete == FALSE)) 
							new_orthologs->add_to_list(*orthologs->get_nth_element(i)->item());
					else
							completed_orthologs->add_to_list(*orthologs->get_nth_element(i)->item());
				}
			
				cout<<"New orthologs "<<new_orthologs->get_list_length()<<" and "<<completed_orthologs->get_list_length()<<" completed orthologs\n";
				
				delete orthologs;
				orthologs=new_orthologs;
				new_orthologs=new List<Ortholog_pair>();
				
				if (k<1) {
				find_ortholog_self_matches(the_graph, orthologs, ks_thresh); 
				} else {
						//third round: find *not* one-to-one anchors
				find_more_anchors(the_graph, orthologs, ks_thresh);
				}
				clear_ortho_edges(the_graph, orthologs);
				
				delete[] curr_genes1;
				delete[] curr_genes2;
				remove_orphans(the_graph, reduced_p1, num_reduced1, num_curr1, curr_genes1, num_orph1, orphans1);
				remove_orphans(the_graph, reduced_p2, num_reduced2, num_curr2, curr_genes2, num_orph2, orphans2);
				assign_positions(num_curr1, curr_genes1);
				assign_positions(num_curr2, curr_genes2);
				cout<<"There are "<<orthologs->get_list_length()<<" ortholog pairs left to analyze\n";
				
			} while(num_completed > 0);
		}
		
		number_regions(curr_genes1, num_curr1);
        outfile_ortho=output_header + "_found_orthologs.txt";
		outfile.open(outfile_ortho.c_str());
		cnt_ortho1=cnt_ortho2=0;
		
		for(i=0; i<completed_orthologs->get_list_length(); i++) {
			anchor=completed_orthologs->get_nth_element(i)->item();
			cnt_ortho1+=anchor->ortho_pair[0]->get_num_genes();
			cnt_ortho2+=anchor->ortho_pair[1]->get_num_genes();
			
			
			outfile<<anchor->ortho_pair[0]->region_size<<"\t"
					<<anchor->ortho_pair[0]->get_chromosome()<<"\t"<<anchor->ortho_pair[0]->get_position()<<"\t"
					<<anchor->ortho_pair[1]->get_chromosome()<<"\t"<<anchor->ortho_pair[1]->get_position()<<"\t";
			
			for(j=0; j<anchor->ortho_pair[0]->get_num_genes(); j++)
				outfile<<anchor->ortho_pair[0]->get_gene(j)->name<<"\t";
			for(j=0; j<anchor->ortho_pair[1]->get_num_genes()-1; j++)
				outfile<<anchor->ortho_pair[1]->get_gene(j)->name<<"\t";
			outfile<<anchor->ortho_pair[1]->get_gene(anchor->ortho_pair[1]->get_num_genes()-1)->name<<"\n";

		}
		
		for(i=0; i<orthologs->get_list_length(); i++) {
			anchor=orthologs->get_nth_element(i)->item();
			cnt_ortho1+=anchor->ortho_pair[0]->get_num_genes();
			cnt_ortho2+=anchor->ortho_pair[1]->get_num_genes();
			
			outfile<<anchor->ortho_pair[0]->region_size<<"\t"
				<<anchor->ortho_pair[0]->get_chromosome()<<"\t"<<anchor->ortho_pair[0]->get_position()<<"\t"
				<<anchor->ortho_pair[1]->get_chromosome()<<"\t"<<anchor->ortho_pair[1]->get_position()<<"\t";
			
			for(j=0; j<anchor->ortho_pair[0]->get_num_genes(); j++)
				outfile<<anchor->ortho_pair[0]->get_gene(j)->name<<"\t";
			for(j=0; j<anchor->ortho_pair[1]->get_num_genes()-1; j++)
				outfile<<anchor->ortho_pair[1]->get_gene(j)->name<<"\t";
			outfile<<anchor->ortho_pair[1]->get_gene(anchor->ortho_pair[1]->get_num_genes()-1)->name<<"\n";
			
		}
		
		outfile.close();
		outfile_orph1=output_header+"_sp1_orphans.txt";
		outfile.open(outfile_orph1.c_str());
		for(i=0; i<num_orph1; i++) {
			for(j=0; j<orphans1[i]->get_num_genes()-1; j++) 
					outfile<<orphans1[i]->get_gene(j)->name<<"\t";
			outfile<<orphans1[i]->get_gene(orphans1[i]->get_num_genes()-1)->name<<"\n";
		}
		outfile.close();
		
        outfile_orph2=output_header+"_sp2_orphans.txt";
        outfile.open(outfile_orph2.c_str());
		for(i=0; i<num_orph2; i++) {
			for(j=0; j<orphans2[i]->get_num_genes()-1; j++) 
				outfile<<orphans2[i]->get_gene(j)->name<<"\t";
			outfile<<orphans2[i]->get_gene(orphans2[i]->get_num_genes()-1)->name<<"\n";
		}
		outfile.close();
		
		outfile_amb1=output_header+"_sp1_ambig.txt";
		outfile.open(outfile_amb1.c_str());
		for(i=0; i<the_graph->get_num_nodes(); i++) {
			if (the_graph->get_node(i)->get_element()->get_organism() == 0) {
				if ((the_graph->get_node(i)->get_element()->in_ortho_list == FALSE) && 
					(the_graph->get_node(i)->get_element()->orphan == FALSE)) {
					node1=the_graph->get_node(i);
					for(j=0; j<node1->get_element()->get_num_genes()-1; j++) 
						outfile<<node1->get_element()->get_gene(j)->name<<"\t";
					outfile<<node1->get_element()->get_gene(node1->get_element()->get_num_genes()-1)->name<<"\n";						
				}
			}
			
		}
		
		outfile.close();
		
		
        outfile_amb2=output_header+"_sp2_ambig.txt";
        outfile.open(outfile_amb2.c_str());
		for(i=0; i<the_graph->get_num_nodes(); i++) {
			if (the_graph->get_node(i)->get_element()->get_organism() == 1) {
				if ((the_graph->get_node(i)->get_element()->in_ortho_list == FALSE) && 
					(the_graph->get_node(i)->get_element()->orphan == FALSE)) {
					node1=the_graph->get_node(i);
					for(j=0; j<node1->get_element()->get_num_genes()-1; j++) 
						outfile<<node1->get_element()->get_gene(j)->name<<"\t";
					outfile<<node1->get_element()->get_gene(node1->get_element()->get_num_genes()-1)->name<<"\n";						
				}
			}
			
		}
		
		outfile.close();
		
		cout<<"Species 1 has "<<cnt_ortho1<<" of "<<num_genes1<<" in orthology relationships\n";
		cout<<"Species 2 has "<<cnt_ortho2<<" of "<<num_genes2<<" in orthology relationships\n";

		
		delete[] reduced_p1;
		delete[] reduced_p2;
		delete the_graph;
		return(0);
	}
	else {
		cerr<<"Usage: map_orthology <genefile1> <genefile2> <parafile1> <parafile2> <orthofile>\n";
		return(-1);
	}
}



Ensembl_gene& Ensembl_gene::operator=(Ensembl_gene &assign_from)
{
	strcpy(name, assign_from.name);
	chrom=assign_from.chrom;
	position=assign_from.position;
	return(*this);
}

Orthology_node::Orthology_node()
{
	num_genes=0;
	genes=0;
	omit=FALSE;
	left_complete=right_complete=in_ortho_list=orphan=FALSE;
	no_order=TRUE;
	my_ortholog=0;
}


Orthology_node::Orthology_node(Ensembl_gene *new_gene)
{
	num_genes=1;
	genes=new Ensembl_gene * [1];
	omit=FALSE;
	
	genes[0] = new Ensembl_gene;
	(*genes[0])=(*new_gene);
	chrom=new_gene->chrom;
	left_complete=right_complete=in_ortho_list=orphan=FALSE;
	no_order=TRUE;
	my_ortholog=0;
}



Orthology_node & Orthology_node::operator=(Orthology_node &assign_from)
{
	int i;
	
	if (num_genes != assign_from.get_num_genes()) {
		if (num_genes !=0) {
			for(i=0; i<num_genes; i++)
				delete genes[i];
			delete[] genes;
		}
		genes = new Ensembl_gene * [assign_from.get_num_genes()];
		num_genes=assign_from.get_num_genes();
		for(i=0; i<num_genes; i++)
			genes[i] = new Ensembl_gene;
	}
	
	chrom=assign_from.get_chromosome();
	position=assign_from.get_position();
	gene_pos=assign_from.get_gene_pos();
	organism=assign_from.get_organism();
	no_order=assign_from.no_order;
	my_ortholog=0;
	
	for(i=0; i<num_genes; i++)
		(*genes[i]) = (*assign_from.get_gene(i));
	return(*this);
}


void Orthology_node::add_gene(Ensembl_gene *new_gene)
{
	int i;
	
	Ensembl_gene **temp=0;
	
	if (num_genes > 0)
		temp=genes;
	
	num_genes++;
	genes=new Ensembl_gene * [num_genes];
	
	for(i=0; i<num_genes-1; i++)
		genes[i] = temp[i];
	
	genes[num_genes-1] = new Ensembl_gene;
	
	(*genes[num_genes-1])=(*new_gene);
	
	if (temp != 0)
		delete[] temp;
}


Ensembl_gene * Orthology_node::get_gene(int gene_num)
{
	if ((gene_num >=0) && (gene_num < num_genes)) return(genes[gene_num]);
	else return(0);
}


Orthology_node::~Orthology_node()
{
	int i;
	
	if (num_genes >0) {
		for(i=0; i<num_genes; i++)
			delete genes[i];
		delete[] genes;
	}
}

Ortholog_pair::Ortholog_pair()
{
	ortho_pair[0]=ortho_pair[1]=0; }

Ortholog_pair & Ortholog_pair::operator=(Ortholog_pair &assign_from)
{
	ortho_pair[0]=assign_from.ortho_pair[0]; 
	ortho_pair[1]=assign_from.ortho_pair[1]; 
	return(*this);	
}


BOOL read_genefile(char *filename, Orthology_node *&initial_genes, int &num_genes)
{
	int i, pos;
	ifstream genefile;
	Ensembl_gene *curr_gene;
	Linked_list<Ensembl_gene> *readgenes, *startgenes;

	num_genes=0;

	genefile.open(filename);
	
	if (genefile.fail())
	{
		cerr<<"Can't find file "<<filename<<endl;
			return(FALSE);
	}
	
	curr_gene = new Ensembl_gene;
	genefile>>curr_gene->name>>curr_gene->chrom>>pos;
	if (pos != -1) curr_gene->no_order=FALSE;
	
	readgenes=new Linked_list<Ensembl_gene> (0, curr_gene);
	num_genes++;
	startgenes=readgenes;
	
		
	while (!genefile.eof()) {
		
		curr_gene = new Ensembl_gene;
		genefile>>curr_gene->name>>curr_gene->chrom>>pos;
		if (pos != -1) curr_gene->no_order=FALSE;

	  if (strlen(curr_gene->name) > 0) {
	    readgenes->next=new Linked_list<Ensembl_gene> (readgenes, curr_gene);
	    readgenes=readgenes->next;
	   
	    num_genes++;
	  }
	}

	genefile.close();

	readgenes=startgenes;
	initial_genes = new Orthology_node[num_genes];
	
	for(i=0; i<num_genes; i++) {
		initial_genes[i].add_gene(readgenes->get_element());
		initial_genes[i].set_chromosome(readgenes->get_element()->chrom);
		initial_genes[i].no_order=readgenes->get_element()->no_order;
		startgenes=readgenes->next;
		delete readgenes->get_element();
		delete readgenes;
		readgenes=startgenes;
	}
	
	return(TRUE);
}


void assign_positions(int num_genes, Orthology_node **gene_list)
{
	int i, j, last_chrom;
	
	if (gene_list[0]->no_order == FALSE)
		gene_list[0]->set_position(0);
	else
		gene_list[0]->set_position(-1);
	
	gene_list[0]->set_gene_pos(0);
	last_chrom=gene_list[0]->get_chromosome();
	for(j=0; j<gene_list[0]->get_num_genes(); j++) 
		gene_list[0]->get_gene(j)->position=gene_list[0]->get_position();
	
	
	for(i=1; i<num_genes; i++) {
		gene_list[i]->set_gene_pos(i);
		if (gene_list[i]->no_order == TRUE)
			gene_list[i]->set_position(-1);
		else {
			if (gene_list[i]->get_chromosome() == last_chrom) gene_list[i]->set_position(gene_list[i-1]->get_position()+1);
			else gene_list[i]->set_position(0);
		}
		for(j=0; j<gene_list[i]->get_num_genes(); j++) 
			gene_list[i]->get_gene(j)->position=gene_list[i]->get_position();
			
		last_chrom=gene_list[i]->get_chromosome();
	}
	
}




BOOL read_parafile(char *filename,  Graph<Orthology_node> *the_graph, double ks_thresh)
{
	int i, j, len;
	char name[50], paraname[50], sat[10];
	double ks, ka;
	ifstream genefile;

	genefile.open(filename);
	if (genefile.fail())
		return(FALSE);

	
	if (!genefile.eof()) {
		genefile>>name>>paraname>>ks>>ka>>len>>sat;

		while (!genefile.eof())
		{
			if (strcmp(name, paraname) != 0) {
				i=0;
				while ( (i<the_graph->get_num_nodes()) && (strcmp(the_graph->get_node(i)->get_element()->get_gene(0)->name, name) != 0)) i++;
				
				j=0;
				while ( (j<the_graph->get_num_nodes()) && (strcmp(the_graph->get_node(j)->get_element()->get_gene(0)->name, paraname) != 0)) j++;
				
				if (strcmp(sat, "YES") == 0) ks=-1;
				
				if ((i<the_graph->get_num_nodes()) && (j<the_graph->get_num_nodes())) {		
					if ((ks != -1) && (ks <= ks_thresh)) {
						the_graph->get_node(i)->add_edge(the_graph->get_node(j));
						the_graph->get_node(i)->get_edge(the_graph->get_node(j))->set_double_weight(ks);
					}
				}
				else
				  cerr<<"Error paralogs "<<name<<" and "<<paraname<<" not found in graph\n";
			}
				

			genefile>>name>>paraname>>ks>>ka>>len>>sat;
		}
	}
	genefile.close();
	return(TRUE);
}


void merge_nodes (Orthology_node *node1, Orthology_node *node2)
{
	int i, num1, num2;
	Ensembl_gene **genelist1, **genelist2;
	
	genelist1=new Ensembl_gene * [node1->get_num_genes()];
	genelist2=new Ensembl_gene * [node2->get_num_genes()];
	
	num1=node1->get_num_genes();
	num2=node2->get_num_genes();
	
	for(i=0; i<num1; i++)
		genelist1[i]=node1->get_gene(i);
	
	for(i=0; i<num2; i++)
		genelist2[i]=node2->get_gene(i);

	for(i=0; i<num1; i++) 
		node2->add_gene(genelist1[i]);
	
	for(i=0; i<num2; i++) 
		node1->add_gene(genelist2[i]);

	
	delete[] genelist1;
	delete[] genelist2;
}

void collapse_tandems(int num_genes, int &new_num_objs, Orthology_node *&new_list, Graph<Orthology_node> *the_graph)
{
	int i, j, k, l, cnt;
	BOOL do_merge;
	Edge<Orthology_node> *edge_obj;
	Node<Orthology_node> *node1, *node2, *edgenode;
	
	new_num_objs=num_genes;
	
	for(i=0; i<the_graph->get_num_nodes(); i++) {
		node1=the_graph->get_node(i);
		
		j=0;
		while(j<node1->get_num_edges()) {
			j=0;
			node1->reset_edges();
			do_merge=FALSE;
			while((do_merge==FALSE) && (j<node1->get_num_edges())) {
				//for(j=0; j<node1->get_num_edges(); j++) {
				node2=node1->get_next_edge();
			
				do_merge=FALSE;
				if ((node1->get_element()->omit_obj() == FALSE) && (node2->get_element()->omit_obj() == FALSE)) {
					for(k=0; k<node1->get_element()->get_num_genes(); k++) {
						for(l=0; l<node2->get_element()->get_num_genes(); l++) {
							if ((node1->get_element()->get_position() != -1) && (node1->get_element()->get_chromosome() == node2->get_element()->get_chromosome())) {
								if ((node1->get_element()->get_gene(k)->position == node2->get_element()->get_gene(l)->position+1) || 
									(node1->get_element()->get_gene(k)->position == node2->get_element()->get_gene(l)->position-1)) {
									do_merge=TRUE;
								}
							}
						}
					}
				
					if (do_merge == TRUE) {							
						merge_nodes(node1->get_element(), node2->get_element());
						
						new_num_objs--;
						node2->get_element()->set_omit();
						
						node2->reset_edges();
						for(k=0; k<node2->get_num_edges(); k++) {
							edge_obj=node2->get_next_edge_obj();
							edgenode=edge_obj->get_other_node(node2);
							if ((node1->edge_exists(edgenode) == FALSE ) && 
								(edgenode->get_element()->omit_obj() == FALSE) && (node1 != edgenode)) {
								node1->add_edge(edgenode);
								node1->get_edge(edgenode)->set_double_weight(edge_obj->get_double_weight());
							}
						}
					
					}
				}
				j++;
			}
		}
	}
	
	cnt=0;
	new_list = new Orthology_node[new_num_objs];
	for(i=0; i<the_graph->get_num_nodes(); i++) {
		if (the_graph->get_node(i)->get_element()->omit_obj() == FALSE) {
			new_list[cnt++]=(*the_graph->get_node(i)->get_element());
			for(j=0; j<the_graph->get_node(i)->get_element()->get_num_genes(); j++) {
				if (strcmp("ENSECAG00000012544", the_graph->get_node(i)->get_element()->get_gene(j)->name) == 0) 
					cout<<"After tandem: "<<the_graph->get_node(i)->get_element()->get_gene(j)->name<<endl;
			}
		}
	}

}


BOOL read_orthofile(char *filename, Graph<Orthology_node> *the_graph, Graph<Orthology_node> *single_gene_graph, BOOL use_ka)
{
	int i, j,k, l, m, len, cnt;
	char name1[50], name2[50], sat[10];
	double ks, ka;
	BOOL found1, found2;
	ifstream genefile;
	Node<Orthology_node> *node1, *node2, *snode1, *snode2;
    std::map<std::string, int> Node_hash, SNode_hash;
    
    Node_hash.clear();
    SNode_hash.clear();
    
    for(i=0; i<the_graph->get_num_nodes(); i++) {
        for(j=0; j<the_graph->get_node(i)->get_element()->get_num_genes(); j++) {
            Node_hash[the_graph->get_node(i)->get_element()->get_gene(j)->name]=i;
        }
    }
	
    for(i=0; i<single_gene_graph->get_num_nodes(); i++) {
        SNode_hash[single_gene_graph->get_node(i)->get_element()->get_gene(0)->name]=i;
    }
	genefile.open(filename);
	if (genefile.fail())
		return(FALSE);
	
	
	if (!genefile.eof()) {
		genefile>>name1>>name2>>ks>>ka>>len>>sat;
		
		while (!genefile.eof())
		{
			//i=0;
			found1=FALSE;
            
            if (!(Node_hash.find(name1) == Node_hash.end())) {
                found1=TRUE;
                node1=the_graph->get_node(Node_hash[name1]);
            }
			//while ( (i<the_graph->get_num_nodes()) && (found1 == FALSE)) {
			//	node1=the_graph->get_node(i);
			//	for(j=0; j<node1->get_element()->get_num_genes(); j++) {
			//		if (strcmp(node1->get_element()->get_gene(j)->name, name1) == 0) found1=TRUE;
			//	}
			//	if (found1 ==FALSE)
			//		i++;
			//}
			
			//i=0;
			found2=FALSE;
            
            if (!(Node_hash.find(name2) == Node_hash.end())) {
                found2=TRUE;
                node2=the_graph->get_node(Node_hash[name2]);
            }
            
			//while ( (i<the_graph->get_num_nodes()) && (found2 == FALSE)) {
			//	node2=the_graph->get_node(i);
			//	for(j=0; j<node2->get_element()->get_num_genes(); j++) {
			//		if (strcmp(node2->get_element()->get_gene(j)->name, name2) == 0) found2=TRUE;
			//	}
			//	if (found2 ==FALSE)
			//		i++;
			//}
			
			if (strcmp(sat, "YES") == 0) ks=-1;
			
			if ((found1 == TRUE) && (found2 == TRUE)) {			  
				node1->add_edge(node2);
                //cout<<"Adding edge "<<node1->get_element()->get_gene(0)->name<<" --> "<<node2->get_element()->get_gene(0)->name<<endl;
                if (use_ka == FALSE)
                    node1->get_edge(node2)->set_double_weight(ks);
                else
                    node1->get_edge(node2)->set_double_weight(ka);
			}
			else
				cerr<<"Error: homologs |"<<name1<<"| and |"<<name2<<"| " <<found1<<", "<<found2<<" not found in graph\n";
			
			
			//Save the appropriate single node matches to save Ks
			//i=0;
			found1=FALSE;
			
            if(!(SNode_hash.find(name1)==SNode_hash.end())) {
                found1=TRUE;
                snode1=single_gene_graph->get_node(SNode_hash[name1]);
            }
            
            //while ( (i<single_gene_graph->get_num_nodes()) && (found1 == FALSE)) {
			//	snode1=single_gene_graph->get_node(i);
			//	if (strcmp(snode1->get_element()->get_gene(0)->name, name1) == 0) found1=TRUE;
			//	else i++;
			//}
			//i=0;
			found2=FALSE;
            
            if(!(SNode_hash.find(name2)==SNode_hash.end())) {
                found2=TRUE;
                snode2=single_gene_graph->get_node(SNode_hash[name2]);
            }
            
			//while ( (i<single_gene_graph->get_num_nodes()) && (found2 == FALSE)) {
			//	snode2=single_gene_graph->get_node(i);
			//	if (strcmp(snode2->get_element()->get_gene(0)->name, name2) == 0) found2=TRUE;
			//	else i++;
			//}
			
			if ((found1 == TRUE) && (found2 == TRUE)) {			  
				snode1->add_edge(snode2);
                if (use_ka == FALSE)
                    snode1->get_edge(snode2)->set_double_weight(ks);
                else
                    snode1->get_edge(snode2)->set_double_weight(ka);
			}
			
			genefile>>name1>>name2>>ks>>ka>>len>>sat;
		}
	}
	genefile.close();
	
	//Set average Ks values
	cout<<"Setting edge weights\n"<<flush;
	for(i=0; i<the_graph->get_num_nodes(); i++) {
		node1=the_graph->get_node(i);
		node1->reset_edges();
		
		for (j=0; j<node1->get_num_edges(); j++) {
			node2=node1->get_next_edge();
			
			ks=0;
			cnt=0;
			//cout<<"Setting "<<node1->get_element()->get_gene(0)->name<<": "<<node2->get_element()->get_gene(0)->name<<endl<<flush;
			for(k=0; k<node1->get_element()->get_num_genes(); k++) {
				m=0;
				found1=FALSE;
				while ( (m<single_gene_graph->get_num_nodes()) && (found1 == FALSE)) {
					snode1=single_gene_graph->get_node(m);
					if (strcmp(snode1->get_element()->get_gene(0)->name, node1->get_element()->get_gene(k)->name) == 0) found1=TRUE; 
					else m++;
				}
				
				if (found1 == FALSE) { snode1=0; cerr<<"Error: can't find gene "<<node1->get_element()->get_gene(k)->name<<endl;}
				//else {cout<<"Found single node: "<<snode1->get_element()->get_gene(0)->name<<endl;}
				
				for(l=0; l<node2->get_element()->get_num_genes(); l++) {
					m=0;
					found1=FALSE;
					while ( (m<single_gene_graph->get_num_nodes()) && (found1 == FALSE)) {
						snode2=single_gene_graph->get_node(m);
						if (strcmp(snode2->get_element()->get_gene(0)->name, node2->get_element()->get_gene(l)->name) == 0) found1=TRUE; 
						else m++;
					}
					if (found1 == FALSE) { snode2=0; cerr<<"Error: can't find gene "<<node2->get_element()->get_gene(l)->name<<endl;}
					//else {cout<<"Found single node: "<<snode2->get_element()->get_gene(0)->name<<endl;}
					
					
					if ((snode1 == 0) || (snode2 == 0)) {
						cerr<<"ERROR can't find "<<node1->get_element()->get_gene(k)->name<<": "<<snode1<<" or "
						<<node2->get_element()->get_gene(l)->name<<": "<<snode2<<endl;
					}
					
					if (snode1->edge_exists(snode2) == TRUE) {
						if (snode1->get_edge(snode2)->get_double_weight() != -1) {
							ks += snode1->get_edge(snode2)->get_double_weight();
							cnt++;
						}
					}
					//else { cout<<"ERROR missing edge: "<<snode1->get_element()->get_gene(0)->name<<"\t"<<snode2->get_element()->get_gene(0)->name<<endl;}
				}
			}
			
			if (ks != 0) {
				node1->get_edge(node2)->set_double_weight(ks/(double)cnt);
				//cout<<node1->get_element()->get_gene(0)->name<<": "<<node2->get_element()->get_gene(0)->name<<": "<<node1->get_edge(node2)->get_double_weight()<<endl;
			}
		}
		
		
	}
	
	return(TRUE);
	
	
}


void find_ortholog_self_matches(Graph<Orthology_node> * the_graph, List<Ortholog_pair> *orthologs, double ks_thresh)
{
	int i;
	Node<Orthology_node> *node1, *node2;
	Ortholog_pair a_pair;
	
	for(i=0; i<the_graph->get_num_nodes(); i++) {
		node1=the_graph->get_node(i);
		node1->set_uncovered();
	}

	
	for(i=0; i<the_graph->get_num_nodes(); i++) {
		node1=the_graph->get_node(i);
		node1->set_covered();
		
		//one-to-one relations
		if (node1->get_num_edges() == 1) {
			node1->reset_edges();
			node2=node1->get_next_edge();
			if (node2->is_covered() == FALSE) {
				if (node2->get_num_edges() == 1) {
					if ((node1->get_element()->left_complete == FALSE) && (node1->get_element()->right_complete == FALSE)) {
						if ((node1->get_element()->in_ortho_list == FALSE) && (node2->get_element()->in_ortho_list == FALSE)){
							if ((node1->get_edge(node2)-> get_double_weight() <= ks_thresh) && (node1->get_edge(node2)-> get_double_weight() != -1) ) {
								a_pair.ortho_pair[0]=node1->get_element();
								a_pair.ortho_pair[1]=node2->get_element();
								node1->get_element()->in_ortho_list=TRUE;
								node2->get_element()->in_ortho_list=TRUE;
								node1->get_element()->set_ortholog(node2->get_element());
								node2->get_element()->set_ortholog(node1->get_element());
								orthologs->add_to_list(a_pair);
							}
						}
					}
				}
			}
		}
	}
	
	
	
}


//new function dealing with nodes that have only one edge that is below threshold
void find_more_anchors(Graph<Orthology_node> * the_graph, List<Ortholog_pair> *orthologs, double ks_thresh)
{
	int i, j, cnt;
	double min_weight, next_weight;
	Node<Orthology_node> *node1, *node2, *curr_min;
	Ortholog_pair a_pair;
	
	
	for(i=0; i<the_graph->get_num_nodes(); i++) {
		node1=the_graph->get_node(i);
		node1->set_uncovered();
		}
		
		
	for(i=0; i<the_graph->get_num_nodes(); i++) {
		node1=the_graph->get_node(i);
		node1->set_covered();
		if (node1->get_num_edges() > 1) {
			node1->reset_edges();
			node2=node1->get_next_edge();
			min_weight=node1->get_edge(node2)-> get_double_weight();
			curr_min=node2;
			next_weight=-1.0;
		
			for (j=1; j<node1->get_num_edges(); j++) {
				node2=node1->get_next_edge();
			
				if (node1->get_edge(node2)-> get_double_weight() < min_weight) {
					next_weight=min_weight;
					min_weight=node1->get_edge(node2)-> get_double_weight();
					curr_min=node2;
				}
				else {
					if ((next_weight == -1.0) || (node1->get_edge(node2)-> get_double_weight() < next_weight)) {
						next_weight=node1->get_edge(node2)-> get_double_weight();
					}
			
				}   
			}
		
			if ((min_weight < ks_thresh) && (next_weight >= 3.0*min_weight) && (min_weight >= 0.01 ) ) {
				node2=curr_min;
				if ( (node2->is_covered() == FALSE)  &&(node1->get_element()->in_ortho_list == FALSE) && (node2->get_element()->in_ortho_list==FALSE)) {
					a_pair.ortho_pair[0]=node1->get_element();
					a_pair.ortho_pair[1]=node2->get_element();
					node1->get_element()->in_ortho_list=TRUE;
					node2->get_element()->in_ortho_list=TRUE;
					node1->get_element()->set_ortholog(node2->get_element());
					node2->get_element()->set_ortholog(node1->get_element());
					orthologs->add_to_list(a_pair);
					}
			}
			}
			else {
				if (node1->get_num_edges() > 0) {
				node1->reset_edges();
				
				node2=node1->get_next_edge();
				if ((node1->get_edge(node2)-> get_double_weight() <= ks_thresh) && (node1->get_edge(node2)-> get_double_weight() != -1) ) {
					if ( (node2->is_covered() == FALSE)  &&(node1->get_element()->in_ortho_list == FALSE) && (node2->get_element()->in_ortho_list==FALSE)) {		
								a_pair.ortho_pair[0]=node1->get_element();
								a_pair.ortho_pair[1]=node2->get_element();
								node1->get_element()->in_ortho_list=TRUE;
								node2->get_element()->in_ortho_list=TRUE;
								node1->get_element()->set_ortholog(node2->get_element());
								node2->get_element()->set_ortholog(node1->get_element());
								orthologs->add_to_list(a_pair);
								}
							}
							}
			}
	}


	
}

void clear_ortho_edges(Graph<Orthology_node> *the_graph, List<Ortholog_pair> *orthologs)
{
	int i;
	Ortholog_pair *curr_pair;
	Node<Orthology_node> *curr_node, *temp;
	
	for(i=0; i<orthologs->get_list_length(); i++) {
		curr_pair=orthologs->get_nth_element(i)->item();
		
		curr_node=the_graph->get_node(curr_pair->ortho_pair[0]->get_node_num());
		
		curr_node->reset_edges();
		
		while (curr_node->get_num_edges() > 0) {
			temp=curr_node->get_next_edge();
			curr_node->delete_edge(temp);
			curr_node->reset_edges();
		}
		
		curr_node=the_graph->get_node(curr_pair->ortho_pair[1]->get_node_num());
		
		curr_node->reset_edges();
		
		while (curr_node->get_num_edges() > 0) {
			temp=curr_node->get_next_edge();
			curr_node->delete_edge(temp);
			curr_node->reset_edges();
		}
		
	}
	
}

void remove_orphans(Graph<Orthology_node> *the_graph, Orthology_node **full_list, int full_size, int &new_num_objs, Orthology_node **&new_list, int &num_orphans, Orthology_node **&orphan_list)
{
	int i, cnt_real, cnt_orphans;	
	
	
	num_orphans=0;
	new_num_objs=full_size;
	
	for(i=0; i<full_size; i++) {
		if((the_graph->get_node(full_list[i]->get_node_num())->get_num_edges() == 0) && (the_graph->get_node(full_list[i]->get_node_num())->get_element()->in_ortho_list==FALSE)) 
			{num_orphans++; new_num_objs--;}
	}
	
	orphan_list=new Orthology_node * [num_orphans];
	new_list=new Orthology_node * [new_num_objs];
	
	cnt_real=cnt_orphans=0;
	for(i=0; i<full_size; i++) {
		if((the_graph->get_node(full_list[i]->get_node_num())->get_num_edges() == 0) && (the_graph->get_node(full_list[i]->get_node_num())->get_element()->in_ortho_list==FALSE))  {
			orphan_list[cnt_orphans]=full_list[i];
			orphan_list[cnt_orphans]->orphan=TRUE;
			cnt_orphans++;
		}
		else
			new_list[cnt_real++]=full_list[i];
				
	}
	//cout<<"Found "<<cnt_real<<" genes\n";
}



void number_regions(Orthology_node **gene_list, int num_objs)
{
	int i, end;
	
	gene_list[0]->region_size=1;
	
	for(i=1; i<num_objs; i++) {
		if (gene_list[i]->left_complete==TRUE)
			gene_list[i]->region_size=gene_list[i-1]->region_size+1;
		else
			gene_list[i]->region_size=1;
	}
	

	end=num_objs-1;
	
	i=end;
	while(i > 0) {
		if (gene_list[i]->left_complete == TRUE) 
				gene_list[i-1]->region_size=gene_list[i]->region_size;
		i--;
			
	}
	
	i=0;

	for(i=0; i<num_objs; i++) {
	if (gene_list[i]->get_ortholog() != 0) 
		gene_list[i]->get_ortholog()->region_size=gene_list[i]->region_size;
	}

}



void rank_edges(Graph<Orthology_node> *the_graph, Graph<Orthology_node> *single_gene_graph)
{
	int i, j, k, l, m;
	double val;
	BOOL found1;
	Node<Orthology_node> *curr_node, *edge_node, *snode1, *snode2;
	
	
	for(i=0; i<the_graph->get_num_nodes(); i++) {
		curr_node=the_graph->get_node(i);
		curr_node->best_edge=0;
		curr_node->best_edge_score=0;
		
		if (curr_node->get_num_edges() > 0) {
			curr_node->reset_edges();
			edge_node=curr_node->get_next_edge();
			
			val=0;
			for(k=0; k<curr_node->get_element()->get_num_genes(); k++) {
				m=0;
				found1=FALSE;
				while ( (m<single_gene_graph->get_num_nodes()) && (found1 == FALSE)) {
					snode1=single_gene_graph->get_node(m);
					if (strcmp(snode1->get_element()->get_gene(0)->name, curr_node->get_element()->get_gene(k)->name) == 0) found1=TRUE; 
					else m++;
				}
				
				for(l=0; l<edge_node->get_element()->get_num_genes(); l++) {
					m=0;
					found1=FALSE;
					while ( (m<single_gene_graph->get_num_nodes()) && (found1 == FALSE)) {
						snode2=single_gene_graph->get_node(m);
						if (strcmp(snode2->get_element()->get_gene(0)->name, edge_node->get_element()->get_gene(l)->name) == 0) found1=TRUE; 
						else m++;
					}
					
					if (snode1->edge_exists(snode2) == TRUE) {
						if (snode1->get_edge(snode2)->get_double_weight() != -1) {
							if (snode1->get_edge(snode2)->get_double_weight() != 0) 
								val += 1.0/snode1->get_edge(snode2)->get_double_weight();
							else
								val=MAX_DOUBLE;
						}
					}
				}
			}
			
			curr_node->best_edge=curr_node->get_edge(edge_node);
			curr_node->best_edge_score=val;
			
			
			for(j=1; j<curr_node->get_num_edges(); j++) {
				edge_node=curr_node->get_next_edge();
				val=0;
				for(k=0; k<curr_node->get_element()->get_num_genes(); k++) {
					m=0;
					found1=FALSE;
					while ( (m<single_gene_graph->get_num_nodes()) && (found1 == FALSE)) {
						snode1=single_gene_graph->get_node(m);
						if (strcmp(snode1->get_element()->get_gene(0)->name, curr_node->get_element()->get_gene(k)->name) == 0) found1=TRUE; 
						else m++;
					}
					
					for(l=0; l<edge_node->get_element()->get_num_genes(); l++) {
						m=0;
						found1=FALSE;
						while ( (m<single_gene_graph->get_num_nodes()) && (found1 == FALSE)) {
							snode2=single_gene_graph->get_node(m);
							if (strcmp(snode2->get_element()->get_gene(0)->name, edge_node->get_element()->get_gene(l)->name) == 0) found1=TRUE; 
							else m++;
						}
						
						if (snode1->edge_exists(snode2) == TRUE) {
							if (snode1->get_edge(snode2)->get_double_weight() != -1) {
								if (snode1->get_edge(snode2)->get_double_weight() != 0) 
									val += 1.0/snode1->get_edge(snode2)->get_double_weight();
								else
									val=MAX_DOUBLE;
							}
						}
					}
				}
				
				if (val > curr_node->best_edge_score) {
					curr_node->best_edge=curr_node->get_edge(edge_node);
					curr_node->best_edge_score=val;
				}
			}
			
		}
	}
	
}

