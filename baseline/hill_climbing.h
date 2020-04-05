#include <vector>
#include "SparseMatrix.h"

using namespace std;
class HillClimbing{
public:
	int type;
	vector<vector<int>*>* node_to_strat;//(S_v)
	vector<vector<double>*>* prob_strat_to_node;//(->at(i)->at(j), stra i to node j)
	SparseMatrix* gf;
	HillClimbing(SparseMatrix* gf, vector<vector<pair<int, double>*>*>* prob_activate=NULL, int type=0){
		this->type = type;
		this->gf = gf;
		if(prob_activate == NULL){
			this->prob_strat_to_node = NULL;
			this->node_to_strat = NULL;
			return;
		}
		int d = prob_activate->size();
		int n = gf->getSize();
		prob_strat_to_node = new vector<vector<double>*>;
		for(int i = 0; i < d; i++){
			vector<double>* tmp = new vector<double>;
			tmp->resize(n, 1.0);
			prob_strat_to_node->push_back(tmp);
		}
		node_to_strat = new vector<vector<int>*>;
		for(int i = 0; i < n; i++){
			vector<int>* tmp = new vector<int>;
			node_to_strat->push_back(tmp);
		}
		for(int i = 0; i < prob_activate->size(); i++){
			for(int j = 0; j < prob_activate->at(i)->size(); j++){
				prob_strat_to_node->at(i)->at(prob_activate->at(i)->at(j)->first) = prob_activate->at(i)->at(j)->second;
				node_to_strat->at(prob_activate->at(i)->at(j)->first)->push_back(i);
			}
		}
	}
	double function(double x);
	pair<vector<double>*, double>* run(int k, int theta, double delta, vector<vector<int>*>* RR_sets=NULL);
	void iter_greedy(vector<vector<int>*>* RR_sets, vector<double>* s, vector<double>* value, double del, int n);
	void iter_greedy_2(vector<vector<int>*>* RR_sets, vector<double>* s, vector<double>* value, double del, int n, vector<vector<pair<int, int>*>*>* strat_to_node_and_RR);
	vector< vector<int>*>* generate_RR_sets(int theta);
	void _generate_RR_set(int start_node, vector<bool>* active, vector<int>* RR_set);
};

class CIMM{
public:
	HillClimbing* hill_climbing;
	SparseMatrix* gf;
	CIMM(SparseMatrix* gf, HillClimbing* hill_climbing){
		this->hill_climbing = hill_climbing;
		this->gf = gf;
	};
	vector<double>* run(int k, double delta, double eps, double l);
	vector< vector<int>*>* generate_RR_sets(int theta, vector<vector<int>*>* RR_sets=NULL);
	void _generate_RR_set(int start_node, vector<bool>* active, vector<int>* RR_set);
	double _lambda2(double eps, int k, int n, double l, double delta, int d);
	double _lambda(double eps, int k, double delta, int n, double l, int d);
};