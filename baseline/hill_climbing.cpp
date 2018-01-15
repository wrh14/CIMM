#include "hill_climbing.h"
#include "math.h"

using namespace std;

double HillClimbing::function(double x)
{
	return 2 * x - x * x;
}

pair<vector<double>*, double>* HillClimbing::run(int k, int theta, double delta, vector<vector<int>*>* RR_sets){
	pair<vector<double>*, double>* res = new pair<vector<double>*, double>;
	res->first = new vector<double>;
	res->second = 0.0;

	bool RR_sets_null = (RR_sets == NULL);
	if(RR_sets_null){
		RR_sets = generate_RR_sets(theta);
	}

	vector<double>* s = new vector<double>;
	s->resize(RR_sets->size(), 1.0);
	int n = gf->getSize();

	if(type == 0){
		res->first->resize(n, 0.0);
		for(int i = 0; i < double(k) / delta; i++){
			iter_greedy(RR_sets, s, res->first, delta, n);
		}
	}
	else if(type == 1){
		int d = prob_strat_to_node->size();
		res->first->resize(d, 0.0);
		vector<vector<pair<int, int>*>*>* strat_to_node_and_RR = new vector<vector<pair<int, int>*>*>;
		for(int i = 0; i < d; i++){
			strat_to_node_and_RR->push_back(new vector<pair<int, int>*>);
		}
		for(int i = 0; i < RR_sets->size(); i++){
			vector<int>* RR_set = RR_sets->at(i);
			for(int j = 0; j < RR_set->size(); j++){
				int node = RR_set->at(j);
				for(int l = 0; l < node_to_strat->at(node)->size(); l++){
					strat_to_node_and_RR->at(node_to_strat->at(node)->at(l))->push_back(new pair<int, int>(i, node));
				}
			}
		}
		for(int i = 0; i < double(k) / delta; i++){
			iter_greedy_2(RR_sets, s, res->first, delta, n, strat_to_node_and_RR);
		}

		for(int i = 0; i < d; i++){
			vector<pair<int, int>*>* tmp_vec_pair = strat_to_node_and_RR->at(i);
			for(int j = 0; j < tmp_vec_pair->size(); j++){
				delete tmp_vec_pair->at(j);
			}
			delete strat_to_node_and_RR->at(i);
		}
		delete strat_to_node_and_RR;
	}

	for(int i = 0; i < s->size(); i++){
		res->second += 1 - s->at(i);
	}

	res->second *= double(n) / theta;

	if(RR_sets_null){
		for(int i = 0; i < RR_sets->size(); i++){
			delete RR_sets->at(i);
		}
		delete RR_sets;
	}
	delete s;
	return res;
}

void HillClimbing::iter_greedy(vector<vector<int>*>* RR_sets, vector<double>* s, vector<double>* value, double del, int n){
	vector<double>* delta = new vector<double>;

	delta->resize(n, 0.0);

	for(int i = 0; i < RR_sets->size(); i++){
		vector<int>* RR_set = RR_sets->at(i);
		for(int j = 0; j < RR_set->size(); j++){
			int node = RR_set->at(j);
			delta->at(node) += s->at(i) * (function(value->at(node) + del) - function(value->at(node))) / (1 - function(value->at(node)));
		}
	}

	int opt;
	double max_delta = -1;
	for(int i = 0; i < delta->size(); i++){
		if(delta->at(i) > max_delta){
			max_delta = delta->at(i);
			opt = i;
		}
	}

	double factor = (1 - function(value->at(opt) + del)) / (1 - function(value->at(opt)));
	for(int i = 0; i < RR_sets->size(); i++){
		vector<int>* RR_set = RR_sets->at(i);
		for(int j = 0; j < RR_set->size(); j++){
			int node = RR_set->at(j);
			if(node == opt){
				s->at(i) = s->at(i) * factor;
				break;
			}
		}
	}

	value->at(opt) += del;
	delete delta;
}

void HillClimbing::iter_greedy_2(vector<vector<int>*>* RR_sets, vector<double>* s, vector<double>* value, double del, int n, vector<vector<pair<int, int>*>*>* strat_to_node_and_RR){
	int d = strat_to_node_and_RR->size();
	vector<double>* delta = new vector<double>;
	delta->resize(d, 0.0);
	for(int i = 0; i < d; i++){
		int label = -1;
		double cur_delta = 1;
		vector<pair<int, int>*>* R_i = strat_to_node_and_RR->at(i);
		for(int j = 0; j < R_i->size(); j++){
			pair<int, int>* RR_and_node = R_i->at(j);
			int l = RR_and_node->first;
			int node = RR_and_node->second;
			if(l != label && label != -1){
				delta->at(i) += (1- cur_delta) * s->at(label);
				label = l;
				cur_delta = s->at(l);
			}
			else if(label == -1){
				label = l;
				cur_delta = s->at(l);
			}
			cur_delta *= pow(prob_strat_to_node->at(i)->at(node), value->at(i) + del) / pow(prob_strat_to_node->at(i)->at(node), value->at(i));
		}
		if(label != -1){
			delta->at(i) += (1- cur_delta) * s->at(label);
		}
	}

	int opt;
	double max_delta = -1;
	for(int i = 0; i < delta->size(); i++){
		if(delta->at(i) > max_delta){
			max_delta = delta->at(i);
			opt = i;
		}
	}

	for(int i = 0; i < RR_sets->size(); i++){
		for(int j = 0; j < RR_sets->at(i)->size(); j++){
			int node = RR_sets->at(i)->at(j);
			for(int l = 0; l < node_to_strat->at(node)->size(); l++){
				if(node_to_strat->at(node)->at(l) == opt){
					s->at(i) *= pow(prob_strat_to_node->at(opt)->at(node), value->at(opt) + del) / pow(prob_strat_to_node->at(opt)->at(node), value->at(opt));
				}
			}
		}
	}

	value->at(opt) += del;
	delete delta;
}

void HillClimbing::_generate_RR_set(int start_node, vector<bool>* active, vector<int>* RR_set){
	active->at(start_node) = true;
	RR_set->push_back(start_node);
	for(set<Edge>::iterator iter = gf->getIterFColBegin(start_node); iter != gf->getIterFColEnd(start_node); iter++){
		int target = iter->src;
		if(active->at(target)){
			continue;
		}
		double prob = iter->value;
		bool if_activated = ((rand() / (RAND_MAX+1.0)) < prob);
		if(if_activated){
			_generate_RR_set(target, active, RR_set);
		}
	}
}

vector<vector<int>*>* HillClimbing::generate_RR_sets(int theta){
	vector<vector<int>*>* res = new vector<vector<int>*>;

	int n = gf->getSize();

	for (int it = 0; it < theta; it++){
        vector<bool>* active = new vector<bool>;
        active->resize(n, false);
		vector<int>* RR_set = new vector<int>;
		int start_node = rand()%n;
		_generate_RR_set(start_node, active, RR_set);
		delete active;
		res->push_back(RR_set);
	}
	return res;
}

double CIMM::_lambda(double eps, int k, double delta, int n, double l, int d){
	return (2 + 2.0 / 3.0 * eps) * (double(k) / delta * log(d) + l * log(n) + log(log2(n))) * double(n) / (eps * eps);
}

double CIMM::_lambda2(double eps, int k, int n, double l, double delta, int d){
	double alpha = sqrt(l * log(n) + log(2));
	double beta = sqrt((1 - 1.0 / exp(1)) * (double(k) / delta * log(d) + l * log(n) + log2(n)));
	double alpha_beta = ((1 - 1.0 / exp(1))) * alpha + beta;
	return 2 * n * (alpha_beta * alpha_beta) / (eps * eps / 2.0);
}

vector<double>* CIMM::run(int k, double delta, double eps, double l){
	vector<vector<int>*>* RR_sets = new vector<vector<int>*>;
	int n = gf->getSize();
	double LB = 1.0;
	eps = eps * sqrt(2);
	l = l + log(2) / log(n);
	int max_round = max(max(log2((double)n), 1.0) - 1.0, 1.0);

	int d;
	if(this->hill_climbing->type == 0){
		d = n;
	}
	else if(this->hill_climbing->type == 1){
		d = this->hill_climbing->prob_strat_to_node->size();
	}
	else{
		d = -1;
	}

	for(int i = 1; i <  max_round + 1; i++){
		cout << "round " << i << ": ";
		double y = double(n) / pow(2, i);
		double lambda = _lambda(eps, k, delta, n, l, d);
		double theta = lambda / y;
		cout << "theta=" << theta << endl;
		int j = RR_sets->size();
		while(j < theta){
			generate_RR_sets(1, RR_sets);
			j = RR_sets->size();
		}
		pair<vector<double>*, double>* res = hill_climbing->run(k, RR_sets->size(), delta, RR_sets);
		if(res->second >= (1 + eps) * y){
			LB = res->second / (1 + eps);
			break;
		}
	}

	cout << "Final Step: ";

	int theta = int(_lambda2(eps, k, n, l, delta, d) / LB);
	vector<vector<int>*>* RR_sets_2 = new vector<vector<int>*>;//Regenerate a RR set.
	int j = 0;
	while(j < theta){
		RR_sets_2 = generate_RR_sets(1, RR_sets_2);
		j = RR_sets_2->size();
	}

	cout << "RR_sets' size=" << RR_sets_2->size() << endl;

	pair<vector<double>*, double>* res = hill_climbing->run(k, RR_sets_2->size(), delta, RR_sets_2);
	for(int i = 0; i < RR_sets->size(); i++){
		delete RR_sets->at(i);
	}
	for(int i = 0; i < RR_sets_2->size(); i++){
		delete RR_sets_2->at(i);
	}
	delete RR_sets;
	delete RR_sets_2;

	return res->first;
}

void CIMM::_generate_RR_set(int start_node, vector<bool>* active, vector<int>* RR_set){
	active->at(start_node) = true;
	RR_set->push_back(start_node);
	for(set<Edge>::iterator iter = gf->getIterFColBegin(start_node); iter != gf->getIterFColEnd(start_node); iter++){
		int target = iter->src;
		if(active->at(target)){
			continue;
		}
		double prob = iter->value;
		bool if_activated = ((rand() / (RAND_MAX+1.0)) < prob);
		if(if_activated){
			_generate_RR_set(target, active, RR_set);
		}
	}
}

vector<vector<int>*>* CIMM::generate_RR_sets(int theta, vector<vector<int>*>* RR_sets){
	if(RR_sets == NULL){
		RR_sets = new vector<vector<int>*>;
	}

	int n = gf->getSize();

	for (int it = 0; it < theta; it++){
        vector<bool>* active = new vector<bool>;
        active->resize(n, false);
		vector<int>* RR_set = new vector<int>;
		int start_node = rand()%n;
		_generate_RR_set(start_node, active, RR_set);
		delete active;
		RR_sets->push_back(RR_set);
	}
	return RR_sets;
}
