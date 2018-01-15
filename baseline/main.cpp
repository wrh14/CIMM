#include "cim_baseline.h"
#include "hill_climbing.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "math.h"
#include <ctime>
#include "event_timer.h"
using namespace std;

SparseMatrix* inputNetwork(string filename); //This file should have the nodeNum in the first line, edgeNum in the second line,and edges with their weight. Input is an actual.ifSelf is a bool value to indicate whether this network has a self ring.
vector<vector<pair<int, double>*>*>* inputLocToNode(string dataset);
double spread(SparseMatrix* gf, vector<double>* value, int sample_num_1, int sample_num_2);
double spread_2(SparseMatrix* gf, vector<double>* value, int sample_num_1, int sample_num_2, vector<vector<pair<int, double>*>*>* prob_activate);
void alo_hill_climbing();
void algo_ud();
void algo_cd();
void algo_hd();
void algo_cimm();
void algo_cimm_2();
void algo_mc_greedy();
void algo_mc_greedy_2();
void algo_cimm_value();

int main(){
	cout << "1. algo_cimm; 2. algo_ud; 3. algo_cd; 4. algo_hd; 5. algo_cimm_2; 6. algo_mc_greedy; 7. algo_mc_greedy_2; 8. algo_cimm_value" << endl;
	int algo = 0;
	cin >> algo;
	if(algo == 1){
		algo_cimm();
	}
	else if(algo == 2){
		algo_ud();
	}
	else if(algo == 3){
		algo_cd();
	}
	else if(algo == 4){
		algo_hd();
	}
	else if(algo == 5){
		algo_cimm_2();
	}
	else if(algo == 6){
		algo_mc_greedy();
	}
	else if(algo == 7){
		algo_mc_greedy_2();
	}
	else if(algo == 8){
		algo_cimm_value();
	}
	return 0;
}

void algo_cimm_value(){
	cout << "dataset:" << endl;
	string dataset;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);
	double l = 1.0;
	HillClimbing* hill_climbing = new HillClimbing(graph);
	CIMM* cimm = new CIMM(graph, hill_climbing);

	srand(time(0));
	ofstream file, file_time;

	double eps = 0.5;

	file.open(save_address + "_cimm_value_eps=0.5", ios::out);
	double delta = 1;
	int k = 50;
	vector<double>* value = cimm->run(k, delta, eps, l);
	for(int i = 0; i < value->size(); i++){
		file << i << " " << value->at(i) << endl;
	}
	delete value;
	file.close();
}

void algo_cimm(){
	cout << "dataset:" << endl;
	string dataset;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);
	double l = 1.0;
	HillClimbing* hill_climbing = new HillClimbing(graph);
	CIMM* cimm = new CIMM(graph, hill_climbing);

	srand(time(0));
	ofstream file, file_time;

	double eps = 5;

	file.open(save_address + "_cimm_eps=5_new", ios::out);
	file_time.open(time_address + "_cimm_eps=5_new", ios::out);
	double delta = 0.1;
	for(int k = 5; k < 51; k=k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value = cimm->run(k, delta, eps, l);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		delete pctimer;
		cout << "k= " << k << " delta= " << delta << " eps= " << eps << " time= " << running_time << endl;
		file_time << "k= " << k << " delta= " << delta << " eps= " << eps << " time= " << running_time << endl;
		int sample_num_1 = 10000;
		int sample_num_2 = 1;
		double num_spread = spread(graph, value, sample_num_1, sample_num_2);
		delete value;
		file << "k= " << k << " delta= " << delta << " eps= " << eps << " spread= " << num_spread << endl;
		cout << "k= " << k << " delta= " << delta << " eps= " << eps << " spread= " << num_spread << endl;
	}
	file.close();
	file_time.close();

	eps = 10;

	file.open(save_address + "_cimm_eps=10_new", ios::out);
	file_time.open(time_address + "_cimm_eps=10_new", ios::out);
	for(int k = 5; k < 51; k=k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value = cimm->run(k, delta, eps, l);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		delete pctimer;
		cout << "k= " << k << " delta= " << delta << " eps= " << eps << " time= " << running_time << endl;
		file_time << "k= " << k << " delta= " << delta << " eps= " << eps << " time= " << running_time << endl;
		int sample_num_1 = 10000;
		int sample_num_2 = 1;
		double num_spread = spread(graph, value, sample_num_1, sample_num_2);
		delete value;
		file << "k= " << k << " delta= " << delta << " eps= " << eps << " spread= " << num_spread << endl;
		cout << "k= " << k << " delta= " << delta << " eps= " << eps << " spread= " << num_spread << endl;
	}
	file.close();
	file_time.close();

}

void algo_cimm_2(){
	cout << "dataset:" << endl;
	string dataset, dataset_2;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset_2 = "../data_2/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);
	//input data of location
	double delta = 1;
	double l = 1.0;

	vector<vector<pair<int, double>*>*>* prob_activate = inputLocToNode(dataset_2);

	HillClimbing* hill_climbing = new HillClimbing(graph, prob_activate, 1);

	CIMM* cimm = new CIMM(graph, hill_climbing);

	srand(time(0));
	ofstream file, file_time;

	double eps = 0.5;

	file.open(save_address + "_cimm_2_eps=0.5_new", ios::out);
	file_time.open(time_address + "_cimm_2_eps=0.5_new", ios::out);
	for(int k = 5; k < 51; k=k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value = cimm->run(k, delta, eps, l);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		delete pctimer;
		cout << dataset << " :k= " << k << " delta= " << delta << " eps= " << eps << " time= " << running_time << endl;
		file_time << "k= " << k << " delta= " << delta << " eps= " << eps << " time= " << running_time << endl;
		int sample_num_1 = 10000;
		int sample_num_2 = 1;
		double num_spread = spread_2(graph, value, sample_num_1, sample_num_2, prob_activate);
		delete value;
		file << "k= " << k << " delta= " << delta << " eps= " << eps << " spread= " << num_spread << endl;
		cout << dataset << " :k= " << k << " delta= " << delta << " eps= " << eps << " spread= " << num_spread << endl;
	}
	file.close();
	file_time.close();
	 eps = 1;

	file.open(save_address + "_cimm_2_eps=1_new", ios::out);
	file_time.open(time_address + "_cimm_2_eps=1_new", ios::out);
	for(int k = 5; k < 51; k=k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value = cimm->run(k, delta, eps, l);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		delete pctimer;
		cout << dataset << " :k= " << k << " delta= " << delta << " eps= " << eps << " time= " << running_time << endl;
		file_time << "k= " << k << " delta= " << delta << " eps= " << eps << " time= " << running_time << endl;
		int sample_num_1 = 10000;
		int sample_num_2 = 1;
		double num_spread = spread_2(graph, value, sample_num_1, sample_num_2, prob_activate);
		delete value;
		file << "k= " << k << " delta= " << delta << " eps= " << eps << " spread= " << num_spread << endl;
		cout << dataset << " :k= " << k << " delta= " << delta << " eps= " << eps << " spread= " << num_spread << endl;
	}
	file.close();
	file_time.close();
}

void algo_mc_greedy(){
	cout << "dataset:" << endl;
	string dataset;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);
	//input data of location

	OriginalGreedy* original_greedy = new OriginalGreedy(graph);

	srand(time(0));
	ofstream file, file_time;

	double delta = 0.1;

	file.open(save_address + "_mc", ios::out);
	file_time.open(time_address + "_mc", ios::out);

	double total_time = 0.0;

	vector<double>* value = new vector<double>;
	vector<int>* iter_number = new vector<int>;
	priority_queue<pair<double, int>, vector<pair<double, int>>, less<pair<double, int>>>* q = new priority_queue<pair<double, int>, vector<pair<double, int>>, less<pair<double, int>>>;
	double num_spread = 0;
	double num_spread_2 = 0;
	for(int i = 0; i < 10; i++){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		num_spread = num_spread_2;
		num_spread_2 = original_greedy->run(5, delta, i * 5 / delta, num_spread, value, iter_number, q);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		total_time += running_time;
		delete pctimer;
		cout << "k= " << (i + 1) * 5 << " time= " << total_time << endl;
		file_time << "k= " << (i + 1) * 5 << " time= " << total_time << endl;
		file << "k= " << (i + 1) * 5 << " spread= " << num_spread_2 << endl;
		cout << "k= " << (i + 1) * 5 << " spread= " << num_spread_2 << endl;
	}
	file.close();
	file_time.close();
}

void algo_mc_greedy_2(){
	cout << "dataset:" << endl;
	string dataset, dataset_2;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset_2 = "../data_2/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);
	//input data of location
	double delta = 1;

	vector<vector<pair<int, double>*>*>* prob_activate = inputLocToNode(dataset_2);

	OriginalGreedy* original_greedy = new OriginalGreedy(graph, prob_activate, 1);

	srand(time(0));
	ofstream file, file_time;

	file.open(save_address + "_mc_2", ios::out);
	file_time.open(time_address + "_mc_2", ios::out);
	double total_time = 0.0;

	vector<double>* value = new vector<double>;
	vector<int>* iter_number = new vector<int>;
	priority_queue<pair<double, int>, vector<pair<double, int>>, less<pair<double, int>>>* q = new priority_queue<pair<double, int>, vector<pair<double, int>>, less<pair<double, int>>>;
	double num_spread = 0;
	for(int i = 0; i < 10; i++){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		num_spread = original_greedy->run(5, delta, i * 5 / delta, num_spread, value, iter_number, q);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		total_time += running_time;
		delete pctimer;
		cout << "k= " << (i + 1) * 5 << " time= " << total_time << endl;
		file_time << "k= " << (i + 1) * 5 << " time= " << total_time << endl;
		file << "k= " << (i + 1) * 5 << " spread= " << num_spread << endl;
		cout << "k= " << (i + 1) * 5 << " spread= " << num_spread << endl;
	}
	file.close();
	file_time.close();
}

void algo_ud(){
	cout << "dataset:" << endl;
	string dataset;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);

	int n = graph->getSize();
	double c = 0.1;
	srand(time(0));
	UD* ud = new UD();

	ofstream file, file_time;
	file.open(save_address + "_ud_eps=5e-1", ios::out);
	file_time.open(time_address + "_ud_eps=5e-1", ios::out);
	CIMM* cimm = new CIMM(new HillClimbing());
	double eps = 0.5;
	for(int k = 5; k < 51; k = k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value_opt = cimm->run(graph, k, 1.0, eps, 1.0);

		int sample_num_1 = 1000;
		int sample_num_2 = 100;
		/*
		vector<double>* value_opt = cimm->run(k, 1.0, eps, 1.0);
		double opt = spread(graph, value_opt, sample_num_1, sample_num_2);
*/
		pctimer->SetTimeEvent("opt");

		//int theta = 2.0 * n * (1 - 1.0 / exp(1)) * (k * log(n) + log(n) + log(2)) / (opt * eps * eps);
		int theta = int(n * log(n));
		pair<double, vector<int>*>* res_ud = ud->run(graph, k, c, theta);
		vector<double>* value = new vector<double>;
		value->resize(graph->getSize(), 0.0);

		int k_rec = 0;
		for(int i = 0; i < res_ud->second->size(); i++){
			k_rec += res_ud->first;
			value->at(res_ud->second->at(i)) = res_ud->first;
		}

		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		double opt_running_time = pctimer->TimeSpan("start", "opt");
		delete pctimer;
		cout << "k= " << k << " c= " << c << " eps= " << eps << " time= " << running_time << " opt_running_time= " << opt_running_time << endl;
		file_time << "k= " << k << " c= " << c << " eps= " << eps << " time= " << running_time << " opt_running_time= " << opt_running_time << endl;

		double num_spread = spread(graph, value, sample_num_1, sample_num_2);
		delete value;
		file << "k=" << k << " c= " << c << " eps=" << eps << " spread = " << num_spread << endl;
		cout << "k=" << k << " c= " << c << " eps=" << eps << " spread = " << num_spread << endl;
	}

	file.close();

	file.open(save_address + "_ud_eps=1e-1", ios::out);
	file.open(time_address + "_ud_eps=1e-1", ios::out);
	eps = 0.1;
	for(int k = 5; k < 51; k=k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value_opt = cimm->run(graph, k, 1.0, eps, 1.0);
		int sample_num_1 = 1000;
		int sample_num_2 = 100;
		double opt = spread(graph, value_opt, sample_num_1, sample_num_2);

		pctimer->SetTimeEvent("opt");


		int theta = 2.0 * n * (1 - 1.0 / exp(1)) * (k * log(n) + log(n) + log(2)) / (opt * eps * eps);
		pair<double, vector<int>*>* res_ud = ud->run(graph, k, c, theta);
		cout << res_ud->second->size() << endl;
		vector<double>* value = new vector<double>;
		value->resize(graph->getSize(), 0.0);

		int k_rec = 0;
		for(int i = 0; i < res_ud->second->size(); i++){
			k_rec += res_ud->first;
			value->at(res_ud->second->at(i)) = res_ud->first;
		}

		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		double opt_running_time = pctimer->TimeSpan("start", "opt");
		delete pctimer;
		cout << "k= " << k << " c= " << c << " eps= " << eps << " time= " << running_time << "opt_running_time= " << opt_running_time << endl;
		file_time << "k= " << k << " c= " << c << " eps= " << eps << " time= " << running_time << "opt_running_time= " << opt_running_time << endl;

		sample_num_1 = 1000;
		sample_num_2 = 100;
		double num_spread = spread(graph, value, sample_num_1, sample_num_2);
		delete value;
		file << "k=" << k << " c= " << c << " eps=" << eps << " spread = " << num_spread << endl;
		cout << "k=" << k << " c= " << c << " eps=" << eps << " spread = " << num_spread << endl;
	}
	file.close();
}

void algo_cd(){
	cout << "dataset:" << endl;
	string dataset;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);

	int n = graph->getSize();
	double c = 0.1;
	srand(time(0));
	CD* cd = new CD();
	CIMM* cimm = new CIMM(graph, new HillClimbing(graph));

	ofstream file, file_time;
	file.open(save_address + "_cd_eps=5e-1", ios::out);
	file_time.open(time_address + "_cd_eps=5e-1", ios::out);

	double eps = 0.5;
	for(int k = 5; k < 51; k = k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		int sample_num_1 = 1000;
		int sample_num_2 = 100;
		/*vector<double>* value_opt = cimm->run(k, 1.0, eps, 1.0);
		double opt = spread(graph, value_opt, sample_num_1, sample_num_2);*/

		pctimer->SetTimeEvent("opt");

		//int theta = 2.0 * n * (1 - 1.0 / exp(1)) * (k * log(n) + log(n) + log(2)) / (opt * eps * eps);
		int theta = int(n * log(n));
		vector<double>* value = cd->run(graph, k, c, theta);

		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		double opt_running_time = pctimer->TimeSpan("start", "opt");
		delete pctimer;
		cout << "k= " << k << " c= " << c << " eps= " << eps << " time= " << running_time << " opt_running_time= " << opt_running_time << endl;
		file_time << "k= " << k << " c= " << c << " eps= " << eps << " time= " << running_time << " opt_running_time= " << opt_running_time << endl;

		double num_spread = spread(graph, value, sample_num_1, sample_num_2);
		delete value;
		//delete value_opt;
		file << " k=" << k << " c= " << c << " eps=" << eps << " spread = " << num_spread << endl;
		cout << " k=" << k << " c= " << c << " eps=" << eps << " spread = " << num_spread << endl;
	}

	file.close();
}

void algo_hd(){
	cout << "dataset:" << endl;
	string dataset;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);
	double unit_delta = 0.1;
	double l = 1.0;
	HeuristicsDegree* heuristicsDegree = new HeuristicsDegree();

	srand(time(0));
	ofstream file, file_time;

	int M = 100;

	file.open(save_address + "_hd_M=100", ios::out);
	file_time.open(time_address + "_hd_M=100", ios::out);
	for(int k = 5; k < 51; k=k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value = heuristicsDegree->run(graph, k, M);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		delete pctimer;
		cout << "k= " << k << " M=" << M << " time= " << running_time << endl;
		file_time << "k= " << k << " M=" << M << " time= " << running_time << endl;
		int sample_num_1 = 1000;
		int sample_num_2 = 100;
		double num_spread = spread(graph, value, sample_num_1, sample_num_2);
		delete value;
		file << "k= " << k << " M=" << M << " spread= " << num_spread << endl;
		cout << "k= " << k << " M=" << M << " spread= " << num_spread << endl;
	}
	file.close();
	file_time.close();

	M = 200;
	file.open(save_address + "_hd_M=200", ios::out);
	file_time.open(time_address + "_hd_M=200", ios::out);
	for(int k = 5; k < 51; k=k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value = heuristicsDegree->run(graph, k, M);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		delete pctimer;
		cout << "k= " << k << " M=" << M << " time= " << running_time << endl;
		file_time << "k= " << k << " M=" << M << " time= " << running_time << endl;
		int sample_num_1 = 1000;
		int sample_num_2 = 100;
		double num_spread = spread(graph, value, sample_num_1, sample_num_2);
		delete value;
		file << "k= " << k << " M=" << M << " spread= " << num_spread << endl;
		cout << "k= " << k << " M=" << M << " spread= " << num_spread << endl;
	}
	file.close();
	file_time.close();
}


SparseMatrix* inputNetwork(string filename){
	SparseMatrix* result;

    int nodeNum;
	int edgeNum;

	ifstream infile(filename.c_str());
	cout<<"Begin to read the file: "<<filename<<endl;

	string s,s2;
	infile >> nodeNum >> edgeNum;
	getline(infile,s);
	cout<<"NodeNum is "<<nodeNum<<endl;
	cout << "EdgeNum is "<< edgeNum<<endl;
	result = new SparseMatrix(nodeNum);

	int tempNode1, tempNode2;
	double temp_1;

	for(int i=0; i<edgeNum * 2; i++){
    	getline(infile,s);
    	istringstream line (s);
    	line >> tempNode1 >> tempNode2 >> temp_1;
    	result->insert(Edge(tempNode1 - 1, tempNode2 - 1, temp_1));
	}

    cout << "Finish the reading! Real edgeNum is" << result->getNonZeroNum()<<endl;
	return result;
}

vector<vector<pair<int, double>*>*>* inputLocToNode(string filename){
	vector<vector<pair<int, double>*>*>* prob_activate = new vector<vector<pair<int, double>*>*>;
	ifstream infile(filename.c_str());
	int strat_num;
	string s,s2;
	infile >> strat_num;
	getline(infile, s);

	for(int i = 0; i < strat_num; i++){
		prob_activate->push_back(new vector<pair<int, double>*>);
	}

	getline(infile, s);
	while(infile){
		int strat_id;
		int node_id;
		double prob;
		istringstream line (s);
		line >> strat_id >> node_id >> prob;
		prob_activate->at(strat_id)->push_back(new pair<int, double>(node_id, prob));
		getline(infile, s);
	}

	for(int i = 0; i < strat_num; i++){
		prob_activate->push_back(new vector<pair<int, double>*>);
	}
	return prob_activate;
}

double _spread(SparseMatrix* gf, int start_node, vector<bool>* active){
	double res = 1;
	active->at(start_node) = true;
	for(set<Edge>::iterator iter = gf->getIterFRowBegin(start_node); iter != gf->getIterFRowEnd(start_node); iter++){
		int target = iter->dst;
		if(active->at(target)) continue;
		double prob = iter->value;
		bool if_activated = ((rand() / (RAND_MAX+1.0)) < prob);
		if(if_activated){
			res += _spread(gf, target, active);
		}
	}
	return res;
}

double function_v(double x){
	return 2 * x - x * x;
}

double spread(SparseMatrix* gf, vector<double>* value, int sample_num_1, int sample_num_2){
	double res = 0;
	for(int sample = 0; sample < sample_num_1; sample++){
		vector<int>* seed = new vector<int>;
		for(int i = 0; i < value->size(); i++){
			if((rand() / (RAND_MAX+1.0)) < function_v(value->at(i))){
				seed->push_back(i);
			}
		}
		for(int sample_2 = 0; sample_2 < sample_num_2; sample_2++){
			vector<bool>* active = new vector<bool>;
			active->resize(gf->getSize(), false);
			for(int i = 0; i < seed->size(); i++){
				if(active->at(seed->at(i))) continue;
				res += _spread(gf, seed->at(i), active);
			}
			delete active;
		}
		delete seed;
	}

	res = double(res) / (sample_num_1 * sample_num_2);
	return res;
}

/*
int main(){
	string s = "../data_example/dm-wc.txt";
	SparseMatrix* graph = inputNetwork(s);

	double k = 0;

	string s_2 = "./rr_cimm.txt";
	vector<double>* value = new vector<double>;
	value->resize(graph->getSize(), 0.0);
	ifstream infile(s_2.c_str());
	int nodeNum;
	string s_3, s_4;
	infile >> nodeNum;
	getline(infile,s_3);
	cout<<"NodeNum is "<<nodeNum<<endl;
	for(int i=0; i<nodeNum; i++){
    	getline(infile,s_3);
    	double tempValue = 0;
    	istringstream line (s_3);
    	line >> s_4 >> tempValue;
    	k = k + tempValue;
    	value->at(i) = tempValue;
	}

	cout << "k = " << k << endl;

	srand(time(0));

	int sample_num_1 = 100;
	int sample_num_2 = 1000;
	cout << "spread = " << spread(graph, value, sample_num_1, sample_num_2) << endl;
}
*/
/*
int main(){
	string s = "../data_example/dm-wc.txt";
	SparseMatrix* graph = inputNetwork(s);

	double k = 0;

	string s_2 = "./greedy.txt";
	vector<double>* value = new vector<double>;
	value->resize(graph->getSize(), 0.0);
	ifstream infile(s_2.c_str());
	int nodeNum;
	string s_3;
	infile >> nodeNum;
	getline(infile,s_3);
	cout<<"NodeNum is "<<nodeNum<<endl;
	for(int i=0; i<nodeNum; i++){
    	getline(infile,s_3);
    	double tempValue = 0;
    	int node;
    	istringstream line (s_3);
    	line >> node >> tempValue;
    	k = k + tempValue;
    	value->at(node - 1) = 1;
	}

	cout << "k = " << k << endl;

	srand(time(0));

	int sample_num_1 = 100;
	int sample_num_2 = 1000;
	cout << "spread = " << spread(graph, value, sample_num_1, sample_num_2) << endl;
}*/
void algo_hd(){
	cout << "dataset:" << endl;
	string dataset;
	cin >> dataset;
	string save_address = "../result/" + dataset;
	string time_address = "../time/" + dataset;
	dataset = "../data/" + dataset;
	SparseMatrix* graph = inputNetwork(dataset);
	double unit_delta = 0.1;
	double l = 1.0;
	HeuristicsDegree* heuristicsDegree = new HeuristicsDegree();

	srand(time(0));
	ofstream file, file_time;

	int M = 100;

	file.open(save_address + "_hd_M=100", ios::out);
	file_time.open(time_address + "_hd_M=100", ios::out);
	for(int k = 5; k < 51; k=k+5){
		EventTimer* pctimer = new EventTimer();
		pctimer->SetTimeEvent("start");
		vector<double>* value = heuristicsDegree->run(graph, k, M);
		pctimer->SetTimeEvent("end");
		double running_time = pctimer->TimeSpan("start", "end");
		delete pctimer;
		cout << "k= " << k << " M=" << M << " time= " << running_time << endl;
		file_time << "k= " << k << " M=" << M << " time= " << running_time << endl;
		int sample_num_1 = 1000;
		int sample_num_2 = 100;
		double num_spread = spread(graph, value, sample_num_1, sample_num_2);
		delete value;
		file << "k= " << k << " M=" << M << " spread= " << num_spread << endl;
		cout << "k= " << k << " M=" << M << " spread= " << num_spread << endl;
double spread_2(SparseMatrix* gf, vector<double>* value, int sample_num_1, int sample_num_2, vector<vector<double>*>* prob_activate){
	//TODO
	double res = 0;
	vector<double>* prob = new vector<double>;
	prob->resize(gf->getSize(), 1.0);
	for(int i = 0; i < prob_activate->size(); i++){
		vector<pair<int, double>*>* tmp_vec_pair = prob_activate->at(i);
		for(int j = 0; j < tmp_vec_pair->size(); j++){
			prob->at(tmp_vec_pair->at(j)->first) *= pow(tmp_vec_pair->at(j)->second, value->at(i));
		}
	}
	for(int i = 0; i < gf->getSize(); i++){
		prob->at(i) = 1 - prob->at(i);
	}
	for(int sample = 0; sample < sample_num_1; sample++){
		vector<int>* seed = new vector<int>;
		for(int i = 0; i < prob->size(); i++){
			if((rand() / (RAND_MAX+1.0)) < prob->at(i)){
				seed->push_back(i);
			}
		}
		for(int sample_2 = 0; sample_2 < sample_num_2; sample_2++){
			vector<bool>* active = new vector<bool>;
			active->resize(gf->getSize(), false);
			for(int i = 0; i < seed->size(); i++){
				if(active->at(seed->at(i))) continue;
				res += _spread(gf, seed->at(i), active);
			}
			delete active;
		}
		delete seed;
	}

	res = double(res) / (sample_num_1 * sample_num_2);
	return res;
}
