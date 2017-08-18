#include "option.h"
#include "hypergraph.hpp"
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
using namespace std;
bool CheckAnyTrue(vector<bool> speed_sample, vector<bool> cov_sample)
{
	for(unsigned int i = 0; i < speed_sample.size();i++){
		if(speed_sample[i]) return true;
	}
	for(unsigned int i = 0; i < cov_sample.size();i++){
		if(cov_sample[i]) return true;
	}
	return false;
}

void MIM (vector<Graph> &mtg, vector<HyperGraph> & mthg, int t, long long num, bool mo, vector<double> speed_thres, 
			vector<double> cov_thres, vector<int> speed_req, vector<int> cov_req, double time_threshold, double dec,
			vector<int> & seeds, unsigned int n, int k, long long int &totalDegree, vector<double> speed_ratio, vector<double>cov_ratio)
{
		// indicates if more samples are required for each speed/cov level
		vector<bool> speed_sample{};
		vector<bool> cov_sample{};
		// the coverage for each speed/cov level by the seed set
		vector<int> speed_deg{};
		vector<int> cov_deg{};
		for(unsigned int i = 0; i < speed_thres.size();i++){
			speed_sample.push_back(true);
			speed_deg.push_back(0);
		}
		for(unsigned int i = 0; i < cov_thres.size();i++){
			cov_sample.push_back(true);
			cov_deg.push_back(0);
		}
		long long int curSamples = num;
		while (CheckAnyTrue(speed_sample, cov_sample)){
			//cout<<"Samples "<<curSamples<<endl;
			seeds.clear();
			for(unsigned int i = 0; i < speed_thres.size();i++){
				speed_deg[i]=0;
			}
			for(unsigned int i = 0; i < cov_thres.size();i++){
				cov_deg[i]=0;
			}
			//cout<<"??????"<<endl;
			addHyperedge(mtg,mthg,t,curSamples,mo, speed_thres, cov_thres, speed_sample, cov_sample, time_threshold,dec);


			buildSeedSet(mthg[0],seeds,n,k,speed_ratio, cov_ratio, speed_deg, cov_deg);
			//cout<<"speed deg: "<<speed_deg[0]<<" speed req:"<<speed_req[0]<<endl;
			//cout<<"cov deg: "<<cov_deg[0]<<" cov req:"<<cov_req[0]<<endl;
			curSamples *= 2;
			for(unsigned int i = 0; i < speed_thres.size();i++){
				if(speed_deg[i] >= speed_req[i]/* || k == 1*/) speed_sample[i] = false;
			}
			for(unsigned int i = 0; i < cov_thres.size();i++){
				if(cov_deg[i] >= cov_req[i]/* || k == 1*/) cov_sample[i] = false;
			}
		}
		// calculate the estimation of f
		totalDegree = 0;
		//cout<<"sped covered: "<<speed_deg[0]<<"sped total: "<<mthg[0].getNumEdge_sped(0)<<endl;
		for(unsigned int i = 0; i < speed_thres.size();i++){
			totalDegree += ((double)speed_deg[i])/((double)mthg[0].getNumEdge_sped(i))*(double)n;
		}
		//cout<<"tot deg aft sped： "<<totalDegree<<endl;
		//cout<<"cov covered: "<<cov_deg[0]<<"cov total: "<<mthg[0].getNumEdge_cov(0)<<endl;
		for(unsigned int i = 0; i < cov_thres.size();i++){
			totalDegree += ((double)cov_deg[i])/((double)mthg[0].getNumEdge_cov(i))*(double)n;
		}
		//cout<<"tot deg aft cov： "<<totalDegree<<endl;
		
}

void MTAP(vector<Graph> &mtg, vector<HyperGraph> & mthg, double delta, double epsilon, int t, bool mo, vector<double> speed_thres, 
			vector<double> cov_thres, double time_threshold, double dec,
			vector<int> & seeds, unsigned int n, vector<double> speed_ratio, vector<double>cov_ratio)
{
	// the parameters
	double numSets = speed_ratio.size() + cov_ratio.size();
	if(numSets <= 1) numSets = 2;
	double e = exp(1);
	double sigma = sqrt(log(2*numSets/delta));
	double lognk, tau,phi,gamma;
	int k =0;
	double stopping_condition = 0;
	// requirement for each speed/cov level;
	vector<int> speed_req{};
	vector<int> cov_req{};
	for(unsigned int i = 0; i < speed_ratio.size();i++){
		stopping_condition += ((double) n) * speed_ratio[i];
	}
	for(unsigned int i = 0; i < cov_ratio.size();i++){
		stopping_condition += ((double) n) * cov_ratio[i];
	}
	
	long long int totalDegree = 0;
		
	long long int totalDegree_prev = 0;
	int k_prev = 0;

	while(totalDegree < (1-0.01)*stopping_condition - 1.1)
	{
		int inc = 1;
		if(totalDegree>0) inc = ((1-0.01)*stopping_condition - 1.1 - totalDegree) * (k-k_prev)/(totalDegree-totalDegree_prev);
		if(inc < 1)
		{
			inc = 1;
		}
		k_prev = k;
		k += inc;
		totalDegree_prev = totalDegree;
		totalDegree = 0;
		speed_req.clear();
		cov_req.clear();
		// update the parameters for a new k
		lognk = lgamma(n+1) - lgamma(k+1) - lgamma(n-k-1);
		tau = sqrt((1-1/e) * (lognk + log(2*numSets/delta)));
		phi = ((1-1/e)*sigma+tau)/epsilon;
		gamma = 2*(phi*phi + log((2*numSets*numSets)/(delta*(numSets - 1))));	
		for(unsigned int i = 0; i < speed_ratio.size();i++){
			if(speed_ratio[i] == 0) speed_req.push_back(0);
			else speed_req.push_back(gamma);
		}	
		for(unsigned int i = 0; i < cov_ratio.size();i++){
			if(cov_ratio[i] == 0) cov_req.push_back(0);
			else cov_req.push_back(gamma);
		}

		//cout<<"n: "<<n<<"k: "<<k<<" lognk: "<<lognk<<" tau: "<<tau<<" phi: "<<phi<<endl;
		//cout <<"gamma: "<<gamma<<endl;
		MIM (mtg, mthg, t, gamma, mo, speed_thres, cov_thres, speed_req, cov_req, time_threshold, dec,
				seeds, n, k, totalDegree, speed_ratio, cov_ratio);
		//cout<<"total: "<<totalDegree<<" req: "<<stopping_condition<<" k: "<<k<<" gamma: "<<gamma<<endl;
	}
}

void CalculateLipschitzFunctionValue(bool if_final, vector<int> &seeds, double time, double & speed_value, double & cov_value, vector<Graph> &mtg, vector<HyperGraph> & mthg, double delta, double epsilon, int t, bool mo, double time_threshold,  vector<double> delay_decs,
			 unsigned int n, vector<double> speed_ratio, vector<double> cov_ratio, vector<double> speedup_rates)
{
	// Calculated for each time t
	vector<double> speed_thres{};   
	vector<double> cov_thres{};
	vector<double> speedup_time_ints{};	
	vector<double> speedup_times {time};
	for(unsigned int i = 0; i < speedup_times.size();i++){
		if(i==0) speedup_time_ints.push_back(speedup_times[i]);
		else speedup_time_ints.push_back(speedup_times[i] - speedup_times[i-1]);
	}
	speedup_time_ints.push_back(time_threshold - speedup_times[speedup_times.size()-1]);
	for(unsigned int i = 0; i < speedup_times.size();i++){
		speed_thres.push_back(mthg[0].calc_modified_threshold(speedup_times[i], speedup_time_ints, delay_decs));
	}
	for(unsigned int i = 0; i < cov_ratio.size();i++){
		cov_thres.push_back(mthg[0].calc_modified_threshold(time_threshold, speedup_time_ints, delay_decs));
	}
	vector<double> ignore_speed_ratio(speed_ratio.size(),0);
	vector<double> ignore_cov_ratio(cov_ratio.size(),0);
	
	if(if_final)
	{
		//cout<<"??"<<endl;
		MTAP (mtg, mthg, delta,epsilon, t, mo, speed_thres, cov_thres, time_threshold, delay_decs[0],
				seeds, n, speed_ratio, cov_ratio);
		return;
	}
	if(time == 0)
	{
		speed_value = double(n) * speed_ratio[0];
	}
	else if(time == time_threshold)
	{
		speed_value = 0;
	}
	else
	{
		vector<int> loc_seeds;
		MTAP (mtg, mthg, delta,epsilon, t, mo, speed_thres, cov_thres, time_threshold, delay_decs[0],
				loc_seeds, n, speed_ratio, ignore_cov_ratio);
		speed_value = loc_seeds.size();
		mthg[0].clearEdges();
	}
	vector<int> loc_seeds;
	MTAP (mtg, mthg, delta,epsilon, t, mo, speed_thres, cov_thres, time_threshold, delay_decs[0],
				loc_seeds, n, ignore_speed_ratio, cov_ratio);
	//cout<<"!!"<<endl;
	if(time == time_threshold)
	{
		seeds.insert(seeds.begin(), loc_seeds.begin(),loc_seeds.end());
	}
	cov_value = loc_seeds.size();
	cout<<"time:"<<time<<"resultcov: "<<cov_value<<" resultspeed: "<<speed_value<<endl;
	mthg[0].clearEdges();
	
}

// Only consider one time speed up for now. (two speed ratio/delay_decs, one speedup_rates)
// return the selected time
double Lipschitz(vector<Graph> &mtg, vector<HyperGraph> & mthg, double delta, double epsilon, int t, bool mo, double time_threshold,  vector<double> delay_decs,
			vector<int> & seeds, vector<int> & base_seeds, unsigned int n, vector<double> speed_ratio, vector<double> cov_ratio, vector<double> speedup_rates)
{
	//Stores all values related to nodes selected by the algorithm
	vector<double> times;
	// decreasing with t
	vector<double> speed_values;
	// increasing with t
	vector<double> cov_values;
	// Lipschitz constant for interval (one less than the times vector)
	vector<double> ls;
	vector<double> Rs;
	vector<double>::iterator it;
	times.push_back(0);
	times.push_back(time_threshold);
	double speed_value,cov_value,l,R;
	//cout<<"??"<<endl;
	CalculateLipschitzFunctionValue(false, seeds, 0, speed_value, cov_value, mtg, mthg, delta, epsilon, t, mo, time_threshold, delay_decs,
			 n, speed_ratio, cov_ratio, speedup_rates);
	//cout<<".."<<endl;
	speed_values.push_back(speed_value);
	cov_values.push_back(cov_value);
	CalculateLipschitzFunctionValue(false, base_seeds,time_threshold, speed_value, cov_value, mtg, mthg, delta, epsilon, t, mo, time_threshold, delay_decs,
			 n, speed_ratio, cov_ratio, speedup_rates);
	speed_values.push_back(speed_value);
	cov_values.push_back(cov_value);
	l = max(speed_values[0]-speed_values[1], cov_values[1]-cov_values[0]);
	ls.push_back(l);
	R = (speed_values[0]+cov_values[0]+speed_values[1]+cov_values[1])/2-l*(times[1]-times[0])/2;
	Rs.push_back(R);
	int i = 1;

	while(times[i]-times[i-1]>1/ls[i-1])
	{
		//cout<<"!!"<<endl;
		double new_time = (times[i]+times[i-1])/2 + (speed_values[i-1]+cov_values[i-1]-speed_values[i]-cov_values[i])/2/ls[i-1];
		CalculateLipschitzFunctionValue(false, seeds, new_time, speed_value, cov_value, mtg, mthg, delta, epsilon, t, mo, time_threshold, delay_decs,
			 n, speed_ratio, cov_ratio, speedup_rates);
		it = times.begin();
		advance(it,i);
		times.insert(it,new_time);
		it = speed_values.begin();
		advance(it,i);
		speed_values.insert(it,speed_value);
		it = cov_values.begin();
		advance(it,i);
		cov_values.insert(it,cov_value);
		ls[i-1] = max(speed_values[i-1]-speed_values[i], cov_values[i]-cov_values[i-1]);
		Rs[i-1] = (speed_values[i-1]+cov_values[i-1]+speed_values[i]+cov_values[i])/2-l*(times[i]-times[i-1])/2;
		l = max(speed_values[i]-speed_values[i+1], cov_values[i+1]-cov_values[i]);
		R = (speed_values[i]+cov_values[i]+speed_values[i+1]+cov_values[i+1])/2-l*(times[i+1]-times[i])/2;
		it = ls.begin();
		advance(it,i);
		ls.insert(it,l);
		it = Rs.begin();
		advance(it,i);
		Rs.insert(it,R);
		int min_ind = 0;
		int min_val = 2*n;
		for(unsigned int j = 0; j < Rs.size(); j++)
		{
			if(Rs[j] < min_val)
			{
				min_val = Rs[j];
				min_ind = j;
			}
		}
		i = min_ind + 1;
		cout<<"time i: "<<times[i]<<" time i-1: "<<times[i-1]<<" lip: "<<ls[i-1]<<endl;
	}
	
	int min_ind = 0;
	int min_val = 2*n;
	for(unsigned int j = 0; j < times.size(); j++)
	{
		if(speed_values[j]+cov_values[j] < min_val)
		{
			min_val = speed_values[j]+cov_values[j];
			min_ind = j;
		}
	}
	cout<<"time: "<<times[min_ind]<<" value: "<<min_val<<endl;
	CalculateLipschitzFunctionValue(true, seeds, times[min_ind], speed_value, cov_value, mtg, mthg, delta, epsilon, t, mo, time_threshold, delay_decs,
			 n, speed_ratio, cov_ratio, speedup_rates);
	return times[min_ind]; 
}

// Only consider one time speed up for now. (two speed ratio/delay_decs, one speedup_rates)
// return the selected time
void MTAPCompOnly(vector<Graph> &mtg, vector<HyperGraph> & mthg, double delta, double epsilon, int t, bool mo, double time_threshold,  vector<double> delay_decs,
			vector<int> & seeds, vector<int> & base_seeds, unsigned int n, vector<double> speed_ratio, vector<double> cov_ratio, vector<double> speedup_rates)
{
	//Stores all values related to nodes selected by the algorithm
	vector<double> times;
	// decreasing with t
	vector<double> speed_values;
	// increasing with t
	vector<double> cov_values;
	// Lipschitz constant for interval (one less than the times vector)
	vector<double> ls;
	vector<double> Rs;
	vector<double>::iterator it;
	times.push_back(0);
	times.push_back(time_threshold);
	double speed_value,cov_value,l,R;

	CalculateLipschitzFunctionValue(false, base_seeds,time_threshold, speed_value, cov_value, mtg, mthg, delta, epsilon, t, mo, time_threshold, delay_decs,
			 n, speed_ratio, cov_ratio, speedup_rates);
	speed_values.push_back(speed_value);
	cov_values.push_back(cov_value);
	
}


int main(int argc, char ** argv)
{
	srand(time(NULL));
	
	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}
	//outerr.open("errlog.txt", std::ios_base::app);
	char * inFile = op.getPara("-i");
	if (inFile == NULL){
		inFile = (char*)"network.bin";
	}
	//cout<<inFile<<endl;
	char * model = op.getPara("-m");
	if (model == NULL)
			model = (char *) "LT";
	
	// number of threads
	int t = 1;
	char * tmp = op.getPara("-t");
	if (tmp != NULL){
		t = atoi(tmp);
	}

	vector<Graph> mtg;
	for(int ind = 0; ind < t; ++ind){
		Graph tmpG;
		if (strcmp(model, "LT") == 0){
			tmpG.readGraphLT(inFile);
		} else if (strcmp(model, "IC") == 0){
			tmpG.readGraphIC(inFile);
		} else {
			printf("Incorrect model option!");
			return -1;
		}
		mtg.push_back(tmpG);
	}

	int n = mtg[0].getSize();
	
	double speedup = 2;
	tmp = op.getPara("-s");
	if (tmp != NULL){
		speedup = 1 + 0.5*atof(tmp);
	}
	//cout<<"speedup"<<speedup<<endl;
	double speedreq = 100;

	double speedmulti = 100;
	double covmulti = 1;
	if(n > 20000 && n <200000)
	{
		speedmulti = 500;
		covmulti = 10;
	}
	else if(n>200000)
	{
		speedmulti = 5000;
		covmulti = 100;
	}
	tmp = op.getPara("-sr");
	if (tmp != NULL){
		speedreq = atof(tmp)*speedmulti;
	}
	double covreq = 0;
	tmp = op.getPara("-cr");
	if (tmp != NULL){
		covreq = (atof(tmp) * 200 + 800)*covmulti;
	}
	
	tmp = op.getPara("-epsilon");
	float epsilon = 0.1;
	if (tmp != NULL){
		epsilon = atof(tmp);
	}

	float delta = 1.0/n;
	//double k = 100;
	
	tmp = op.getPara("-delta");
	if (tmp != NULL){
		delta = atof(tmp);
	}

/* 	tmp = op.getPara("-k");
	if (tmp != NULL){
		k = atof(tmp);
	} */
	int mo = 0;
	if (strcmp(model, "IC") == 0){
			mo = 1;
	}
	char * outPath = op.getPara("-o");
	if (outPath == NULL){
		outPath = (char *)"Result/";
	}

	string pathString = string(outPath);
	string outFile = pathString + "/s_"+to_string(speedup)+"sr_"+to_string(speedreq)+"cr_"+to_string(covreq)+".result";
	std::ifstream ifile(outFile.c_str());
	
	
	if(!(bool)ifile)
	{
		cout <<outFile<<endl;
		//double iter=k;
		int speedup_time_count = 1;
		int cov_req_count = 1;
		vector<HyperGraph> mthg;
		for(int ind = 0; ind < t; ++ind){
			HyperGraph hg(n, speedup_time_count, cov_req_count);
			mthg.push_back(hg);
		}
		
		//int numNodes = mtg[0].getSize();
		double time_threshold = 10;
	/* 	vector<double> speedup_times {3,7};
		vector<double> speedup_rates {1, 1.2, 2};
		vector<double> delay_decs {0.5, 0.6, 1};
		vector<double> speedup_time_ints {3,4,3}; */
		
		//coverage req. 
		vector<double> speed_ratio {speedreq/(double)n};
		vector<double> cov_ratio{covreq/(double)n};
		//cout <<speed_ratio[0]<<endl;
		vector<double> speedup_rates {1};
		speedup_rates.push_back(speedup);
		vector<double> delay_decs{};
		
		for(unsigned int i = 0; i < speedup_rates.size();i++){
			delay_decs.push_back(speedup_rates[i]/speedup_rates[speedup_rates.size()-1]);
		}
		
		


		vector<int> seeds;
		vector<int> base_seeds;
		for(int a=0;a<10;a++)
		{
			clock_t start = clock();
			//cout<<"init time: "<<float(init - start)/CLOCKS_PER_SEC<<endl;
			mthg[0].clearAllEdges();
			//double selected_time = Lipschitz(mtg, mthg, delta, epsilon, t, mo, time_threshold, delay_decs, seeds, base_seeds, n, speed_ratio, cov_ratio, speedup_rates);
			MTAPCompOnly(mtg, mthg, delta, epsilon, t, mo, time_threshold, delay_decs, seeds, base_seeds, n, speed_ratio, cov_ratio, speedup_rates);
			cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << "s "<<endl;
		}
		/* ofstream out(outFile.c_str());
		out << (float)(clock()-start)/CLOCKS_PER_SEC <<endl;
		out << inFile << endl;
		out << speedup << endl;
		out << speedreq << endl;
		out << covreq << endl;
		
		out<<selected_time<<endl;
		for (unsigned int i = 0; i < seeds.size(); ++i){
			//cout << seeds[i] << " ";
			out << seeds[i] << " ";
		} 
		out << endl;
		out << seeds.size()<<endl;
		
		// seed set calculated without the speedup
		for (unsigned int i = 0; i < base_seeds.size(); ++i){
			//cout << seeds[i] << " ";
			out << base_seeds[i] << " ";
		} 
		out << endl;
		out << base_seeds.size()<<endl; */
		
		//cout << endl;
		
	/* 	cout << "Seed Nodes: ";

		for (unsigned int i = 0; i < seeds.size(); ++i){
			cout << seeds[i] << " ";
			//out << seeds[i] << endl;
		} */
		//cout << endl;
		//cout << "Influence: " << totalDegree << endl;
		//cout<<"Seed set size: "<<seeds.size()<<endl;
		//cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << "s scenario: " <<outFile<< endl;
		
	/* 	cout << "Seed Nodes: ";
		ofstream out(outFile);
		for (unsigned int i = 0; i < seeds.size(); ++i){
			cout << seeds[i] << " ";
			out << seeds[i] << endl;
		}
		out.close();
		cout << endl << endl;
		//printf("Benefit Gain: %0.2lf\n",totalDegree*g.getTotalBenefit()/(double)hg.getNumEdge());
		cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << "s" << endl;
		cout << "Memory: " << getCurrentMemoryUsage() << " MB" << endl; */
		//vector<NodeTime> visit_mark(numNodes+1,{0,0.0});
		//vector<bool> visit(numNodes+1,false);
		
	/* 	hg.pollingIC_base(g, visit, visit_mark, time_threshold, delay_decs[0]);
		
		vector<NodeTime> sample = hg.getBaseEdge(0);
		cout<<"Base Sample"<<endl;
		for (unsigned int i = 0; i < sample.size(); ++i){
			cout << sample[i].nodeID << " ";
			cout << sample[i].influenceTime << endl;
		} */
		
		//vector<int> visit_mark1(numNodes+1,0);

		
		
	/* 	cout<<"thres: "<< time_threshold <<endl;
		cout<<" new thres cov: "<< new_threshold_cov<<endl;
		cout<<" new thres sped1: "<< new_threshold_sped1<<endl;
		
		hg.pollingIC_coverage(visit_mark1,new_threshold_cov, 0);
		hg.pollingIC_speedup(visit_mark1, new_threshold_sped1, 0);

		vector<int> sample1 = hg.getEdge_cov(0,0);
		cout<<"sample for coverage"<<endl;
		for (unsigned int i = 0; i < sample1.size(); ++i){
			cout << sample1[i] << endl;
		}
		vector<int> sample2 = hg.getEdge_sped(0,0);
		cout<<"sample for speed up"<<endl;
		for (unsigned int i = 0; i < sample2.size(); ++i){
			cout << sample2[i] << endl;
		}
		
		cout<<"orig: "<<sample.size()<<endl;
		cout<<"cov: "<<sample1.size()<<endl;
		cout<<"sped1: "<<sample2.size()<<endl; */
	}
	return 0;
}

