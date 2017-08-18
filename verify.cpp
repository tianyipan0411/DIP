#include "option.h"
#include "hypergraph.hpp"
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>

#include <sstream>

using namespace std;


void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
int main(int argc, char ** argv)
{
	srand(time(NULL));
	
	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}

	char * inFile = op.getPara("-i");
	if (inFile == NULL){
		inFile = (char*)"network.bin";
	}
	
	char * outPath = op.getPara("-o");
	if (outPath == NULL){
		outPath = (char *)"Result/";
	}
	
	string pathString = string(outPath);
	string line;
    ifstream inputFile (inFile);
	double runningTime=0.0, speed_rate=0.0, speed_req=0.0,cov_req=0.0, speed_time=0.0;
	int seed_size=0, orig_seed_size=0;
	int intersection = 0;
	set<int> seeds;
	set<int> base_seeds;
	string graph_filename;
	double tot_sims = 10000;
	if (inputFile.is_open())
	{
		getline (inputFile,line);
		runningTime = stod(line);
		getline (inputFile,line);
		graph_filename = line;
		getline (inputFile,line);
		speed_rate = stod(line);
		//cout<<"speed:"<<speed_rate<<endl;
		getline (inputFile,line);
		speed_req = stod(line);
		//cout<<"speed_req:"<<speed_req<<endl;
		getline (inputFile,line);
		cov_req = stod(line);
		//cout<<"cov:"<<cov_req<<endl;
		getline (inputFile,line);
		speed_time = stod(line);
		getline (inputFile,line);
		vector<string> tmp_seeds = split(line,' ');
		for(unsigned int i=0; i<tmp_seeds.size();i++)
		{
			seeds.insert(stoi(tmp_seeds[i]));
		}
		getline (inputFile,line);
		seed_size = stoi(line);
		getline (inputFile,line);
		tmp_seeds.clear();
		tmp_seeds = split(line,' ');
		for(unsigned int i=0; i<tmp_seeds.size();i++)
		{
			int tmp_seed = stoi(tmp_seeds[i]);
			if(seeds.count(tmp_seed)) intersection++;
			base_seeds.insert(tmp_seed);
		}	
		getline (inputFile,line);
		orig_seed_size = stoi(line);
		inputFile.close();
	}
	string outFile = "Final/" + pathString + "/s_"+to_string(speed_rate)+"sr_"+to_string(speed_req)+"cr_"+to_string(cov_req)+".txt";
	std::ifstream ifile(outFile.c_str());
	
	if(!(bool)ifile)
	{
		cout<<outFile<<endl;
		ofstream out(outFile.c_str());
		out<<"Dataset RunningTime SpeedRate SpeedReq CovReq SeedSize SeedSizeBase Intersection ActCovNodes ActCovNodesBase SpeedTime ActSpeedTime ActSpeedTimeBase TotSuccess TotSuccessBase TotSims"<<endl;
		Graph g;
		g.readGraphIC(graph_filename.c_str());
		HyperGraph hg(g.getSize(), 1, 1);
		double tot_inf_nodes = 0.0;
		double tot_inf_base_nodes = 0.0;
		double tot_sped_time = 0.0;
		double tot_sped_time_base = 0.0;
		int success_count = 0;
		int success_count_base = 0;
		for (set<int>::iterator i = base_seeds.begin(); i != base_seeds.end(); i++) {
			int element = *i;
			//cout<<element<<endl;
			if(seeds.size()<base_seeds.size())
			{
				if (!seeds.count(element))
				{
					seeds.insert(element);
				}
			}
		}		

		for(int i = 1; i < tot_sims; i++)
		{
			int inf_nodes,inf_nodes_base;
			double sped_time,sped_time_base;
			if(hg.verifySeedSet(g, seeds, speed_rate, speed_req, cov_req, inf_nodes, sped_time))
			{
				success_count++;
			}
			tot_inf_nodes += inf_nodes;
			tot_sped_time += sped_time;
			if(hg.verifySeedSet(g, base_seeds, speed_rate, speed_req, cov_req, inf_nodes_base,sped_time_base))
			{
				success_count_base++;
			}
			tot_inf_base_nodes += inf_nodes_base;
			tot_sped_time_base += sped_time_base;
			if(i%1000==0)
			{
				cout<<outFile<<" "<<i<<endl;
			}
		}
		out<<pathString<<" "<<runningTime<<" "<<speed_rate<<" "<<speed_req<<" "<<cov_req<<" "<<seeds.size()<<" "<<base_seeds.size()<<" "<<intersection<<" "<<tot_inf_nodes/tot_sims<<" "<< tot_inf_base_nodes/tot_sims<<" "<<speed_time
			<<" "<<tot_sped_time/tot_sims<<" "<<tot_sped_time_base/tot_sims<<" "<<success_count<<" "<<success_count_base<<" "<<tot_sims<<endl;
	}	
	return 0;
}





