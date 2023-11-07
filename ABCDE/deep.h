#pragma once
#include "pch.h"


class Deep :public Model
{
private:
	pt::ptree propTree;
public:
	Deep();
	Deep(const string& param);
	Deep& operator=(const Deep&);
	void act_with_config_file();
	double run(int iter, int element_number, int seed);
	double parse_result(string output);
	void prepare_tmp_deep_ini_file(Distribution::Thetha thetha, vector<int>& dtype, int seed);
	void create_tmp_deep_ini_file();
	string config_file;
	string tmp_config_file;
	string deep_exe;
	double error;
	vector<string> keys;
	vector<int> index_in_keys;
	int count_snp;
	int index_n, index_l, index_score, index_seed;
};
