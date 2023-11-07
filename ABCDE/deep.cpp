#include "pch.h"

Deep::Deep() {}

Deep::Deep(const string& param)
{
	config_file = param;
	act_with_config_file();
}

Deep& Deep::operator=(const Deep& other)
{
	propTree = other.propTree;
	config_file = other.config_file;
	keys = other.keys;
	index_in_keys = other.index_in_keys;
	config_file = other.config_file;
	deep_exe = other.deep_exe;
	index_score = other.index_score;
	count_snp = other.count_snp;
	index_n = other.index_n;
	index_l = other.index_l;
	index_seed = other.index_seed;
	return *this;
}



void Deep::act_with_config_file()
{
    boost::property_tree::ini_parser::read_ini(config_file, propTree);
	deep_exe = propTree.get<std::string>("abcde.name_exe_file");
	index_score = stoi(propTree.get<std::string>("abcde.index_score"));
	count_snp = stoi(propTree.get<std::string>("abcde.count_snp"));
	index_n = stoi(propTree.get<std::string>("abcde.index_n"));
	index_l = stoi(propTree.get<std::string>("abcde.index_l"));
	index_seed = stoi(propTree.get<std::string>("abcde.index_seed"));
	vector<string> str_list;//, str_index;
	string s = propTree.get<std::string>("abcde.keys");
	boost::split(str_list, s, boost::is_any_of(";"));
	for (int i = 0; i < str_list.size(); i++)
	{
		keys.push_back(str_list[i]);
	}
	s = propTree.get<std::string>("abcde.index_in_keys");
	boost::split(str_list, s, boost::is_any_of(";"));
	for (int i = 0; i < str_list.size(); i++)
	{
		index_in_keys.push_back(stoi(str_list[i]));
	}
}

double Deep::run(int iter, int element_number, int seed)
{
	double res;
	namespace bp = boost::process;
	bp::ipstream is;
	string output;
	std::vector<std::string> data;
	std::string line;
	cout << bp::search_path(deep_exe).string() + " --default-name=" + tmp_config_file + " -d" << endl;
	bp::child c(bp::search_path(deep_exe).string() + " --default-name=" + tmp_config_file + " -d", bp::std_out > is);
	while (c.running() && std::getline(is, line) && !line.empty())
	{
		cout << line + " iteration = " << iter << " element number = " << element_number << " seed = " << seed << endl;
		data.push_back(line);
	}
	c.wait();
	output = data.back();
	res = parse_result(output);
	if(res >= 1.e+12)
	{
	    string error_name = "error " + to_string(iter) + " " + to_string(element_number) + " " + to_string(seed) + ".ini";
	    write_ini(error_name , propTree);
	}
	std::remove(tmp_config_file.c_str());
	return res;
}

double Deep::parse_result(string output)
{
	const char* pattern = ":[-+]?[0-9]*\\.?[0-9]+";
	boost::regex re(pattern);
	int i = 1;
	boost::sregex_iterator it(output.begin(), output.end(), re);
	boost::sregex_iterator end;
	for (; it != end; ++it)
	{
		if(i == index_score)
		{
			return stod(it->str().erase(0, 1));
		}
		i++;
	}
}

void Deep::prepare_tmp_deep_ini_file(Distribution::Thetha thetha, vector<int>& dtype, int seed)
{
	vector<string> delimeters = {" ", ";", "," };
	string str;
	vector<string> split_str;
	int index = 0;
	int add_int;
	string delimeter;
	for (auto& key : keys)
	{
		str = propTree.get<std::string>(key);
		for (auto& d : delimeters)
		{
			boost::split(split_str, str, boost::is_any_of(d));
			if (split_str.size() > 1)
			{
				delimeter = d;
				break;
			}
		}
		boost::split(split_str, str, boost::is_any_of(delimeter));
		if (dtype[index] == 0)
		{
			add_int = (int)thetha.param[index];
			split_str[index_in_keys[index]] = to_string(add_int);
		}
		else
		    split_str[index_in_keys[index]] = to_string(thetha.param[index]);
		
		index += 1;
		string output;
		for (int i = 0; i < split_str.size(); i++)
		{
			output += split_str[i];
			if(i < split_str.size()-1)
			    output += delimeter;
		}
		propTree.put(key, output);
	}
	//add seed
	string key = "default_model.command";
	str = propTree.get<std::string>(key);
	split_str.clear();
	boost::split(split_str, str, boost::is_any_of(" "));
	split_str[index_seed] = to_string(seed);		
	string output;
	for (int i = 0; i < split_str.size(); i++)
	{
		output += split_str[i];
		if (i < split_str.size() - 1)
			output += ' ';
	}
	propTree.put(key, output);
	//end add seed
	string command = propTree.get<std::string>(key);
	split_str.clear();
	boost::split(split_str, command, boost::is_any_of(" "));
	int n = stoi(split_str[index_n]);
	int l = stoi(split_str[index_l]);	
	key = "default_model.partsizes";
	string partsizes = propTree.get<std::string>(key);
	split_str.clear();
	boost::split(split_str, partsizes, boost::is_any_of(";"));
	split_str[0] = to_string(n * l);
	split_str[1] = to_string(n + n * count_snp);
	string output_partsizes;
	for (int i = 0; i < split_str.size(); i++)
	{
		output_partsizes += split_str[i];
		if (i < split_str.size() - 1)
			output_partsizes += ';';
	}
	propTree.put(key, output_partsizes);
	write_ini(tmp_config_file, propTree);
}
	
void Deep::create_tmp_deep_ini_file()
{
	const char* name = tmpnam(NULL);
	tmp_config_file = name;
}
