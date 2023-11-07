#include "pch.h"


int main(int argc, char* argv[])
{
	
	int numTask, rank;


	Parametrs param;
	param.process_program_options(argc, argv);

	cout << "Abcde";
	cout << param.config_file;
	Abcde abcde(param.config_file);
	cout << "ABCDE OK";

	Deep deep(param.config_file);

	cout << "DEEP OK";
	Solution solution(abcde, deep, param);
	solution.run_manager();

	cout << "Solution Ok";

	return 0;
}
