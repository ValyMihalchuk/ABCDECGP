#include "pch.h"

Solution::Solution(const Abcde& _main_model, const Deep& _aux_model, const Parametrs& _param) {
	main_model = _main_model;
	aux_model = _aux_model;
	param = _param;
	alpha = 0.0;
}

inline void Solution::copy_posterior(Distribution::Posterior& posterior_to, Distribution::Posterior& posterior_from)
{
	for (int i = 0; i < main_model.count_iter; i++)
	{
		posterior_to.thetha[i] = posterior_from.thetha[i];
		posterior_to.w[i] = posterior_from.w[i];
		posterior_to.error[i] = posterior_from.error[i];
	}
}

void Solution::run_manager()
{
	manager.read_log_file(main_model.posterior, main_model.new_posterior, main_model.norm_error, main_model.count_iter, main_model.count_opt_param);
	switch (manager.state)
	{
	case Run_manager::STATE::INIT:
		run_init(manager.iter + 1, manager.index_thetha);
		break;
	case Run_manager::STATE::RUN_APPROXIMATE:
		run_approximate(manager.iter + 1, manager.index_thetha);
		break;
	case Run_manager::STATE::RUN:
		run(manager.iter + 1, manager.index_thetha);
		break;
	}
}

void Solution::run_init(int iter, int index_thetha)
{

	int size = 1;

		ofstream out("log_iteration.txt", std::ios::app);
		out << "INIT" << endl;
		vector<Distribution::Thetha> all_thetha;
		for (int j = 0; j < size; j++)
		{
			vector<vector<double>> _param;
			for (int i = 0; i < main_model.count_iter / size; i++)
			{
				main_model.curr_thetha = main_model.generate_vector_param(Distribution::NORM_WITH_PARAM);
				_param.push_back(main_model.curr_thetha.param);
				all_thetha.push_back(main_model.curr_thetha);
			}

		}
		for (int i = 0; i < main_model.count_iter / size; i++)
		{
			double error;
			main_model.curr_thetha = all_thetha[i];
			out << "iteration = " << -1 << endl;
			out << "element number = " << i << endl;
			for (int s = 0; s < main_model.count_opt_param; s++)
				out << main_model.curr_thetha.param[s] << endl;
			aux_model.create_tmp_deep_ini_file();
			int seed = main_model.generator.generate_seed();
			aux_model.prepare_tmp_deep_ini_file(main_model.bounds(main_model.curr_thetha), main_model.dtype, seed);
			error = aux_model.run(-1, i, seed);
			if (i == 0)
				main_model.norm_error = error;
			main_model.posterior.thetha[i] = main_model.curr_thetha;
			main_model.posterior.w[i] = 1.0 / main_model.count_iter;
			main_model.posterior.error[i] = error / main_model.norm_error;
			main_model.new_posterior.thetha[i] = main_model.curr_thetha;
			main_model.new_posterior.w[i] = 1.0 / main_model.count_iter;
			main_model.new_posterior.error[i] = error / main_model.norm_error;
			main_model.posterior.thetha[i].delta = main_model.new_posterior.thetha[i].delta = main_model.generator.prior_distribution(Distribution::TYPE_DISTR::EXPON, 0.005);
			out << "delta = " << main_model.posterior.thetha[i].delta << endl;
			out << "error = " << error / main_model.norm_error << endl;
		}

		for (int i = 0; i < main_model.count_iter; i++)
			manager.create_log_file(manager.state, main_model.posterior, main_model.new_posterior, main_model.norm_error, -1, i, main_model.count_opt_param);
		main_model.set_sample_dist_param();
		main_model.get_index_best();
		print_log(-1);
		manager.state = Run_manager::STATE::RUN_APPROXIMATE;
		manager.change_state(manager.state);
		out.close();
		run_approximate(0, 0);

}

void Solution::run_approximate(int iter, int index_thetha)
{

	int size = 1;

		ofstream out("log_iteration.txt", std::ios::app);
		out << "RUN_APPROXIMATE" << endl;
		for (int t = iter; t < main_model.start_iter; t++)
		{
			vector<Distribution::Thetha> all_thetha;
			for (int j = 0; j < size; j++)
			{
				vector<vector<double>> _param;

				for (int i = 0; i < main_model.count_iter / size; i++)
				{
					if (main_model.crossing_mode == Abcde::CROSSING_MODE::ALL)
					{
						double choice = main_model.generator.prior_distribution(Distribution::TYPE_DISTR::RANDOM, 0.0, 1.0);
						if (choice < 0.05)
						{
							main_model.curr_thetha = main_model.mutation(i + j * main_model.count_iter / size);
						}
						else
						{
							main_model.curr_thetha = main_model.crossover(i + j * main_model.count_iter / size);
						}
					}
					else if (main_model.crossing_mode == Abcde::CROSSING_MODE::ONLY_CROSSOVER)
						main_model.curr_thetha = main_model.crossover(i + j * main_model.count_iter / size);
					else if (main_model.crossing_mode == Abcde::CROSSING_MODE::ONLY_MUTATION)
						main_model.curr_thetha = main_model.mutation(i + j * main_model.count_iter / size);
					

					_param.push_back(main_model.curr_thetha.param);
					all_thetha.push_back(main_model.curr_thetha);
				}

			}
			for (int i = 0; i < main_model.count_iter / (size); i++)
			{
				double error;
				main_model.curr_thetha = all_thetha[i];
				out << "iteration = " << t << endl;
				out << "element number = " << i << endl;
				for (int s = 0; s < main_model.count_opt_param; s++)
					out << main_model.curr_thetha.param[s] << endl;
				out << "delta = " << main_model.curr_thetha.delta << endl;
				int seed = main_model.generator.generate_seed();
				aux_model.create_tmp_deep_ini_file();
				aux_model.prepare_tmp_deep_ini_file(main_model.bounds(main_model.curr_thetha), main_model.dtype, seed);
				error = aux_model.run(t, i, seed);
				out << "error = " << error / main_model.norm_error << endl;
				alpha = main_model.get_statistics(Parametrs::MODE::INIT, error / main_model.norm_error, i);
				out << "original alpha = " << alpha << endl;
				alpha = min(1.0, alpha);
				out << "alpha = " << alpha << endl;
				if (main_model.accept_alpha(alpha))
				{
					out << "accept alpha" << endl;
					main_model.new_posterior.thetha[i] = main_model.curr_thetha;
					main_model.new_posterior.error[i] = error / main_model.norm_error;
				}
				else
					out << "not accept" << endl;
			}
	
			main_model.get_index_best();

			main_model.update_posterior();//перерасчет весов
			copy_posterior(main_model.posterior, main_model.new_posterior); // перестановка
			main_model.set_sample_dist_param();

			for (int i = 0; i < main_model.count_iter; i++)
				manager.create_log_file(manager.state, main_model.posterior, main_model.new_posterior, main_model.norm_error, t, i, main_model.count_opt_param);
			print_log(t);
		}
		double s = 0.0;
		for (int i = 0; i < main_model.count_iter; i++)
		{
			s += main_model.posterior.thetha[i].delta;
		}
		main_model.posterior.delta_one = s / main_model.count_iter;
		main_model.new_posterior.delta_one = main_model.posterior.delta_one;
		manager.state = Run_manager::STATE::RUN;
		manager.change_state(manager.state);
		manager.change_delta(main_model.posterior.delta_one);
		out.close();
		run(0, 0);

}

void Solution::run(int iter, int index_thetha)
{

	int size = 1;

		ofstream out("log_iteration.txt", std::ios::app);
		out << "RUN" << endl;
		for (int t = iter; t < main_model.t; t++)
		{
			vector<Distribution::Thetha> all_thetha;
			for (int j = 0; j < size; j++)
			{
				vector<vector<double>> _param;

				for (int i = 0; i < main_model.count_iter / size; i++)
				{
					if (main_model.crossing_mode == Abcde::CROSSING_MODE::ALL)
					{
						double choice = main_model.generator.prior_distribution(Distribution::TYPE_DISTR::RANDOM, 0.0, 1.0);
						if (choice < 0.05)
						{
							main_model.curr_thetha = main_model.mutation(i + j * main_model.count_iter / size);
						}
						else
						{
							main_model.curr_thetha = main_model.crossover(i + j * main_model.count_iter / size);
						}
					}
					else if (main_model.crossing_mode == Abcde::CROSSING_MODE::ONLY_CROSSOVER)
						main_model.curr_thetha = main_model.crossover(i + j * main_model.count_iter / size);
					else if (main_model.crossing_mode == Abcde::CROSSING_MODE::ONLY_MUTATION)
						main_model.curr_thetha = main_model.mutation(i + j * main_model.count_iter / size);
					
					_param.push_back(main_model.curr_thetha.param);
					all_thetha.push_back(main_model.curr_thetha);
				}

			}
			for (int i = 0; i < main_model.count_iter / size; i++)
			{
				double error;
				main_model.curr_thetha = all_thetha[i];
				int seed = main_model.generator.generate_seed();
				aux_model.create_tmp_deep_ini_file();
				aux_model.prepare_tmp_deep_ini_file(main_model.bounds(main_model.curr_thetha), main_model.dtype, seed);
				error = aux_model.run(t, i, seed);
				out << "iteration = " << t << endl;
				out << "element number = " << i << endl;
				for (int s = 0; s < main_model.count_opt_param; s++)
					out << main_model.curr_thetha.param[s] << endl;
				out << "delta = " << main_model.curr_thetha.delta << endl;

				alpha = main_model.get_statistics(Parametrs::MODE::AUX, error / main_model.norm_error, i);
				out << "error = " << error / main_model.norm_error << endl;

				out << "original alpha = " << alpha << endl;
				alpha = min(1.0, alpha);
				out << "alpha = " << alpha << endl;
				if (main_model.accept_alpha(alpha))
				{
					out << "accept alpha" << endl;
					main_model.new_posterior.thetha[i] = main_model.curr_thetha;
					main_model.new_posterior.error[i] = error / main_model.norm_error;
				}
				else
					out << "not accept" << endl;
			}
			
			main_model.get_index_best();

			main_model.update_posterior();//перерасчет весов
			copy_posterior(main_model.posterior, main_model.new_posterior);//перестановка
			main_model.set_sample_dist_param();

			for (int i = 0; i < main_model.count_iter; i++)
				manager.create_log_file(manager.state, main_model.posterior, main_model.new_posterior, main_model.norm_error, t, i, main_model.count_opt_param);
			print_log(t);
		}
		manager.state = Run_manager::STATE::END;
		manager.change_state(manager.state);
		out.close();

}

void Solution::print_log(int iter)
{
	ofstream logfile("log_result.txt", std::ios::app);
	logfile << "iteration = " << iter << endl;
	for (int i = 0; i < main_model.count_iter; i++)
	{
		logfile << "element number = " << i << " ";
		for (int j = 0; j < main_model.count_opt_param; j++)
		{
			logfile << "param[" << j << "] = " << main_model.posterior.thetha[i].param[j] << " ";
		}
		logfile << "w = " << main_model.posterior.w[i] << " ";
		logfile << "error = " << main_model.posterior.error[i] << " ";
		logfile << "delta = " << main_model.posterior.thetha[i].delta << " ";
		logfile << endl;
	}
	logfile.close();
}
