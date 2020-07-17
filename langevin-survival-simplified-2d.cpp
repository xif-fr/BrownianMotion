#include <fmt/core.h>
#include "pysimul-common.h"
#include <random>
#include <vector>
#include <algorithm>

// This is the boiled-down version of langevin-survival.cpp with the following flags :
// LANGEVIN_OVERDAMPED, TARGET_2D_CYL, FPT_INTERVAL, ENABLE_POISSON_RESET

constexpr size_t N_targets = 1;
const size_t pysimul_N = N_targets;

void* comp_thread (void* _data) {
	simul_thread_info_t& _thread = *(simul_thread_info_t*)_data;
	
	std::random_device _rd;
	std::mt19937 rng (_rd());
	std::normal_distribution<> normal_distrib_gen (0, 1); // mean, std
	std::uniform_real_distribution<> unif01 (0, 1);

	double t = 0;
	size_t n_trajectories = 0;
	_register_var(_thread, "n_trajectories", &n_trajectories);
	constexpr double Δt = 1e-8;
	_register_const(_thread, "Delta_t", Δt);
	uint8_t pause = 0;
	_register_var(_thread, "pause", &pause);

	// Diffusion
	double D;
	_register_var(_thread, "D", &D);

	// Gaussian distribution of initial position
	double init_pos_sigma;
	_register_var(_thread, "x0sigma", &init_pos_sigma);
	
	// Poissonian resetting
	double reset_rate;
	_register_var(_thread, "reset_rate", &reset_rate);
	const double proba_reset_step = Δt * reset_rate;

	// First passage time at the target @ x=survdist_time_pos
	std::vector<double> first_time;
	constexpr std::array<double, N_targets> first_times_xtarg = { 0.2 };
	_register_var(_thread, "first_times-0", &first_time);
	auto _first_times_xtarg = first_times_xtarg;
	_register_Narray(_thread, "first_times_xtarg", _first_times_xtarg);

	// 2D circular target with "tolerence radius" Rtol around x=x_targ,y=0
	double Rtol;
	_register_var(_thread, "2D-Rtol", &Rtol);
	auto check_is_in_target = [&Rtol,&first_times_xtarg] (pt2_t x) -> bool {
		return !(x - pt2_t{ first_times_xtarg[0], 0 }) < Rtol*Rtol;
	};

	while (not _thread.do_quit) {
		
		pysimul_mutex_unlock(&_thread);
		if (pause)
			::usleep(1000000);
		pysimul_mutex_lock(&_thread);
		
		// Let's compute one entire trajectory
		
		t = 0;
		pt2_t x;
		
		auto init_pos = [&] () -> void {
			reinit:
			x = pt2_t{
				.x = init_pos_sigma * normal_distrib_gen(rng),
				.y = init_pos_sigma * normal_distrib_gen(rng)
			};
			if (check_is_in_target(x)) {	// we must not reinit the particle onto the target
		//		goto reinit;
				t = 0;	// equivalent to forget the trajectory; introduces a bias towards smaller FPT
			}
		};
		init_pos();
		
		while (true) {
			
			// poissonian resetting
			if (unif01(rng) < proba_reset_step) 
				init_pos();
			
			// Langevin equation implementation
			x = x + sqrt(2 * D * Δt) * vec2_t{ .x = normal_distrib_gen(rng), .y = normal_distrib_gen(rng) };
			
			// target reached -> first passage time
			if (check_is_in_target(x)) {
				first_time.push_back(t); // register the FPT
				break;
			}
			
			t += Δt;
		}
		
		n_trajectories++;
	}
	
	pysimul_mutex_unlock(&_thread);
	return nullptr;
}
