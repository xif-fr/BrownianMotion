#include "pysimul-common.h"
#include <random>
#include <vector>
#include <algorithm>
#include "vecs2.hpp"
#include <unistd.h>

// This is the boiled-down version of langevin-survival.cpp with the following flags :
// LANGEVIN_OVERDAMPED, TARGET_2D_CYL, FPT_INTERVAL
// Additionally supports anisotropic initial position distribution

const size_t pysimul_N = 0;

void* comp_thread (void* _data) {
	simul_thread_info_t& _thread = *(simul_thread_info_t*)_data;
	
	std::random_device _rd;
	std::mt19937 rng (_rd());
	std::normal_distribution<> normal_distrib_gen (0, 1); // mean, std
	std::uniform_real_distribution<> unif01 (0, 1);

	double t = 0;
	size_t n_trajectories = 0;
	_register_var(_thread, "n_trajectories", &n_trajectories);
	constexpr double Δt = 1e-7;
	_register_const(_thread, "Delta_t", Δt);
	uint8_t pause = 0;
	_register_var(_thread, "pause", &pause);

	// Diffusion
	double D;
	_register_var(_thread, "D", &D);
	const double Ɣ = sqrt(2 * D * Δt);

	// Gaussian distribution of initial position
	double init_pos_sigma_x = 0., init_pos_sigma_y = 0.;
	_register_var(_thread, "x0sigma_x", &init_pos_sigma_x);
	_register_var(_thread, "x0sigma_y", &init_pos_sigma_y);

	#define ENABLE_PERIODICAL_RESET

	// Poissonian resetting
	#ifdef ENABLE_POISSON_RESET
	double reset_rate;
	_register_var(_thread, "reset_rate", &reset_rate);
	const double proba_reset_step = Δt * reset_rate;
	#endif
	
	// Periodical resetting
	#ifdef ENABLE_PERIODICAL_RESET
	double reset_period = 1;
	_register_var(_thread, "reset_period", &reset_period);
	#endif

	// First passage time at the target
	std::vector<double> first_times;
	double xtarg;
	_register_var(_thread, "xtarg", &xtarg);
	_register_var(_thread, "first_times", &first_times);

	// 2D circular target with "tolerence radius" Rtol around x=x_targ,y=0
	double Rtol;
	_register_var(_thread, "2D-Rtol", &Rtol);
	auto check_is_in_target = [&Rtol,&xtarg] (pt2_t x) -> bool {
		return !(x - pt2_t{ xtarg, 0 }) < Rtol*Rtol;
	};

	while (not _thread.do_quit) {
		
		pysimul_mutex_unlock(&_thread);
		if (pause)
			::usleep(1000000);
		pysimul_mutex_lock(&_thread);
		
		// Let's compute one entire trajectory
		
		t = 0;
		pt2_t x;
		#ifdef ENABLE_PERIODICAL_RESET
		uint32_t n_period = 1;
		#endif
		
		auto init_pos = [&] () -> void {
			x = pt2_t{
				.x = init_pos_sigma_x * normal_distrib_gen(rng),
				.y = init_pos_sigma_y * normal_distrib_gen(rng)
			};
			// note : the particle can be reinitialized on the target, we must not prevent that
		};
		init_pos();
		
		while (true) {
			
			// poissonian resetting
			#ifdef ENABLE_POISSON_RESET
			if (unif01(rng) < proba_reset_step) 
				init_pos();
			#endif
			#ifdef ENABLE_PERIODICAL_RESET
			if (t > reset_period*n_period) {
				init_pos();
				n_period++;
			}
			#endif
			
			// Langevin equation implementation
			x = x + Ɣ * vec2_t{ .x = normal_distrib_gen(rng), .y = normal_distrib_gen(rng) };
			
			// target reached -> first passage time
			if (check_is_in_target(x)) {
				first_times.push_back(t); // register the FPT
				break;
			}
			
			t += Δt;
		}
		
		n_trajectories++;
	}
	
	pysimul_mutex_unlock(&_thread);
	return nullptr;
}
