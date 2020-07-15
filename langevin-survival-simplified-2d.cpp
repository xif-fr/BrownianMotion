#include <fmt/core.h>
#include "pysimul-common.h"
#include <random>
#include <vector>
#include <algorithm>

// This is the boiled-down version of langevin-survival.cpp with the following flags :
// LANGEVIN_OVERDAMPED, TARGET_2D_CYL, FPT_INTERVAL, ENABLE_POISSON_RESET

constexpr size_t N_targets = 7;
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
	constexpr double Δt = 1e-6;
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
	std::array<std::vector<double>,N_targets> first_times;
	constexpr std::array<double, N_targets> first_times_xtarg = { 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40 };
	for (size_t i = 0; i < N_targets; i++) {
		auto s = fmt::format("first_times-{}", i);
		_register_var(_thread, s.c_str(), &first_times[i]);
	}
	auto _first_times_xtarg = first_times_xtarg;
	_register_Narray(_thread, "first_times_xtarg", _first_times_xtarg);

	// 2D circular target with "tolerence radius" Rtol around x=x_targ,y=0
	double Rtol;
	_register_var(_thread, "2D-Rtol", &Rtol);
	auto check_is_in_target = [&] (pt2_t x, double Rtarg) -> bool {
		return !(x - pt2_t{ Rtarg, 0 }) < Rtol*Rtol;
	};

	while (not _thread.do_quit) {
		
		pysimul_mutex_unlock(&_thread);
		if (pause)
			::usleep(1000000);
		pysimul_mutex_lock(&_thread);
		
		// Let's compute one entire trajectory
		
		t = 0;
		pt2_t x;
		
		bool all_targets_reached = false;
		std::array<bool,N_targets> targets_reached;
		targets_reached.fill(false);
		auto set_target_reached = [&] (uint8_t i_targ) {
			targets_reached[i_targ] = true;
			all_targets_reached = (std::count(targets_reached.begin(), targets_reached.end(), true) == N_targets);
		};

		auto init_pos = [&] () -> void {
			x = pt2_t{
				.x = init_pos_sigma * normal_distrib_gen(rng),
				.y = init_pos_sigma * normal_distrib_gen(rng)
			};
			for (uint8_t i = 0; i < N_targets; i++) {
				if (check_is_in_target(x, first_times_xtarg[i])) 
					set_target_reached(i); // when initialized into the target, we forget about it
			}
		};
		init_pos();
		
		while (not all_targets_reached) {
			
			// poissonian resetting
			if (unif01(rng) < proba_reset_step) 
				init_pos();
			
			// Langevin equation implementation
			x = x + sqrt(2 * D * Δt) * vec2_t{ .x = normal_distrib_gen(rng), .y = normal_distrib_gen(rng) };
			
			// first passage time distribution for each target at x=first_times_xtarg[i]
			for (uint8_t i = 0; i < N_targets; i++) {
				if (not targets_reached[i]) {
					if (check_is_in_target(x, first_times_xtarg[i])) {
						set_target_reached(i);
						first_times[i].push_back(t); // register the FPT
					}
				}
			}
			
			t += Δt;
		}
		
		n_trajectories++;
	}
	
	pysimul_mutex_unlock(&_thread);
	return nullptr;
}
