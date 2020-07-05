#include <fmt/core.h>
#include "pysimul-common.h"
#include <random>
#include <vector>

constexpr size_t N_targets = 10;
const size_t pysimul_N = N_targets;

void* comp_thread (void* _data) {
	simul_thread_info_t& _thread = *(simul_thread_info_t*)_data;

	#define LANGEVIN_OVERDAMPED
	#ifndef LANGEVIN_OVERDAMPED
	double part_m = 10;
	_register_var(_thread, "part_m", &part_m);
	#endif
	double γ = 600;
	_register_var(_thread, "gamma", &γ);
	double T = 10;
	_register_var(_thread, "T", &T);
	
	std::random_device _rd;
	std::mt19937 rng (_rd());
	std::normal_distribution<> normal_distrib_gen (0, 1); // mean, std
	std::uniform_real_distribution<> unif01 (0, 1);
	
//	#define ENABLE_SURVIVAL_PROBABILITIES_INTERVAL
	#if defined(ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE) || defined(ENABLE_SURVIVAL_PROBABILITIES_INTERVAL)
	// Survival probability as a function of time at a given position
	constexpr size_t survdist_Ndt = 100;
	constexpr double survdist_max_t = 5.0;
	_register_const(_thread, "survdist_max_t", survdist_max_t);
	constexpr double survdist_time_pos = 0.3; // let's put the target @ x=survdist_time_pos
	_register_const(_thread, "survdist_time_pos", survdist_time_pos);
	std::array<uint64_t,survdist_Ndt> survdist_time_acc;
	survdist_time_acc.fill(0);
	uint64_t survdist_time_samples = 0;
	_register_distrib(_thread, "survdist_time", survdist_time_acc, &survdist_time_samples);
	
	// Survival probability as a function of target position at a given time
	constexpr size_t survdist_Ndx = 400;
	constexpr double survdist_max_x = 0.8;
	_register_const(_thread, "survdist_max_x", survdist_max_x);
	constexpr double survdist_pos_time = 1.0; // let's sample at t=survdist_pos_time
	_register_const(_thread, "survdist_pos_time", survdist_pos_time);
	std::array<uint64_t,survdist_Ndx> survdist_pos_acc;
	survdist_pos_acc.fill(0);
	uint64_t survdist_pos_samples = 0;
	_register_distrib(_thread, "survdist_pos", survdist_pos_acc, &survdist_pos_samples);
	#endif
	
	// First passage time at the target @ x=survdist_time_pos
	#define FPT_JUMP_ACROSS
	std::array<std::vector<double>,N_targets> first_times;
	constexpr std::array<double, N_targets> first_times_xtarg = { 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50 };
	constexpr double xtarg_tol = 0.001;
	_register_const(_thread, "xtarg_tol", xtarg_tol);
	for (size_t i = 0; i < N_targets; i++) {
		auto s = fmt::format("first_times-{}", i);
		_register_var(_thread, s.c_str(), &first_times[i]);
	}
	auto _first_times_xtarg = first_times_xtarg;
	_register_Narray(_thread, "first_times_xtarg", _first_times_xtarg);
	
	double t = 0;
	size_t step = 0;
	constexpr size_t trajectory_step_max_length = 1000000000;
	size_t n_trajectories = 0;
	_register_var(_thread, "n_trajectories", &n_trajectories);
	
	constexpr double Δt = 1e-6;//50 * 1.5e-6;
	_register_const(_thread, "Delta_t", Δt);
	uint8_t pause = 0;
	_register_var(_thread, "pause", &pause);
	double t_pause = Inf;
	_register_var(_thread, "t_pause", &t_pause);
	
	double init_pos_sigma = 0.1; // gaussian distribution of initial position
	_register_var(_thread, "x0sigma", &init_pos_sigma);
	
	#define ENABLE_PERIODICAL_RESET
	
	// Poissonian resetting
	#ifdef ENABLE_POISSON_RESET
	double reset_rate = 0.001 / Δt;
	_register_var(_thread, "reset_rate", &reset_rate);
	const double proba_reset_step = Δt * reset_rate;
	#endif
	
	// Periodical resetting
	#ifdef ENABLE_PERIODICAL_RESET
	double reset_period = 1;
	_register_var(_thread, "reset_period", &reset_period);
	#endif
	
	while (not _thread.do_quit) {
		
		if (_thread.regular_callback)
			_thread.regular_callback(_thread.id_for_callback, step, t);
		
		::pthread_mutex_unlock(&_thread.mutex_global);
		::usleep(pause ? 1000000 : 100);
		::pthread_mutex_lock(&_thread.mutex_global);
		
		// Let's compute one entire trajectory
		
		step = 0;
		t = 0;
		
		pt2_t x;
		#ifndef LANGEVIN_OVERDAMPED
		vec2_t v;
		#endif
		
		auto init_pos = [&] () -> void {
			x =  pt2_t{
				.x = init_pos_sigma * normal_distrib_gen(rng),
				.y = init_pos_sigma * normal_distrib_gen(rng)
			};
			#ifndef LANGEVIN_OVERDAMPED
			v = (vec2_t)vecO_t{				// n'a aucun impact, en tout cas à b=infini
				.r = sqrt(2*T/part_m),
				.θ = unif01(rng)*2*π
			};
			#endif
		};
		init_pos();
		
		double max_x_reached = 0.;
		bool survdist_pos_done = false;
		uint8_t target_reached = N_targets;
		#if defined(FPT_JUMP_ACROSS) || defined(FPT_INTERVAL) || defined(FPT_DEMISPACE)
		target_reached = 0;
		#endif
		
		#ifdef ENABLE_SURVIVAL_PROBABILITIES_INTERVAL
		std::array<bool,survdist_Ndx> survdist_survived;
		survdist_survived.fill(true);
		#endif
		
		double last_x;
		std::array<bool,N_targets> targets_reached;
		targets_reached.fill(false);
		
		#ifdef ENABLE_PERIODICAL_RESET
		uint32_t n_period = 1;
		#endif
		
		#if defined(ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE) || defined(ENABLE_SURVIVAL_PROBABILITIES_INTERVAL)
		while (t < survdist_max_t or (target_reached < N_targets and step < trajectory_step_max_length)) {
		#else
		while (target_reached < N_targets and step < trajectory_step_max_length) {
		#endif
			
			#ifdef ENABLE_POISSON_RESET
			if (unif01(rng) < proba_reset_step) {
				init_pos();
			}
			#endif
			#ifdef ENABLE_PERIODICAL_RESET
			if (t > reset_period*n_period) {
				init_pos();
				n_period++;
			}
			#endif
			last_x = x.x;
			
			vec2_t f_alea = sqrt(2 * γ * T / Δt) * vec2_t{ .x = normal_distrib_gen(rng),/* .y = normal_distrib_gen(rng)*/ };
			#ifndef LANGEVIN_OVERDAMPED
//			vec2_t f_ext = well_k * (well_center - x);
			vec2_t a = -γ*v + f_alea;// + f_ext;
			#else
			vec2_t v = f_alea / γ;
			#endif
			
			#ifndef LANGEVIN_OVERDAMPED
			v += a / part_m * Δt;
			#endif
			x = x + v * Δt; // Δx ~= sqrt(D.Δt)
			
			// record x position reached
			#if defined(ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE) || defined(FPT_DEMISPACE)
			if (x.x > max_x_reached) {
				max_x_reached = x.x;
				#ifdef FPT_DEMISPACE
				if (max_x_reached >= first_times_xtarg[target_reached] and target_reached < N_targets) {
					first_times[target_reached].push_back(t);
					target_reached++;
				}
				#endif
			}
			#endif
			#if defined(FPT_JUMP_ACROSS) || defined(FPT_INTERVAL)
			for (uint8_t i = 0; i < N_targets; i++) {
				if (not targets_reached[i]) {
					#ifdef FPT_JUMP_ACROSS
					bool a = x.x < first_times_xtarg[i];
					bool b = first_times_xtarg[i] < last_x;
					bool target = (a and b) or (not a and not b);
					#endif
					#ifdef FPT_INTERVAL
					bool target = std::abs(x.x - first_times_xtarg[i]) < xtarg_tol;
					#endif
					if (target) {
						targets_reached[i] = true;
						first_times[i].push_back(t);
						target_reached = std::count(targets_reached.begin(), targets_reached.end(), true);
					}
				}
			}
			#endif
			
			#ifdef ENABLE_SURVIVAL_PROBABILITIES_INTERVAL
			constexpr double dx = survdist_max_x/survdist_Ndx;
			int64_t x_k = std::floor( x.x / dx );
			if (x_k >= 0 and x_k < survdist_Ndx)
				survdist_survived[x_k] = false;
			#endif
			#if defined(ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE) || defined(ENABLE_SURVIVAL_PROBABILITIES_INTERVAL)
			
			// let's check which targets the particle reached at specified t
			// and fill the survival distribution as a function of target position accordingly
			if (t > survdist_pos_time and not survdist_pos_done) {
				#ifdef ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE
				constexpr double dx = survdist_max_x/survdist_Ndx;
				int64_t x_k = std::floor( max_x_reached / dx );
				for (size_t k = 0; k < survdist_Ndx; k++) {
					if (k >= x_k) // survived all targets at < x_k; only compatible with particule resetting at x=0 (no backwards target reaching)
						survdist_pos_acc[k]++;
				}
				#endif
				#ifdef ENABLE_SURVIVAL_PROBABILITIES_INTERVAL
				for (size_t k = 0; k < survdist_Ndx; k++) {
					if (survdist_survived[k])
						survdist_pos_acc[k]++;
				}
				#endif
				survdist_pos_samples++;
				survdist_pos_done = true;
			}
			
			//
			constexpr double dt = survdist_max_t/survdist_Ndt;
			uint64_t t_k = std::floor( t / dt );
			if (t_k < survdist_Ndt) {
				#ifdef ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE
				if (max_x_reached < survdist_time_pos)  // survived targets at survdist_time_pos
					survdist_time_acc[t_k]++;
				#endif
				#ifdef ENABLE_SURVIVAL_PROBABILITIES_INTERVAL
				constexpr size_t x_k = survdist_time_pos/dx;
				if (survdist_survived[x_k])
					survdist_time_acc[t_k]++;
				#endif
				survdist_time_samples++;
			}
			
			#endif

			
			t += Δt;
			step++;
		}
		
		n_trajectories++;
	}
	
	::pthread_mutex_unlock(&_thread.mutex_global);
	return nullptr;
}
