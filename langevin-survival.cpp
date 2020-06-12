#include <fmt/core.h>
#include "pysimul-common.h"
#include <random>
#include <vector>

constexpr size_t N_targets = 1;
const size_t pysimul_N = N_targets;

void* comp_thread (void* _data) {
	simul_thread_info_t& _thread = *(simul_thread_info_t*)_data;

	double part_m = 10;
	_register_var(_thread, "part_m", &part_m);
	double γ = 600;
	_register_var(_thread, "gamma", &γ);
	double T = 10;
	_register_var(_thread, "T", &T);
	
	std::random_device _rd;
	std::mt19937 rng (_rd());
	std::normal_distribution<> normal_distrib_gen (0, 1); // mean, std
	std::uniform_real_distribution<> unif01 (0, 1);
	
//	double well_k = 1e4;
//	pt2_t well_center = {0.5, 0.5};
//	_register_var(_thread, "well_k", &well_k); _register_var(_thread, "well_x", &well_center.x); _register_var(_thread, "well_y", &well_center.y);
	
	#define ENABLE_SURVIVAL_PROBABILITIES
	#ifdef ENABLE_SURVIVAL_PROBABILITIES
	// Survival probability as a function of time
	constexpr size_t survdist_Ndt = 100;
	constexpr double survdist_max_t = 5.0;
	_register_const(_thread, "survdist_max_t", survdist_max_t);
	constexpr double survdist_time_pos = 0.3; // let's put the target @ x=survdist_time_pos
	_register_const(_thread, "survdist_time_pos", survdist_time_pos);
	std::array<uint64_t,survdist_Ndt> survdist_time_acc;
	survdist_time_acc.fill(0);
	uint64_t survdist_time_samples = 0;
	_register_distrib(_thread, "survdist_time", survdist_time_acc, &survdist_time_samples);
	
	// Survival probability as a function of target position
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
	std::array<std::vector<double>,N_targets> first_times;
//	constexpr std::array<double, N_targets> first_times_xtarg = { 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50 };
	constexpr std::array<double, N_targets> first_times_xtarg = { 0.30 };
	for (size_t i = 0; i < N_targets; i++) {
		auto s = fmt::format("first_times-{}", i);
		_register_var(_thread, s.c_str(), &first_times[i]);
	}
	auto _first_times_xtarg = first_times_xtarg;
	_register_Narray(_thread, "first_times_xtarg", _first_times_xtarg);
	
	double t = 0;
	size_t step = 0;
	constexpr size_t trajectory_step_max_length = 5000000;
	size_t n_trajectories = 0;
	_register_var(_thread, "n_trajectories", &n_trajectories);
	
	constexpr double Δt = 50 * 1.5e-6;
	_register_const(_thread, "Delta_t", Δt);
	uint8_t pause = 0;
	_register_var(_thread, "pause", &pause);
	double t_pause = Inf;
	_register_var(_thread, "t_pause", &t_pause);
	
	double init_pos_sigma = 0.1; // gaussian distribution of initial position
	_register_var(_thread, "x0sigma", &init_pos_sigma);
	
	// Poissonian resetting
	#define ENABLE_POISSON_RESET
	#ifdef ENABLE_POISSON_RESET
	double reset_rate = 0.001 / Δt;
	_register_var(_thread, "reset_rate", &reset_rate);
	const double proba_reset_step = Δt * reset_rate;
	#endif
	
	while (not _thread.do_quit) {
		
		if (_thread.regular_callback)
			_thread.regular_callback(_thread.id_for_callback, step, t);
		
		::pthread_mutex_unlock(&_thread.mutex_global);
		::usleep(pause ? 1000000 : 10);
		::pthread_mutex_lock(&_thread.mutex_global);
		
		// Let's compute one entire trajectory
		
		step = 0;
		t = 0;
		
		pt2_t x;
		vec2_t v;
		
		auto init_pos = [&] () -> void {
			x =  pt2_t{
				.x = init_pos_sigma * normal_distrib_gen(rng),
				.y = init_pos_sigma * normal_distrib_gen(rng)
			};
			v = (vec2_t)vecO_t{				// n'a aucun impact, en tout cas à b=infini
				.r = sqrt(2*T/part_m),
				.θ = unif01(rng)*2*π
			};
		};
		init_pos();
		
		double max_x_reached = 0.;
		bool survdist_pos_done = false;
		uint8_t target_reached = 0;
		
		#ifdef ENABLE_SURVIVAL_PROBABILITIES
		while (t < survdist_max_t or (target_reached < N_targets and step < trajectory_step_max_length)) {
		#else
		while (target_reached < N_targets and step < trajectory_step_max_length) {
		#endif
			
			#ifdef ENABLE_POISSON_RESET
			if (unif01(rng) < proba_reset_step) {
				init_pos();
			}
			#endif
			
			vec2_t f_alea = sqrt(2 * γ * T / Δt) * vec2_t{ .x = normal_distrib_gen(rng),/* .y = normal_distrib_gen(rng)*/ };
//			vec2_t f_ext = well_k * (well_center - x);
			vec2_t a = -γ*v + f_alea;// + f_ext;
			
			v += a / part_m * Δt;
			x = x + v * Δt;
			
			// record x position reached
			if (x.x > max_x_reached) {
				max_x_reached = x.x;
				if (max_x_reached >= first_times_xtarg[target_reached] and target_reached < N_targets) {
					first_times[target_reached].push_back(t);
					target_reached++;
				}
			}
			
			#ifdef ENABLE_SURVIVAL_PROBABILITIES
			// let's check which targets the particle reached at specified t
			// and fill the survival distribution as a function of target position accordingly
			if (t > survdist_pos_time and not survdist_pos_done) {
				constexpr double dx = survdist_max_x/survdist_Ndx;
				int64_t x_k = std::floor( max_x_reached / dx );
				for (size_t k = 0; k < survdist_Ndx; k++) {
					if (k >= x_k) // survived all targets at < x_k
						survdist_pos_acc[k]++;
				}
				survdist_pos_samples++;
				survdist_pos_done = true;
			}
			
			constexpr double dt = survdist_max_t/survdist_Ndt;
			uint64_t t_k = std::floor( t / dt );
			if (t_k < survdist_Ndt) {
				if (max_x_reached < survdist_time_pos)  // survived targets at survdist_time_pos
					survdist_time_acc[t_k]++;
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
