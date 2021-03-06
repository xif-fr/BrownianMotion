#include <fmt/core.h>
#include "pysimul-common.h"
#include <math.h>
#include <random>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "vecs2.hpp"
#include <unistd.h>

#define LANGEVIN_OVERDAMPED
#define TARGET_2D_CYL
//#define ENABLE_SURVIVAL_PROBABILITIES_INTERVAL
//#define FPT_JUMP_ACROSS
#define FPT_INTERVAL
#define ENABLE_PERIODICAL_RESET
#define XTARG_ONE_VARIABLE // for langevin-ft-automated-automated.ipynb
#define INPUT_DATA_FILE

#ifdef XTARG_ONE_VARIABLE
constexpr size_t N_targets = 1;
const size_t pysimul_N = 0;
#else
constexpr size_t N_targets = 10;
const size_t pysimul_N = N_targets;
#endif

// transformation from 2D coordinates to 1D to use for the 1D target (usually just x.x)
#define POS1D(point) (point.x)

#ifdef INPUT_DATA_FILE
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

void* comp_thread (void* _data) {
	simul_thread_info_t& _thread = *(simul_thread_info_t*)_data;
	simul_thread_raii_mutex_unlock _mutex_unlock(&_thread);

	// 2D circular target with "tolerence radius" Rtol around x=x_targ,y=0
	// By default, only 1 dimension is considered and targets are simply x=cst lines
	#ifdef TARGET_2D_CYL
	double Rtol = 0.1;
	_register_var(_thread, "2D-Rtol", &Rtol);
	#if defined(FPT_JUMP_ACROSS) || defined(FPT_DEMISPACE) || defined(ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE)
		#error "2D target is only compatible with *_INTERVAL flags"
	#endif
	#endif
	
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
	#if defined(FPT_JUMP_ACROSS) || defined(FPT_INTERVAL) || defined(FPT_DEMISPACE)
	std::array<std::vector<double>,N_targets> first_times;
	// 
	// Only run-time-parametrizable target (for langevin-ft-automated-automated.ipynb) :
	#ifdef XTARG_ONE_VARIABLE
	std::array<double, N_targets> first_times_xtarg = { 0.1 };
	_register_var(_thread, "first_times_xtarg", &first_times_xtarg[0]);
	_register_var(_thread, "first_times", &first_times[0]);
	// 
	// Multiple fixed targets (langevin-ft-automated.ipynb) :
	#else
	constexpr std::array<double, N_targets> first_times_xtarg = { 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20 };
	auto _first_times_xtarg = first_times_xtarg;
	_register_Narray(_thread, "first_times_xtarg", _first_times_xtarg);
	for (size_t i = 0; i < N_targets; i++) {
		auto s = fmt::format("first_times-{}", i);
		_register_var(_thread, s.c_str(), &first_times[i]);
	}
	#endif
	// For 1D : size of the target (in theory 0, in practice must be small but non-zero for FPT_JUMP_ACROSS or FPT_INTERVAL) :
	constexpr double xtarg_tol = 0.001;
	_register_const(_thread, "xtarg_tol", xtarg_tol);
	#endif
	
	double t = 0;
	size_t step = 0;
	constexpr size_t trajectory_step_max_length = 1000000000;
	size_t n_trajectories = 0;
	_register_var(_thread, "n_trajectories", &n_trajectories);
	
	constexpr double Δt = 1e-7;//50 * 1.5e-6;
	_register_const(_thread, "Delta_t", Δt);
	uint8_t pause = 0;
	_register_var(_thread, "pause", &pause);
	
	#ifndef INPUT_DATA_FILE

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

	double init_pos_sigma = 0.1; // gaussian distribution of initial position
	_register_var(_thread, "x0sigma", &init_pos_sigma);
		
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

	#else

	void* faddr_base = nullptr;
	size_t fsize = 0;
	std::string file_path = "traj_data.bin";
	_register_var(_thread, "file_path", &file_path);
	{
		int fd = ::open(file_path.c_str(), O_RDONLY);
		if (fd == -1) { ::perror("can't open data file"); return nullptr; }
		struct stat sb;
		::fstat(fd, &sb);
		fsize = sb.st_size;
		faddr_base = ::mmap(NULL, fsize, PROT_READ, MAP_PRIVATE, fd, /*offset*/0);
		if (faddr_base == MAP_FAILED) { ::perror("can't mmap data file"); return nullptr; }
		::close(fd);
		::puts("successfully opened trajectory file");
	}
	void* faddr_curr = faddr_base;
	auto file_read_point = [&] () -> const pt2_t& {
		int8_t* faddr_next = (int8_t*)faddr_curr + sizeof(double)*2;
		if (faddr_next - (int8_t*)faddr_base >= fsize) 
			throw std::out_of_range("end of file");
		const pt2_t& x = *(const pt2_t*)faddr_curr;
		faddr_curr = (void*)faddr_next;
		return x;
	};
	auto file_close = [&] () {
		::munmap(faddr_base, fsize);
	};

	try {
	#endif
	
	while (not _thread.do_quit) {
		
		if (_thread.regular_callback)
			_thread.regular_callback(_thread.id_for_callback, step, t);
		
		pysimul_mutex_unlock(&_thread);
		if (pause)
			::usleep(1000000);
		pysimul_mutex_lock(&_thread);
		
		// Let's compute one entire trajectory
		
		step = 0;
		t = 0;
		
		#ifndef INPUT_DATA_FILE
		pt2_t x;
		#ifndef LANGEVIN_OVERDAMPED
		vec2_t v;
		#endif
		#endif
		
		double max_x_reached = 0.;
		bool survdist_pos_done = false;
		uint8_t target_reached = N_targets;
		#if defined(FPT_JUMP_ACROSS) || defined(FPT_INTERVAL) || defined(FPT_DEMISPACE)
		target_reached = 0;
		#endif
		
		// when using ENABLE_SURVIVAL_PROBABILITIES_INTERVAL, we check if a target is reached
		// every step, and keep track of reached targets in `survdist_survived`; it is slower
		// than ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE, but the only accurate way when there
		// is resetting with init_pos_sigma≠0
		#ifdef ENABLE_SURVIVAL_PROBABILITIES_INTERVAL
		std::array<bool,survdist_Ndx> survdist_survived;
		survdist_survived.fill(true);
		#endif
		
		pt2_t last_x;
		std::array<bool,N_targets> targets_reached;
		targets_reached.fill(false);

		#ifndef INPUT_DATA_FILE

		auto init_pos = [&] () -> void {
			x = pt2_t{
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

		#else

		// waiting for a reset
		pt2_t x;
		do {
			x = file_read_point();
		} while (not isnan(x.x));

		#endif
		
		#ifdef ENABLE_PERIODICAL_RESET
		uint32_t n_period = 1;
		#endif
		
		#if defined(ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE) || defined(ENABLE_SURVIVAL_PROBABILITIES_INTERVAL)
		while (t < survdist_max_t or (target_reached < N_targets and step < trajectory_step_max_length)) {
		#else
		while (target_reached < N_targets and step < trajectory_step_max_length) {
		#endif
			
			#ifndef INPUT_DATA_FILE

				// resetting
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
				last_x = x;
				
				// Langevin equation implementation
				vec2_t f_alea = sqrt(2 * γ * T / Δt) * vec2_t{ .x = normal_distrib_gen(rng), .y = normal_distrib_gen(rng) };
				#ifndef LANGEVIN_OVERDAMPED
				vec2_t a = -γ*v + f_alea;
				#else
				vec2_t v = f_alea / γ;
				#endif
				#ifndef LANGEVIN_OVERDAMPED
				v += a / part_m * Δt;
				#endif
				x = x + v * Δt; // Δx ~= sqrt(D.Δt), 1e-3 for D=1 and Δt=1e-6

			#else

				last_x = x;
				x = file_read_point();
				if (isnan(last_x.x))
					last_x = x;

				if (isnan(x.x)) { // reset point, skipped but must be took into account for FPT_JUMP_ACROSS
					x = file_read_point();
					last_x = x;
				}

			#endif
			
			// record x position reached
			#if defined(ENABLE_SURVIVAL_PROBABILITIES_DEMISPACE) || defined(FPT_DEMISPACE)
			if (POS1D(x) > max_x_reached) {
				max_x_reached = POS1D(x);
				#ifdef FPT_DEMISPACE
				if (max_x_reached >= first_times_xtarg[target_reached] and target_reached < N_targets) {
					first_times[target_reached].push_back(t);
					target_reached++;
				}
				#endif
			}
			#endif
			
			// first passage time distribution for each target at x=first_times_xtarg[i]
			#if defined(FPT_JUMP_ACROSS) || defined(FPT_INTERVAL)
			for (uint8_t i = 0; i < N_targets; i++) {
				if (not targets_reached[i]) {
					#ifdef FPT_JUMP_ACROSS
						bool a = POS1D(x) < first_times_xtarg[i];
						bool b = first_times_xtarg[i] < POS1D(last_x);
						bool target = (a and b) or (not a and not b);
					#endif
					#ifdef FPT_INTERVAL
						#ifndef TARGET_2D_CYL
						bool target = std::abs(POS1D(x) - first_times_xtarg[i]) < xtarg_tol;
						#else
						bool target = !(x - pt2_t{ first_times_xtarg[i], 0 }) < Rtol*Rtol;
						#endif
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
			// keeping track of reached targets for survival distributions
			constexpr double dx = survdist_max_x/survdist_Ndx;
			#ifdef TARGET_2D_CYL
			int64_t x_k_inf = std::max<int64_t>( std::floor( (POS1D(x)-Rtol) / dx ), 0 );
			int64_t x_k_sup = std::min<int64_t>( std::floor( (POS1D(x)+Rtol) / dx ), survdist_Ndx );
			for (int64_t x_k = x_k_inf; x_k < x_k_sup; x_k++) {
				if (!(x - pt2_t{ x_k*dx, 0 }) < Rtol*Rtol)
					survdist_survived[x_k] = false;
			}
			#else
			int64_t x_k = std::floor( POS1D(x) / dx );
			if (x_k >= 0 and x_k < survdist_Ndx)
				survdist_survived[x_k] = false;
			#endif
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
			
			// let's check if if the particle reached the target at specified x
			// and fill the survival distribution as a function of time accordingly
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

	#ifdef INPUT_DATA_FILE
	} catch (std::out_of_range&) {
		::puts("data exhausted");
		pysimul_mutex_unlock(&_thread);
		pause = 1;
		file_close();
		while (not _thread.do_quit) 
			::usleep(1000000);
	}
	#endif

	return nullptr;
}
