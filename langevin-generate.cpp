#include <fmt/core.h>
#include <random>
#include <vector>
#include <algorithm>
#include <utility>
#include <math.h>
#define Inf std::numeric_limits<double>::infinity()
#define NaN std::numeric_limits<double>::signaling_NaN()
#include <inttypes.h>
#include "vecs2.hpp"
#include <unistd.h>
#include <fcntl.h>

#define LANGEVIN_OVERDAMPED
#define ENABLE_PERIODICAL_RESET
//#define RESET_WITH_TRAPPING
//#define SPLIT_FILES

/*********************************************************************************************
/ Simulation of a brownian trajectory with resetting, using Langevin's equation in 2D,
/ generating a simple file. Standalone command line programm. Supports periodical resetting
/ of period `reset_period`, and poissonian resetting with rate `proba_reset_step/Δt`.
/
/ Resetting can be either ideal (immediate reinitialization of the particle's position with
/ variance `init_pos_sigma`²) or be done by simulating an optical trap, that is an harmonic
/ well with stiffness `k_trap`, during a time `t_trap`, bringing back the particule to the
/ origin (with variance σ²=T/k), when RESET_WITH_TRAPPING is enabled.
/
/ A single file is saved by default, containing successive x and y positions as 64 bits
/ floating point numbers (and (NaN,NaN) to indicate end of resetting), and ready to be used by
/ `langevin-survival.cpp` with INPUT_DATA_FILE enabled when ideal resetting. If SPLIT_FILES is
/ enabled, a sequence of small (~20MB) files is saved, containing successive x and y positions
/ and a byte indicating if the trapping well is active (x0,y0,trapping0,y1,y1,trapping1,...),
/ ready to be used by `exp-data-diffus-analysis.ipynb` when optical trap resetting.
/ ********************************************************************************************/

#include <signal.h>
volatile bool continue_running = true;
void interupt_handler (int) {
	continue_running = false;
}

int main (int argc, char const* argv[]) {

	std::string fname_base = "langevin-trap-traj-xyc-";
	if (argc > 1)
		fname_base = argv[1];
	#ifdef SPLIT_FILES
	size_t f_num = 1;
	if (argc > 2)
		f_num = atoi(argv[2]);
	std::string fname_format = fname_base+"{}";
	#endif
	int f_fd = -1;
	
	::signal(SIGINT, interupt_handler);
	
	double t = 0;
	size_t step = 0;
	size_t resets = 0;
	
	constexpr double Δt = 1e-5;
	
	#ifndef LANGEVIN_OVERDAMPED
	constexpr double part_m = 0.0001;
	#endif
	constexpr double γ = 0.2;
	constexpr double T = 1;
	
	std::random_device _rd;
	std::mt19937 rng (_rd()); // or set any seed you want
	std::normal_distribution<> normal_distrib_gen (0, 1); // mean, std
	std::uniform_real_distribution<> unif01 (0, 1);
	
	std::function<void()> reset_init = [&] () -> void {};
	std::function<bool()> reset_do = [&] () -> bool { return false; };
	// Poissonian resetting
	#ifdef ENABLE_POISSON_RESET
	const double proba_reset_step = 0.001; // = Δt * reset_rate;
	reset_do = [&] () -> bool {
		return unif01(rng) < proba_reset_step;
	};
	#endif
	// Periodical resetting
	#ifdef ENABLE_PERIODICAL_RESET
	double reset_period = 0.2;
	double t_reset_next = Inf;
	reset_init = [&] () -> void {
		t_reset_next = t + reset_period;
	};
	reset_do = [&] () -> bool {
		return t > t_reset_next;
	};
	#endif
	
	pt2_t x = {0,0};
	#ifndef LANGEVIN_OVERDAMPED
	vec2_t v = {0,0};
	#endif

	bool trapping = false;
	#ifdef RESET_WITH_TRAPPING
	constexpr pt2_t x0_trap = {0, 0};
	constexpr double k_trap = 10;
	constexpr double t_trap = 0.08;
	double t_trap_end = Inf;
	auto init_pos = [&] () -> void {
		trapping = true;
		t_trap_end = t + t_trap;
	};
	#else
	constexpr double init_pos_sigma = 0.31622776601683794; // gaussian distribution of initial position
	auto init_pos = [&] () -> void {
		x = pt2_t{
			.x = init_pos_sigma * normal_distrib_gen(rng),
			.y = init_pos_sigma * normal_distrib_gen(rng)
		};
		#ifndef LANGEVIN_OVERDAMPED
		v = (vec2_t)vecO_t{				// n'a aucun impact, en tout cas à b=infini
			.r = sqrt(2*T/part_m),
			.θ = unif01(rng)*2*M_PI
		};
		#endif
	};
	reset_init();
	#endif
	init_pos();

	#ifndef SPLIT_FILES
	f_fd = ::open(fname_base.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
	if (f_fd == -1) { ::perror("can't create xyc data file"); return 1; }
	#endif
	auto f_write_reset = [&] () {
		pt2_t xNaN = { NaN, NaN };
		::write(f_fd, &xNaN, 2*sizeof(double));
	};
	
	while (continue_running) {
		
		if (step%1000000 == 0) {
			fmt::print("t={:.3f}, {} resets\n", t, resets);
			#ifdef SPLIT_FILES
			if (f_fd != -1)
				::close(f_fd);
			f_fd = ::open(fmt::format(fname_format,f_num).c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
			if (f_fd == -1) { ::perror("can't create xyc data file"); return 1; }
			f_num++;
			#endif
			#ifndef RESET_WITH_TRAPPING
			if (step == 0)  // init ≡ first reset
				f_write_reset();
			#endif
		}
		
		// resetting
		#ifdef RESET_WITH_TRAPPING
		if (trapping) {
			if (t > t_trap_end) {
				trapping = false;
				reset_init();
			}
		} else
		#endif
		{
			if (reset_do()) {
				#ifndef RESET_WITH_TRAPPING
				f_write_reset();
				reset_init();
				#endif
				init_pos();
				resets++;
			}
		}
		
		vec2_t f = sqrt(2 * γ * T / Δt) * vec2_t{ .x = normal_distrib_gen(rng), .y = normal_distrib_gen(rng) };
		#ifdef RESET_WITH_TRAPPING
		if (trapping) {
			f += k_trap * (x0_trap - x);
		}
		#endif
		#ifndef LANGEVIN_OVERDAMPED
		vec2_t a = -γ*v + f;
		v += a / part_m * Δt;
		#else
		vec2_t v = f / γ;
		#endif
		x = x + v * Δt;
		
		::write(f_fd, &x, 2*sizeof(double));
		#ifdef RESET_WITH_TRAPPING
		::write(f_fd, &trapping, 1);
		#endif
		
		t += Δt;
		step++;
	}
	
	fmt::print("\nt_end={:.3e}, steps={}, time step Δt={}, temp T={}, friction γ={} => diffusion D={}\n", t, step, Δt, T, γ, T/γ);
	#ifdef RESET_WITH_TRAPPING
	fmt::print("reset with trapping, σ={:.4f}, T_trap={}, centered at ({},{})\n", sqrt(T/k_trap), t_trap, x0_trap.x, x0_trap.y);
	#else
	fmt::print("ideal reset, σ={}\n", init_pos_sigma);
	#endif
	#ifdef ENABLE_PERIODICAL_RESET
	fmt::print("{} resets (periodical type with T_res={})\n", resets, reset_period);
	#endif
	#ifdef ENABLE_POISSON_RESET
	fmt::print("{} resets (poissonian type with α={})", resets, proba_reset_step/Δt);
	#endif
	fmt::print("--------------------------\n\n");
	
	// For use in `exp-data-diffus-analysis.ipynb`
	#if defined(ENABLE_PERIODICAL_RESET) && defined(SPLIT_FILES) && defined(RESET_WITH_TRAPPING)
	fmt::print("name = '{}'\n", fname_base);
	fmt::print("N = (1,{})\n", f_num-1);
	fmt::print("fps = {:.2f}\n", 1/Δt);
	fmt::print("std_calib = sqrt({})\n", T/k_trap);
	fmt::print("creneau_inv = 1\n");
	fmt::print("fit_diffus_t_end = {}*0.95\n", reset_period);
	fmt::print("D_calib = {}\n", T/γ);
	#ifndef LANGEVIN_OVERDAMPED
	fmt::print("part_m = {}\n", part_m);
	#endif
	#endif
	
	// For use directly in langevin-ft-automated.ipynb
	#if !defined(RESET_WITH_TRAPPING) && !defined(SPLIT_FILES)
	fmt::print("N_traj,{}\n", resets);
	fmt::print("sigma_x,{}\n", init_pos_sigma);
	fmt::print("sigma_y,{}\n", init_pos_sigma);
	fmt::print("D,{}\n", T/γ);
	fmt::print("D_err,{}\n", 0);
	fmt::print("fps,{}\n", 1/Δt);
	fmt::print("reset_period,{}\n", reset_period);
	#endif

	::close(f_fd);

	return 0;
}
