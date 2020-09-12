#include <fmt/core.h>
#include <random>
#include <vector>
#include <algorithm>
#include <utility>
#include <math.h>
#define Inf std::numeric_limits<double>::infinity()
#include <inttypes.h>
#include "vecs2.hpp"
#include <unistd.h>
#include <fcntl.h>

//#define LANGEVIN_OVERDAMPED
#define ENABLE_PERIODICAL_RESET
#define RESET_WITH_TRAPPING

#include <signal.h>
volatile bool continue_running = true;
void interupt_handler (int) {
	continue_running = false;
}

int main (int argc, char const* argv[]) {

	std::string fname_base = "langevin-trap-traj-xyc";
	if (argc > 1)
		fname_base = argv[1];
	size_t f_num = 1;
	if (argc > 2)
		f_num = atoi(argv[2]);
	std::string fname_format = fname_base+"{}";
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
	std::mt19937 rng (_rd());
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
	constexpr double init_pos_sigma = 0.1; // gaussian distribution of initial position
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
	#endif
	init_pos();
	
	while (continue_running) {
		
		if (step%1000000 == 0) {
			fmt::print("t={:.3f}, {} resets\n", t, resets);
			if (f_fd != -1)
				::close(f_fd);
			f_fd = ::open(fmt::format(fname_format,f_num).c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
			if (f_fd == -1) { ::perror("can't create xyc data file"); return 1; }
			f_num++;
		}
		
		// resetting
		if (trapping) {
			#ifdef RESET_WITH_TRAPPING
			if (t > t_trap_end) {
				trapping = false;
				reset_init();
			}
			#endif
		} else {
			if (reset_do()) {
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
		::write(f_fd, &trapping, 1);
		
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
	
	#ifdef ENABLE_PERIODICAL_RESET
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
	
	::close(f_fd);

	return 0;
}
