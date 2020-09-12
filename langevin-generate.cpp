#include <fmt/core.h>
#include <random>
#include <vector>
#include <algorithm>

#define LANGEVIN_OVERDAMPED
#define ENABLE_PERIODICAL_RESET
#define RESET_WITH_TRAPPING

#include <fcntl.h>

void int main (int argc, char const* argv[]) {

	double t = 0;
	size_t step = 0;
	
	constexpr double Δt = 1e-7;
	_register_const(_thread, "Delta_t", Δt);
	
	#ifndef LANGEVIN_OVERDAMPED
	double part_m = 1;
	_register_var(_thread, "part_m", &part_m);
	#endif
	double γ = 1;
	double T = 1;
	
	std::random_device _rd;
	std::mt19937 rng (_rd());
	std::normal_distribution<> normal_distrib_gen (0, 1); // mean, std
	std::uniform_real_distribution<> unif01 (0, 1);
		
	// Poissonian resetting
	#ifdef ENABLE_POISSON_RESET
	const double proba_reset_step = 0.001; // = Δt * reset_rate;
	#endif
	
	// Periodical resetting
	#ifdef ENABLE_PERIODICAL_RESET
	double reset_period = 1;
	#endif



	step = 0;
	t = 0;
	
	pt2_t x;
	#ifndef LANGEVIN_OVERDAMPED
	vec2_t v;
	#endif

	#ifdef RESET_WITH_TRAPPING

	#else
	double init_pos_sigma = 0.1; // gaussian distribution of initial position
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
	#endif
	
	#ifdef ENABLE_PERIODICAL_RESET
	uint32_t n_period = 1;
	#endif
	
	while (true) {
		
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
		
		vec2_t f_alea = sqrt(2 * γ * T / Δt) * vec2_t{ .x = normal_distrib_gen(rng), .y = normal_distrib_gen(rng) };
		#ifndef LANGEVIN_OVERDAMPED
		vec2_t a = -γ*v + f_alea;
		#else
		vec2_t v = f_alea / γ;
		#endif
		#ifndef LANGEVIN_OVERDAMPED
		v += a / part_m * Δt;
		#endif
		x = x + v * Δt;
		
		t += Δt;
		step++;
	}

	return 0;
}
