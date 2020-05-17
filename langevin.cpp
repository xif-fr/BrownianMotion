#include <fmt/core.h>
#include <sys/time.h>
#include <unistd.h>
#include "pysimul-common.h"
#include <random>

constexpr uint16_t N_gas = 800;
const size_t pysimul_N = N_gas;

#ifndef SIMUL_HEADLESS
void sfml_create_window (simul_thread_info_t* _thread) {
	sf::ContextSettings settings;
	settings.antialiasingLevel = 4;
	_thread->win = new sf::RenderWindow (sf::VideoMode(SFMLC01_WINDOW_UNIT,SFMLC01_WINDOW_HEIGHT), "Brownian motion - Langevin equation", sf::Style::Titlebar, settings);
	_thread->win->setVerticalSyncEnabled(false);
	_thread->win->clear(sf::Color::White);
}
void sfml_event_poll (simul_thread_info_t* _thread) {
	sf::Event event;
	while (_thread->win->pollEvent(event))
		_thread->win_evts.push(event);
}
#endif

void* comp_thread (void* _data) {
	simul_thread_info_t& _thread = *(simul_thread_info_t*)_data;

	#ifndef SIMUL_HEADLESS
	sf::RenderWindow& win = *_thread.win;
	sf::Font font;
	if (not font.loadFromFile(FONT_PATH))
		return nullptr;
	win.display();
	#endif

	double part_m = 10;
	_register_var(_thread, "part_m", &part_m);
	
	double γ = 600;
	_register_var(_thread, "gamma", &γ);
	double T = 300;
	_register_var(_thread, "T", &T);
	
	std::mt19937 random_gen (::time(nullptr));
	std::normal_distribution<> normal_distrib_gen (0, 1); // mean, std
	
	double well_k = 1e4;
	pt2_t well_center = {0.5, 0.5};
	_register_var(_thread, "well_k", &well_k); _register_var(_thread, "well_x", &well_center.x); _register_var(_thread, "well_y", &well_center.y);
	
	// Particle position
	std::vector<double> part_x, part_y;
	_register_var(_thread, "part_x", &part_x); _register_var(_thread, "part_y", &part_y);
	// Particule position distribution
	constexpr size_t xdist_Ndr = 200;
	constexpr double xdist_max = 0.04;
	std::array<uint64_t,2*xdist_Ndr+1> xdist_acc;
	std::array<uint64_t,xdist_Ndr+1> rdist_acc;
	xdist_acc.fill(0); rdist_acc.fill(0);
	uint64_t xdist_samples = 0;
	_register_distrib(_thread, "xdist", xdist_acc, &xdist_samples); _register_distrib(_thread, "rdist", rdist_acc, &xdist_samples);
	_register_const(_thread, "xdist_max", xdist_max);
	
	double t = 0;
	_register_var(_thread, "t", &t);
	size_t step = 0;
	_register_var(_thread, "step", &step);
	constexpr size_t display_period = 10;
	
	constexpr double Δt = 1.5e-4, Δt² = Δt*Δt;
	_register_const(_thread, "Delta_t", Δt);

	pt2_t x = pt2_t{0.5,0.5};
	vec2_t v = {1e1,0};
	
	while (not _thread.do_quit) {
		
		if (step%display_period == 0) {
			#ifndef SIMUL_HEADLESS
			while (not _thread.win_evts.empty()) {
				sf::Event& event = _thread.win_evts.front();
				if (event.type == sf::Event::Closed)
					_thread.do_quit = true;
				if (event.type == sf::Event::KeyPressed) {
					switch (event.key.code) {
						case sf::Keyboard::Q:
							_thread.do_quit = true;
							break;
						default:
							break;
					}
				}
				_thread.win_evts.pop();
			}
			#endif
			if (_thread.regular_callback)
				_thread.regular_callback(_thread.id_for_callback, step, t);
		}
		
		vec2_t f_alea = sqrt(2 * γ * T / Δt) * vec2_t{ .x = normal_distrib_gen(random_gen), .y = normal_distrib_gen(random_gen) };
		vec2_t f_ext = well_k * (well_center - x);
		vec2_t a = -γ*v + f_alea + f_ext;
		
		v += a / part_m * Δt;
		x = x + v * Δt;
		
		if (step%display_period == 0) {
			::pthread_mutex_unlock(&_thread.mutex_global);
			
			#ifndef SIMUL_HEADLESS
			win.clear(sf::Color::White);
			
			sf::CircleShape circle = sf::c01::buildCircleShapeCR(x, 4e-2);
			circle.setFillColor(sf::Color(100));
			win.draw(circle);
			
			win.display();
			#else
			usleep(10);
			#endif
			::pthread_mutex_lock(&_thread.mutex_global);
		}
		
		t += Δt;
		step++;
	}
	
	::pthread_mutex_unlock(&_thread.mutex_global);
	return nullptr;
}
