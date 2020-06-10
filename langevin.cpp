#include <fmt/core.h>
#include "pysimul-common.h"
#include <random>

constexpr uint16_t N_gas = 0;
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
	double T = 10;
	_register_var(_thread, "T", &T);
	
	std::random_device _rd;
	std::mt19937 rng (_rd());
	std::normal_distribution<> normal_distrib_gen (0, 1); // mean, std
	std::uniform_real_distribution<> unif01 (0, 1);
	
	#ifdef ENABLE_HARMONIC_WELL
	double well_k = 0e4;
	pt2_t well_center = {0.5, 0.5};
	_register_var(_thread, "well_k", &well_k); _register_var(_thread, "well_x", &well_center.x); _register_var(_thread, "well_y", &well_center.y);
	#endif
	
	// Stats
	#ifdef ENABLE_HIST
	std::vector<double> s_t;
	_register_var(_thread, "sample_t", &s_t);
	// Particle position
	std::vector<double> part_x, part_y;
	_register_var(_thread, "part_x", &part_x); _register_var(_thread, "part_y", &part_y);
	// Particle speed
	std::vector<double> part_vx, part_vy;
	_register_var(_thread, "part_vx", &part_vx); _register_var(_thread, "part_vy", &part_vy);
	#endif
	double part_T_acc = 0; uint64_t part_T_samples = 0;
	_register_var(_thread, "part_T", &part_T_acc); _register_var(_thread, "part_T_samples", &part_T_samples);
	
	// Particule position distribution
	constexpr size_t xdist_Ndr = 200;
	constexpr double xdist_max = 0.1;
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
	constexpr size_t display_period = 200000;
	
	constexpr double Δt = 50 * 1.5e-6;
	_register_const(_thread, "Delta_t", Δt);
	uint8_t pause = 0;
	_register_var(_thread, "pause", &pause);
	double t_pause = Inf;
	_register_var(_thread, "t_pause", &t_pause);
	
	// Poissonian resetting
	double reset_rate = 0.0002 / Δt;
	_register_var(_thread, "reset_rate", &reset_rate);
	const double proba_reset_step = Δt * reset_rate;
	
	pt2_t x = pt2_t{0.5,0.5};
	vec2_t v = {0,0};
	
	while (not _thread.do_quit) {
		
		if (step%display_period == 0) {
			#ifndef SIMUL_HEADLESS
			while (not _thread.win_evts.empty()) {
				sf::Event& event = _thread.win_evts.front();
				if (event.type == sf::Event::Closed)
					_thread.do_quit = true;
				if (event.type == sf::Event::KeyPressed) {
					switch (event.key.code) {
						case sf::Keyboard::Q: _thread.do_quit = true; break;
						case sf::Keyboard::P: pause = !pause; break;
						default: break;
					}
				}
				_thread.win_evts.pop();
			}
			#endif
			if (_thread.regular_callback)
				_thread.regular_callback(_thread.id_for_callback, step, t);
		}
		if (pause) {
			::pthread_mutex_unlock(&_thread.mutex_global);
			usleep(100000);
			::pthread_mutex_lock(&_thread.mutex_global);
		} else {
		
			/******************************/
			
			#ifdef ENABLE_HIST
			s_t.push_back(t);
			part_x.push_back( x.x );
			part_y.push_back( x.y );
			part_vx.push_back( v.x );
			part_vy.push_back( v.y );
			#endif
			
			vec2_t f_alea = sqrt(2 * γ * T / Δt) * vec2_t{ .x = normal_distrib_gen(rng), .y = normal_distrib_gen(rng) };
			vec2_t f_ext = {0,0};
			#ifdef ENABLE_HARMONIC_WELL
			f_ext += well_k * (well_center - x);
			#endif
			vec2_t a = -γ*v + f_alea + f_ext;
			
			v += a / part_m * Δt;
			x = x + v * Δt;
			
			constexpr double kB = 1;
			part_T_acc += part_m/2 * !(v) / kB; //  Ecin / (2 DoF * 1/2 * kB)
			part_T_samples++;
			
			constexpr double dx = xdist_max/xdist_Ndr;
			vec2_t rpos = x - pt2_t{0.5,0.5};
			int64_t x_k = ::lround( rpos.x / dx );
			if (std::abs(x_k) <= xdist_Ndr)
				xdist_acc[ xdist_Ndr + x_k ]++;
			uint64_t x_r = ::lround( rpos.r() / dx );
			if (x_r <= xdist_Ndr)
				rdist_acc[ x_r ]++;
			xdist_samples++;
			
			t += Δt;
			step++;
			
			if (unif01(rng) < proba_reset_step) {
				x = pt2_t{0.5,0.5};
			//	v = {0,0};
			}
			
			if (t > t_pause) {
				pause = true;
				step = 0;
			}
			
			/******************************/
			
		}
		if (step%display_period == 0) {
			::pthread_mutex_unlock(&_thread.mutex_global);
			
			#ifndef SIMUL_HEADLESS
			win.clear(sf::Color::White);
			
			sf::CircleShape circle = sf::c01::buildCircleShapeCR(x, 4e-2);
			circle.setFillColor(sf::Color(100));
			win.draw(circle);
			
			auto text = sf::c01::buildText(font, pt2_t{0.01,1-0.01}, {
				fmt::format(L"t={:.2e}", t),
				fmt::format(L"part_T={:.2e}", part_T_acc/part_T_samples),
			});
			win.draw(text);
			
			sf::VertexArray line (sf::Lines, 2);
			for (size_t i = 0; i < 2*xdist_Ndr+1; i++) {
				line[0] = sf::Vertex( sf::Vector2f(i,SFMLC01_WINDOW_HEIGHT), sf::Color::Blue );
				line[1] = sf::Vertex( sf::Vector2f(i,SFMLC01_WINDOW_HEIGHT-xdist_acc[i]/(float)xdist_samples*5e4), sf::Color::Blue );
				win.draw(line);
			}
			
			win.display();
			#else
			usleep(10);
			#endif
			::pthread_mutex_lock(&_thread.mutex_global);
		}
	}
	
	::pthread_mutex_unlock(&_thread.mutex_global);
	return nullptr;
}
