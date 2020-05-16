#include <fmt/core.h>
#include <deque>
#include <array>
#include <sys/time.h>
#include <unistd.h>
#include "pysimul-common.h"

constexpr uint16_t N_gas = 400;//800;
const size_t pysimul_N = N_gas;

#ifndef SIMUL_HEADLESS
void sfml_create_window (simul_thread_info_t* _thread) {
	sf::ContextSettings settings;
	settings.antialiasingLevel = 4;
	_thread->win = new sf::RenderWindow (sf::VideoMode(SFMLC01_WINDOW_UNIT,SFMLC01_WINDOW_HEIGHT), "Brownian motion - Gas", sf::Style::Titlebar, settings);
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
	
	::srand((int)::time(nullptr));
	
	#ifndef SIMUL_HEADLESS
	sf::RenderWindow& win = *_thread.win;
	sf::Font font;
	if (not font.loadFromFile(FONT_PATH))
		return nullptr;
	win.display();
	#endif
	
	// Potentiel de Lennard-Jones
	
	constexpr double d₀ = 3e-2, d₀² = d₀*d₀, rₐₚₚ² = (4*d₀)*(4*d₀), rlim² = (9*d₀)*(9*d₀);
	constexpr double E₀ = 2;
	
	auto lj_pot = [&] (double r²) -> double {
		double q⁶ = d₀²/r²;
		q⁶ = q⁶*q⁶*q⁶;
		return 4*E₀ * q⁶*(q⁶-1.);
	};
	auto lj_f = [&] (const vec2_t& rd, double r²) -> vec2_t {
		double q⁶ = d₀²/r²; q⁶ = q⁶*q⁶*q⁶;
		if (r² > rₐₚₚ²)
			return rd*(4*E₀*-6.*q⁶/r²);
		else
			return rd*(4*E₀*(12*q⁶-6)*q⁶/r²);
	};
	
	// Potentiel particule brownienne - particule du liquide
	
	constexpr double d_part = 6e-2; // same d₀
	constexpr double E₀_lp = 2;
	constexpr double part_m = 10;//50
	_register_const(_thread, "part_m", part_m);
	_register_const(_thread, "part_d", d_part);
	
	auto ljp_pot = [&] (double r) -> double {
		double δr = r - d_part, δr² = δr*δr;
		double q⁶ = d₀²/δr²; q⁶ = q⁶*q⁶*q⁶;
		return 4*E₀_lp * q⁶*(q⁶-1.);
	};
	auto ljp_f = [&] (const vec2_t& rd, double r) -> vec2_t {
		double δr = r - d_part;
		vec2_t δrd = δr * rd/r;
		double δr² = δr*δr;
		double q⁶ = d₀²/δr²; q⁶ = q⁶*q⁶*q⁶;
		return δrd*(4*E₀_lp*(12*q⁶-6)*q⁶/δr²);
	};
	
	// Potentiel harmonique de confinement de la particule
	
	double well_k = 1e5;
	pt2_t well_center = {0.5, 0.5};
	_register_var(_thread, "well_k", &well_k); _register_var(_thread, "well_x", &well_center.x); _register_var(_thread, "well_y", &well_center.y);
	
	auto well_pot = [&] (pt2_t pos_part) -> double {
		double r² = !(pos_part - well_center);
		return well_k * r² / 2;
	};
	auto well_f = [&] (pt2_t pos_part) -> vec2_t {
		vec2_t rd = well_center - pos_part;
		return well_k * rd;
	};
	
	// Potentiel du conteneur rond
	
	constexpr double cont_r = 0.48;//0.6762;//0.48;
	_register_const(_thread, "cont_r", cont_r);
	constexpr double cont_m = 1e10;
	_register_const(_thread, "cont_m", cont_m);
	constexpr double cont_k = 0.01;
	_register_const(_thread, "cont_k", cont_k);
	
	auto cont_pot = [&] (double r) -> double {
		double f = cont_k/(cont_r-r);
		return cont_k*f*f*f*f*f/5.;
	};
	auto cont_f = [&] (const vec2_t& x, double r) -> vec2_t {
		double f = cont_k/(cont_r-r);
		f = f*f*f; f = f*f;
		vec2_t u = -x/r;
		return f*u;
	};
	
	// Statistiques
	#undef STATS_SPEED_DISTRIB
	#define STATS_FORCE_AUTOCORRELATION
	#define STATS_POSITION_DISTRIB
	
	constexpr size_t stat_period = 50;
	uint8_t enable_fine_stats = 0;
	_register_var(_thread, "enable_fine_stats", &enable_fine_stats);
	std::vector<double> s_t, s_Ecin, s_Epot, s_T, s_P;
	_register_var(_thread, "sample_t", &s_t);
	_register_var(_thread, "Ecin", &s_Ecin); _register_var(_thread, "Epot", &s_Epot);
	_register_var(_thread, "Temp", &s_T); //_register_var(_thread, "Pvir", &s_P);
	double target_T = 1.5e2;
	_register_var(_thread, "target_T", &target_T);
	constexpr double V = π * cont_r*cont_r; // "volume"
	_register_const(_thread, "V", V);
	// Speed distribution
	#ifdef STATS_SPEED_DISTRIB
	constexpr size_t vdist_Ndv = 200;
	constexpr double vdist_max = 1e2;
	std::array<uint64_t,vdist_Ndv> vdist_acc;
	vdist_acc.fill(0);
	uint64_t vdist_samples = 0;
	_register_distrib(_thread, "vdist", vdist_acc, &vdist_samples);
	_register_const(_thread, "vdist_max", vdist_max);
	#endif
	// Particle position
	std::vector<double> part_x, part_y;
	_register_var(_thread, "part_x", &part_x); _register_var(_thread, "part_y", &part_y);
	vec2_t part_x_acc = O⃗;
	// Particle force autocorrelation function & stat_period*Δt-average of speed and force
	#ifdef STATS_FORCE_AUTOCORRELATION
	vec2_t part_f_acc = O⃗, part_v_acc = O⃗;
	std::vector<double> part_fx, part_fy, part_vx, part_vy;
	_register_var(_thread, "part_fx", &part_fx); _register_var(_thread, "part_fy", &part_fy);
	_register_var(_thread, "part_vx", &part_vx); _register_var(_thread, "part_vy", &part_vy);
	constexpr size_t f_autocor_NΔt = 5000;
	std::array<double,f_autocor_NΔt> f_autocor_xx, f_autocor_xy, f_autocor_yy;
	f_autocor_xx.fill(0); f_autocor_xy.fill(0); f_autocor_yy.fill(0);
	uint64_t f_autocor_samples = 0;
	_register_distrib(_thread, "f_autocor_xx", f_autocor_xx, &f_autocor_samples); _register_distrib(_thread, "f_autocor_xy", f_autocor_xy, &f_autocor_samples); _register_distrib(_thread, "f_autocor_yy", f_autocor_yy, &f_autocor_samples);
	std::deque<vec2_t> f_autocor_hist;
	#endif
	// Particule position distribution
	#ifdef STATS_POSITION_DISTRIB
	constexpr size_t xdist_Ndr = 200;
	constexpr double xdist_max = 0.04;
	std::array<uint64_t,2*xdist_Ndr+1> xdist_acc;
	std::array<uint64_t,xdist_Ndr+1> rdist_acc;
	xdist_acc.fill(0); rdist_acc.fill(0);
	uint64_t xdist_samples = 0;
	_register_distrib(_thread, "xdist", xdist_acc, &xdist_samples); _register_distrib(_thread, "rdist", rdist_acc, &xdist_samples);
	_register_const(_thread, "xdist_max", xdist_max);
	#endif
	
	// Contrôle
	
	double t = 0;
	_register_var(_thread, "t", &t);
	size_t step = 0;
	_register_var(_thread, "step", &step);
	constexpr size_t display_period = 200;
	timeval tv_last;
	::gettimeofday(&tv_last,NULL);
	float step_per_s = 0;
	size_t step_last_mes = 0;
	uint8_t pause = 0;
	_register_var(_thread, "pause", &pause);
	// Particule release from central well
	constexpr double release_well_tolerance_center = 3e-4;
	uint8_t release_well = 0;
	_register_var(_thread, "release_well", &release_well);
	double release_t_well = NaN;
	_register_var(_thread, "release_well_t", &release_t_well);
	
	// Particules, intégration et initialisation
	
	std::array<pt2_t,N_gas+2> x;
	std::array<vec2_t,N_gas+2> v, a, a⁻;
	constexpr size_t i_cont = N_gas;   // x[i_cont] is container position
	constexpr size_t i_part = N_gas+1; // x[i_part] is brownian particule position
	
	constexpr double Δt = 1.5e-6, Δt² = Δt*Δt;//1.5e-5
	_register_const(_thread, "Delta_t", Δt);
	#undef KURAEV
	
	constexpr pt2_t part_x0 = pt2_t{0.5,0.5};
	double v₀ = 20;
	constexpr double m = 1;
	
	auto init_particles = [&] () {
		x[i_part] = part_x0;
		v[i_part] = a[i_part] = a⁻[i_part] = O⃗;
		x[i_cont] = pt2_t{0.5,0.5};
		v[i_cont] = a[i_cont] = a⁻[i_cont] = O⃗;
		vec2_t v_tot = O⃗;
		for (size_t i = 0; i < N_gas; i++) {
			v[i] = (vec2_t)vecO_t{ .r = v₀, .θ = rand01()*2*π };
			v_tot += v[i];
			auto dist2others = [&] () {
				double d2 = Inf;
				for (size_t j = 0; j < i; j++)
					d2 = std::min(d2, !(x[i]-x[j]));
				return d2;
			};
			do {
				x[i] = { 2*cont_r*rand01(), 2*cont_r*rand01() };
			} while (dist2others() < d₀²
					 or !(x[i]-x[i_cont]) > cont_r*cont_r*0.9
					 or !(x[i]-x[i_part]) < (d_part+d₀)*(d_part+d₀));
			a[i] = a⁻[i] = O⃗;
		}
		for (size_t i = 0; i < N_gas; i++)
			v[i] -= v_tot / N_gas;
	};
	init_particles();
	
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
						case sf::Keyboard::R:
							init_particles();
							break;
						case sf::Keyboard::P:
							pause = !pause;
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
		
		if (pause) {
			::pthread_mutex_unlock(&_thread.mutex_global);
			usleep(10000);
			::pthread_mutex_lock(&_thread.mutex_global);
			continue;
		}
		
		// Position integration and velocity prediction
		std::array<pt2_t,N_gas+2> x⁺;
		std::array<vec2_t,N_gas+2> a⁺, v⁺;
		for (uint16_t i = 0; i < N_gas+2; i++) {
			#ifdef KURAEV
			// Kuraev method
			x⁺[i] = x[i] + v[i] * Δt + ( 5*a[i] - a⁻[i] ) /8 * Δt²;
			v⁺[i] = v[i] + ( 3*a[i] - a⁻[i] ) /2 * Δt;
			#else
			// Velocity Verlet method
			x⁺[i] = x[i] + v[i] * Δt + a[i] /2 * Δt²;
			#endif
			a⁺[i] = O⃗;
		}
		
		// Forces computation
		for (uint16_t i = 0; i != N_gas; i++) {
			
			vec2_t ri = x⁺[i] - x⁺[i_cont];
			double r = ri.r();
			vec2_t cf = cont_f (ri, r);
			a⁺[i] += cf;
			a⁺[i_cont] += -cf;
			
			for (uint16_t j = i+1; j != N_gas; j++) {
				vec2_t rij = x⁺[i]-x⁺[j];
				double r² = !rij;
				if (r² < rlim²) {
					vec2_t fij = lj_f (rij, r²);
					a⁺[i] += fij;
					a⁺[j] -= fij;
				}
			}
			
			ri = x⁺[i] - x⁺[i_part];
			r = ri.r();
			vec2_t pf = ljp_f (ri, r);
			a⁺[i] += pf;
			a⁺[i_part] += -pf;
			
			a⁺[i_part] += well_f (x⁺[i_part]);
		}
		
		// Berendsen thermostat
		if (target_T > 0 and not s_T.empty()) {
			double T = s_T.back();
			constexpr double τ = 4e-2;
			for (uint16_t i = 0; i < N_gas; i++)
				a⁺[i] -= v[i] /τ * (T/target_T-1);
		}
		
		// Velocity integration
		for (uint16_t i = 0; i < N_gas+2; i++) {
			if (i == i_cont)
				a⁺[i] /= cont_m;
			else if (i == i_part)
				a⁺[i] /= part_m;
			else
				a⁺[i] /= m;
			#ifdef KURAEV
			// Kuraev method
			v[i] += ( 3*a⁺[i] + 6*a[i] - a⁻[i] ) /8 * Δt;
			#else
			// Velocity Verlet method
			v[i] += ( a⁺[i] + a[i] ) /2 * Δt;
			#endif
		}
		a = a⁺;
		a⁻ = a;
		x = x⁺;
		
		// Statistics
		
		if (enable_fine_stats) {
			#ifdef STATS_FORCE_AUTOCORRELATION
			f_autocor_hist.push_back( part_m * a[i_part] );
			if (f_autocor_hist.size() == f_autocor_NΔt) {
				for (size_t k = 0; k < f_autocor_NΔt; k++) {
					f_autocor_xx[k] += f_autocor_hist[0].x * f_autocor_hist[k].x;
					f_autocor_xy[k] += f_autocor_hist[0].x * f_autocor_hist[k].y;
					f_autocor_yy[k] += f_autocor_hist[0].y * f_autocor_hist[k].y;
				}
				f_autocor_samples++;
				f_autocor_hist.pop_front();
			}
			#endif
			#ifdef STATS_POSITION_DISTRIB
			constexpr double dx = xdist_max/xdist_Ndr;
			vec2_t rpos = x[i_part] - well_center;
			int64_t x_k = ::lround( rpos.x / dx );
			if (std::abs(x_k) <= xdist_Ndr)
				xdist_acc[ xdist_Ndr + x_k ]++;
			uint64_t x_r = ::lround( rpos.r() / dx );
			if (x_r <= xdist_Ndr)
				rdist_acc[ x_r ]++;
			xdist_samples++;
			#endif
		}
		#ifdef STATS_FORCE_AUTOCORRELATION
		part_f_acc += part_m * a[i_part];
		part_v_acc += v[i_part];
		#endif
		part_x_acc += x[i_part]-pt2_t{0,0};
		
		if (step%stat_period == 0) {
			s_t.push_back(t);
			
			#ifdef STATS_FORCE_AUTOCORRELATION
			part_fx.push_back( part_f_acc.x / stat_period );
			part_fy.push_back( part_f_acc.y / stat_period );
			part_f_acc = O⃗;
			part_vx.push_back( part_v_acc.x / stat_period );
			part_vy.push_back( part_v_acc.y / stat_period );
			part_v_acc = O⃗;
			#endif
			
			double Epot = 0, Ecin = 0, W = 0;
			
			for (uint16_t i = 0; i != N_gas; i++) {
				Ecin += m/2 * !(v[i]);
				for (uint16_t j = i+1; j != N_gas; j++) {
					vec2_t rij = x[i] - x[j];
					double r² = !rij;
					if (r² < rlim²) {
						Epot += lj_pot (r²);
						vec2_t fij = lj_f (rij, r²);
						W += rij | fij;
					}
				}
				Epot += cont_pot( (x[i] - x[i_cont]).r() );
				Epot += ljp_pot( (x[i] - x[i_part]).r() );
				Epot += well_pot( x[i_part] );
			}
			Ecin += cont_m/2 * !(v[i_cont]);
			Ecin += part_m/2 * !(v[i_part]);
			
			constexpr double kB = 1;
			double T = 2 * Ecin / (kB * 2*(N_gas+2)); // température par équirépartition (Ecin = Ndof 1/2 kB T, avec Ndof = 2 (N_gas+2))
			double P = (2*Ecin - W)/(2*V); // pression par théorème du viriel (est-ce correct ????)
			
			s_Ecin.push_back(Ecin);
			s_Epot.push_back(Epot);
			s_P.push_back(P);
			s_T.push_back(T);
			
			#ifdef STATS_SPEED_DISTRIB
			constexpr double dv = vdist_max/vdist_Ndv;
			for (size_t i = 0; i < N_gas; i++) {
				int64_t k = ::lround( v[i].r() / dv );
				if (k < vdist_Ndv)
					vdist_acc[k]++;
				vdist_samples++;
			}
			#endif
			
			// part_x.push_back(x[i_part].x);	ne pas sous-éch la position, néfaste pour la PSD (repliement)
			// part_y.push_back(x[i_part].y);
			part_x.push_back( part_x_acc.x / stat_period );
			part_y.push_back( part_x_acc.y / stat_period );
			part_x_acc = O⃗;
		}
		
		// Particule release from central well
		if (release_well) {
			double r_from_center = (x[i_part]-well_center).r();
			if (r_from_center < release_well_tolerance_center) {
				release_well = 0;
				release_t_well = t + Δt;
				well_k = 0.;
			}
		}
		
		/********************************************/
		
		if (step%display_period == 0) {
			::pthread_mutex_unlock(&_thread.mutex_global);
			
			timeval tv;
			::gettimeofday(&tv,NULL);
			useconds_t Δµs = (useconds_t)(tv.tv_sec-tv_last.tv_sec)*1000000 + (tv.tv_usec-tv_last.tv_usec);
			if (Δµs > 1e6) {
				step_per_s = (step-step_last_mes) / (Δµs/1e6);
				tv_last = tv;
				step_last_mes = step;
			}
			
			#ifndef SIMUL_HEADLESS
			win.clear(sf::Color::White);
			
			sf::CircleShape cont = sf::c01::buildCircleShapeCR(x[N_gas], cont_r);
			cont.setPointCount(50);
			cont.setOutlineColor(sf::Color::Black);
			cont.setOutlineThickness(1);
			win.draw(cont);
			
			for (uint16_t i = 0; i < N_gas; i++) {
				sf::CircleShape circle = sf::c01::buildCircleShapeCR(x[i], d₀);
				circle.setFillColor(sf::Color(30));
				win.draw(circle);
				circle = sf::c01::buildCircleShapeCR(x[i], d₀/2);
				circle.setFillColor(sf::Color(100));
				win.draw(circle);
			}
			
			sf::CircleShape circle = sf::c01::buildCircleShapeCR(x[i_part], d₀+d_part);
			circle.setFillColor(sf::Color(30));
			win.draw(circle);
			circle = sf::c01::buildCircleShapeCR(x[i_part], d₀/2+d_part);
			circle.setFillColor(sf::Color(100));
			win.draw(circle);
			
			auto text = sf::c01::buildText(font, pt2_t{0.01,1-0.01}, {
				fmt::format(L"Etot={:.5e}, Ecin={:.2e}, Epot={:.2e}", s_Epot.back()+s_Ecin.back(), s_Ecin.back(), s_Epot.back()),
				fmt::format(L"t={:.2e}, sps={:.1e}", t, step_per_s),
				fmt::format(L"T={:.2e}, P={:.2e}", s_T.back(), s_P.back()),
			});
			win.draw(text);
			
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
