#include <utility>
#include <math.h>
double rand01 ();
#include <limits>
#define π M_PI
#define Inf std::numeric_limits<double>::infinity()
#define NaN std::numeric_limits<double>::signaling_NaN()
#include <inttypes.h>
#include "vecs2.hpp"
#include <vector>
#include <map>
#include <pthread.h>
#include <queue>
#include <unistd.h>

#ifndef SIMUL_HEADLESS
#include <SFML/Graphics.hpp>
#include "sfml_c01.hpp"
#endif

extern const size_t pysimul_N;

struct simul_thread_info_t {
	enum var_type_t : int8_t {
		VAR_DOUBLE = 0,
		VAR_SERIES = -1, // std::vector<double>
		VAR_N_ARRAY = -2, // std::array<double,pysimul_N>
		VAR_DISTRIB_ACC = -3, // std::array<number_t,N> + uint64_t
		// 1, 2, 4, 8 -> sizeof integer
	};
	std::map<std::string, std::pair<var_type_t,void*>> vars;
	std::map<std::string, double> vars_pre_init;
	pthread_mutex_t mutex_global = PTHREAD_MUTEX_INITIALIZER;
	pthread_t thread_id;
	#ifndef SIMUL_HEADLESS
	sf::RenderWindow* win;
	std::queue<sf::Event> win_evts;
	#endif
	uint64_t id_for_callback;
	std::function<void(uint64_t id, size_t step, double t)> regular_callback;
	bool do_quit = false;
};

template <typename number_t> void _register_var (simul_thread_info_t& _thread, const char* key, number_t* var) {
	std::map<std::string,double>::iterator it;
	if ((it = _thread.vars_pre_init.find(key)) != _thread.vars_pre_init.end())
		*var = (number_t)::lround(it->second);
	_thread.vars.insert({key, {(simul_thread_info_t::var_type_t)sizeof(number_t), var}});
}
template<> void _register_var<double> (simul_thread_info_t& _thread, const char* key, double* var);
template<> void _register_var<std::vector<double>> (simul_thread_info_t& _thread, const char* key, std::vector<double>* var);

template <size_t N> void _register_Narray (simul_thread_info_t& _thread, const char* key, std::array<double,N>& array) {
	_thread.vars.insert({key, { simul_thread_info_t::VAR_N_ARRAY, array.data() }});
}

template <size_t N, typename number_t> void _register_distrib (simul_thread_info_t& _thread, const char* key, std::array<number_t,N>& dist_acc, uint64_t* dist_samples) {
	_thread.vars.insert({key, {
		simul_thread_info_t::VAR_DISTRIB_ACC,
		new std::tuple<size_t,simul_thread_info_t::var_type_t,void*,uint64_t*>( N, (simul_thread_info_t::var_type_t)sizeof(number_t), dist_acc.data(), dist_samples ) // memory leak
	}});
}
template <size_t N> void _register_distrib (simul_thread_info_t& _thread, const char* key, std::array<double,N>& dist_acc, uint64_t* dist_samples) {
	_thread.vars.insert({key, {
		simul_thread_info_t::VAR_DISTRIB_ACC,
		new std::tuple<size_t,simul_thread_info_t::var_type_t,void*,uint64_t*>( N, simul_thread_info_t::VAR_DOUBLE, dist_acc.data(), dist_samples ) // memory leak
	}});
}

template <typename number_t> void _register_const (simul_thread_info_t& _thread, const char* key, number_t val) {
	_thread.vars.insert({key, {(simul_thread_info_t::var_type_t)sizeof(number_t), new number_t(val)}}); // memory leak
}
template<> void _register_const<double> (simul_thread_info_t& _thread, const char* key, double val);

#ifndef SIMUL_HEADLESS
void sfml_create_window (simul_thread_info_t* _thread);
void sfml_event_poll (simul_thread_info_t* _thread);
#endif
void* comp_thread (void* _data);