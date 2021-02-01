#include <utility>
#include <inttypes.h>
#include <vector>
#include <array>
#include <map>
#include <string>
#include <pthread.h>
#include <queue>
#include <functional>

#ifndef SIMUL_HEADLESS
#include <SFML/Graphics.hpp>
#endif

extern const size_t pysimul_N;

struct simul_thread_info_t {
	enum var_type_t : int8_t {
		VAR_DOUBLE = 0,
		VAR_SERIES = -1, // std::vector<double>
		VAR_N_ARRAY = -2, // std::array<double,pysimul_N>
		VAR_DISTRIB_ACC = -3, // std::array<number_t,N> + uint64_t
		// 1, 2, 4, 8 -> sizeof integer
		VAR_STRING = -4,
	};
	std::map<std::string, std::pair<var_type_t,void*>> vars;
	std::map<std::string, std::pair<var_type_t,void*>> vars_pre_init;
	pthread_mutex_t ticket_mutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_cond_t ticket_cond = PTHREAD_COND_INITIALIZER; // see https://stackoverflow.com/questions/12685112/pthreads-thread-starvation-caused-by-quick-re-locking
	unsigned long ticket_queue_head = 0, ticket_queue_tail = 0;
	pthread_t thread_id;
	#ifndef SIMUL_HEADLESS
	sf::RenderWindow* win;
	std::queue<sf::Event> win_evts;
	#endif
	uint64_t id_for_callback;
	std::function<void(uint64_t id, size_t step, double t)> regular_callback;
	bool do_quit = false;
};

extern "C" {
void pysimul_mutex_lock (simul_thread_info_t*);
void pysimul_mutex_unlock (simul_thread_info_t*);
}

struct simul_thread_raii_mutex_unlock {
	simul_thread_info_t* thread_info;
	simul_thread_raii_mutex_unlock(simul_thread_info_t* thread_info) : thread_info(thread_info) {}
	~simul_thread_raii_mutex_unlock() { pysimul_mutex_unlock(thread_info); }
};

template <typename number_t> void _register_var (simul_thread_info_t& _thread, const char* key, number_t* var) {
	decltype(simul_thread_info_t::vars_pre_init)::iterator it;
	if ((it = _thread.vars_pre_init.find(key)) != _thread.vars_pre_init.end()) {
		if (it->second.first == 8) 
			*var = (number_t)(*(int64_t*)(it->second.second));
		else if (it->second.first == simul_thread_info_t::VAR_DOUBLE)
			*var = (number_t)(*(double*)(it->second.second));
		// todo : delete and remove entry
	}
	_thread.vars.insert({key, {(simul_thread_info_t::var_type_t)sizeof(number_t), var}});
}
template<> void _register_var<double> (simul_thread_info_t& _thread, const char* key, double* var);
template<> void _register_var<float> (simul_thread_info_t& _thread, const char* key, float* var);
template<> void _register_var<std::vector<double>> (simul_thread_info_t& _thread, const char* key, std::vector<double>* var);
template<> void _register_var<std::string> (simul_thread_info_t& _thread, const char* key, std::string* var);

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
