#include "pysimul-common.h"
#include <string.h>

double rand01 () {
	return ::rand()/(double)RAND_MAX;
}

template<> void _register_var<double> (simul_thread_info_t& _thread, const char* key, double* var) {
	std::map<std::string,double>::iterator it;
	if ((it = _thread.vars_pre_init.find(key)) != _thread.vars_pre_init.end())
		*var = it->second;
	_thread.vars.insert({key, {simul_thread_info_t::VAR_DOUBLE, var}});
}

template<> void _register_var<std::vector<double>> (simul_thread_info_t& _thread, const char* key, std::vector<double>* var) {
	_thread.vars.insert({key, {simul_thread_info_t::VAR_SERIES, var}});
}

template<> void _register_const<double> (simul_thread_info_t& _thread, const char* key, double val) {
	_thread.vars.insert({key, {simul_thread_info_t::VAR_DOUBLE, new double(val)}});
}

#ifdef MAIN_SOLO

int main (int argc, char const* argv[]) {
	simul_thread_info_t* _thread = new simul_thread_info_t;
	sfml_create_window(_thread);
	_thread->regular_callback = [&] (uint64_t, size_t, double) { sfml_event_poll(_thread); };
	comp_thread(_thread);
	return 0;
}

#else

extern "C" {

#include "cffi-proto.h"

pysimul_init_data_t pysimul_init () {
	simul_thread_info_t* _thread = new simul_thread_info_t;
	#ifndef SIMUL_HEADLESS
	sfml_create_window(_thread);
	#endif
	_thread->regular_callback = nullptr;
	_thread->id_for_callback = (uint64_t)_thread;
	return { _thread, pysimul_N };
}

void pysimul_var_preinit (simul_thread_info_t* _thread, const char* key, double val) {
	_thread->vars_pre_init[key] = val;
}

void pysimul_start (simul_thread_info_t* _thread) {
	_thread->do_quit = false;
	#ifndef SIMUL_HEADLESS
	_thread->win_evts = std::queue<sf::Event>();
	#endif
	_thread->vars.clear();
	::pthread_mutex_lock(&(_thread->mutex_global)); // locked for the thread in advance
	pthread_attr_t attr;
	::pthread_attr_init(&attr);
	::pthread_create(&(_thread->thread_id), &attr, &comp_thread, _thread);
	::pthread_attr_destroy(&attr);
}

void pysimul_mutex_lock (simul_thread_info_t* _thread) {
	::pthread_mutex_lock(&(_thread->mutex_global));
}

pysimul_getvar_t pysimul_get_var (simul_thread_info_t* _thread, const char* key) {
	auto it = _thread->vars.find(key);
	if (it == _thread->vars.end())
		return { .type = -42 };
	pysimul_getvar_t r = { .type = it->second.first };
	switch (r.type) {
		case simul_thread_info_t::VAR_DISTRIB_ACC /*-3*/: {
			auto& t = *(std::tuple<size_t,simul_thread_info_t::var_type_t,void*,uint64_t*>*)it->second.second;
			r.length = std::get<0>(t);
			r.opt_type = std::get<1>(t);
			r.ival = *std::get<3>(t);
		} break;
		case simul_thread_info_t::VAR_N_ARRAY /*-2*/: break;
		case simul_thread_info_t::VAR_SERIES /*-1*/: r.length = (*(std::vector<double>*)it->second.second).size(); break;
		case simul_thread_info_t::VAR_DOUBLE /*0*/: r.fval = *(double*)it->second.second; break;
		case 1: r.ival = *(int8_t*)it->second.second; break;
		case 2: r.fval = *(int16_t*)it->second.second; break;
		case 4: r.fval = *(int32_t*)it->second.second; break;
		case 8: r.fval = *(int64_t*)it->second.second; break;
	}
	return r;
}

void pysimul_set_var_float (struct simul_thread_info_t* _thread, const char* key, double val) {
	auto v = _thread->vars.at(key);
	if (v.first == simul_thread_info_t::VAR_DOUBLE)
		*(double*)v.second = val;
	else
		throw std::runtime_error("pysimul_set_var_float : var is not double");
}

void pysimul_set_var_integer (struct simul_thread_info_t* _thread, const char* key, int64_t val) {
	auto v = _thread->vars.at(key);
	switch (v.first) {
		case 0: *(double*)v.second = val; break;
		case 1: *(int8_t*)v.second = (int8_t)val; break;
		case 2: *(int16_t*)v.second = (int16_t)val; break;
		case 4: *(int32_t*)v.second = (int32_t)val; break;
		case 8: *(int64_t*)v.second = val; break;
		default:
			throw std::runtime_error("pysimul_set_var_integer : var is not number");
	}
}

void pysimul_set_var_float_in_Narray (simul_thread_info_t* _thread, const char* key, size_t i, double val) {
	if (i >= pysimul_N)
		throw std::out_of_range("pysimul_set_var_float_in_Narray : out of range");
	auto v = _thread->vars.at(key);
	if (v.first == simul_thread_info_t::VAR_N_ARRAY)
		(*(double**)v.second)[i] = val;
	else
		throw std::runtime_error("pysimul_set_var_float_in_Narray : var is not array");
}

void pysimul_get_var_array (simul_thread_info_t* _thread, const char* key, void* numpy_array_cdata) {
	auto v = _thread->vars.at(key);
	switch (v.first) {
		case simul_thread_info_t::VAR_SERIES: {
			std::vector<double>& array = *(std::vector<double>*)v.second;
			::memcpy(numpy_array_cdata, array.data(), array.size()*sizeof(double));
		} break;
		case simul_thread_info_t::VAR_N_ARRAY: {
			::memcpy(numpy_array_cdata, v.second, pysimul_N*sizeof(double));
		} break;
		case simul_thread_info_t::VAR_DISTRIB_ACC: {
			auto& t = *(std::tuple<size_t,simul_thread_info_t::var_type_t,uint64_t*,uint64_t*>*)v.second;
			uint8_t sz = std::get<1>(t)==0 ? sizeof(double) : std::get<1>(t);
			::memcpy(numpy_array_cdata, std::get<2>(t), std::get<0>(t)*sz);
		} break;
		default:
			throw std::runtime_error("pysimul_get_var_array : var is not array");
	}
}

void pysimul_reset_series (simul_thread_info_t* _thread, const char* key) {
	auto v = _thread->vars.at(key);
	if (v.first == simul_thread_info_t::VAR_SERIES)
		(*(std::vector<double>*)v.second).clear();
	else if (v.first == simul_thread_info_t::VAR_DISTRIB_ACC) {
		auto& t = *(std::tuple<size_t,simul_thread_info_t::var_type_t,void*,uint64_t*>*)v.second;
		for (size_t i = 0; i < std::get<0>(t); i++) {
			switch (std::get<1>(t)) {
				case 0: ((double*)std::get<2>(t))[i] = 0; break;
				case 1: ((int8_t*)std::get<2>(t))[i] = 0; break;
				case 2: ((int16_t*)std::get<2>(t))[i] = 0; break;
				case 4: ((int32_t*)std::get<2>(t))[i] = 0; break;
				case 8: ((int64_t*)std::get<2>(t))[i] = 0; break;
				default:
					throw std::runtime_error("pysimul_reset_series : VAR_DISTRIB_ACC var is not number");
			}
		}
		*std::get<3>(t) = 0;
	} else
		throw std::runtime_error("pysimul_reset_series : var is not series or distrib");
}

void pysimul_mutex_unlock (simul_thread_info_t* _thread) {
	::pthread_mutex_unlock(&(_thread->mutex_global));
}

uint8_t pysimul_event_poll (simul_thread_info_t* _thread) {
	#ifndef SIMUL_HEADLESS
	if (_thread->win) {
		pysimul_mutex_lock(_thread);
		sfml_event_poll(_thread);
		pysimul_mutex_unlock(_thread);
		return 0;
	} else
		return 1;
	#else
	return 42;
	#endif
}

void pysimul_end (simul_thread_info_t* _thread) {
	_thread->do_quit = true;
	::pthread_join(_thread->thread_id, nullptr);
	_thread->vars.clear();
}

void pysimul_finish (simul_thread_info_t* _thread) {
	#ifndef SIMUL_HEADLESS
	_thread->win->close();
	delete _thread->win;
	#endif
	delete _thread;
}

uint64_t pysimul_register_regular_callback (simul_thread_info_t* _thread, void (*cb) (uint64_t, size_t, double)) {
	_thread->regular_callback = cb;
	return (uint64_t)_thread;
}

}

#endif