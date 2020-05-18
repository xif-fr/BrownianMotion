struct simul_thread_info_t;

struct pysimul_init_data_t {
	struct simul_thread_info_t* thread_info;
	size_t N;
};

struct pysimul_init_data_t pysimul_init (void);
void pysimul_start (struct simul_thread_info_t*);
void pysimul_mutex_lock (struct simul_thread_info_t*);
void pysimul_mutex_unlock (struct simul_thread_info_t*);
uint8_t pysimul_event_poll (struct simul_thread_info_t*);
void pysimul_end (struct simul_thread_info_t*); // a thread can be re-used after end() with start()
void pysimul_finish (struct simul_thread_info_t*);

struct pysimul_getvar_t {
	int8_t type; // see simul_thread_info_t::var_type_t, -42 if not found
	int8_t opt_type;
	int64_t ival;
	double fval;
	size_t length;
};
struct pysimul_getvar_t pysimul_get_var (struct simul_thread_info_t*, const char* key);
void pysimul_var_preinit (struct simul_thread_info_t* _thread, const char* key, double val);
void pysimul_set_var_float (struct simul_thread_info_t*, const char* key, double val);
void pysimul_set_var_float_in_Narray (struct simul_thread_info_t*, const char* key, size_t i, double val);
void pysimul_reset_series (struct simul_thread_info_t*, const char* key);
void pysimul_set_var_integer (struct simul_thread_info_t*, const char* key, int64_t val);
void pysimul_get_var_array (struct simul_thread_info_t*, const char* key, void* numpy_array_cdata);
uint64_t pysimul_register_regular_callback (struct simul_thread_info_t*, void (*) (uint64_t, size_t, double));

void util_cma (const double* numpy_array_in, double* numpy_array_out, size_t n, uint16_t window);
