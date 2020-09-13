import numpy as np

from pysimul_ffi import ffi, lib as c

pysimul_regular_callbacks = {}

@ffi.def_extern()
def pysimul_regular_callback_f(id, step, t):
	f = pysimul_regular_callbacks[id]
	if f is not None:
		f(step, t)

class PySimul:

	def __init__(self):
		data = c.pysimul_init()
		self._handle = data.thread_info
		self._expl_locked = False
		self.N = data.N
		self.started = False
	
	def __del__(self):
		c.pysimul_finish(self._handle)
	
	def start(self):
		c.pysimul_start(self._handle)
		self.started = True

	def sfml_event_poll(self):
		return c.pysimul_event_poll(self._handle)

	def set_regular_callback(self, callback):
		id = c.pysimul_register_regular_callback(self._handle, c.pysimul_regular_callback_f)
		def cb (step, t):
			nonlocal self
			nonlocal callback
			self._expl_locked = True
			callback(step, t)
			self._expl_locked = False
		pysimul_regular_callbacks[id] = cb

	def end(self):
		c.pysimul_end(self._handle)
		# thread_info should stay valid and reusable for a start() again

	def explicit_lock(self):
		c.pysimul_mutex_lock(self._handle)
		self._expl_locked = True

	def explicit_unlock(self):
		c.pysimul_mutex_unlock(self._handle)
		self._expl_locked = False

	def __getitem__(self, key):
		if not self._expl_locked:
			c.pysimul_mutex_lock(self._handle)
		struct = c.pysimul_get_var(self._handle, key.encode('ascii'))
		if struct.type == -42:
			c.pysimul_mutex_unlock(self._handle)
			raise TypeError("var '"+key+"' not found")
		elif struct.type == -3: # array<number_t,N> + uint64_t
			dist = np.zeros(struct.length, dtype={
				0 : np.float64, 1 : np.int8, 2 : np.int16, 4 : np.int32, 8 : np.int64
			}[struct.opt_type])
			c.pysimul_get_var_array(self._handle, key.encode('ascii'), ffi.cast("void*", dist.ctypes.data))
			value = (dist, struct.ival)
		elif struct.type == -2: # array<double,N>
			value = np.zeros(self.N)
			c.pysimul_get_var_array(self._handle, key.encode('ascii'), ffi.cast("double*", value.ctypes.data))
		elif struct.type == -1: # vector<double>
			value = np.zeros(struct.length)
			c.pysimul_get_var_array(self._handle, key.encode('ascii'), ffi.cast("double*", value.ctypes.data))
		elif struct.type == 0: # double
			value = struct.fval
		else:
			value = struct.ival
		if not self._expl_locked:
			c.pysimul_mutex_unlock(self._handle)
		return value

	def __setitem__(self, key, value):
		if self.started and not self._expl_locked:
			c.pysimul_mutex_lock(self._handle)
		key = key.encode('ascii')
		preinit = not self.started
		if isinstance(value, int):
			c.pysimul_set_var_integer(self._handle, preinit, key, value)
		elif isinstance(value, float):
			c.pysimul_set_var_float(self._handle, preinit, key, value)
		elif isinstance(value, str):
			c.pysimul_set_var_string(self._handle, preinit, key, value.encode('ascii'))
		else:
			c.pysimul_mutex_unlock(self._handle)
			raise TypeError("not a float, int or string")
		if self.started and not self._expl_locked:
			c.pysimul_mutex_unlock(self._handle)

	def set_Narray_value (self, key, index, value):
		if not self._expl_locked:
			c.pysimul_mutex_lock(self._handle)
		c.pysimul_set_var_float_in_Narray(self._handle, key.encode('ascii'), index, value);
		if not self._expl_locked:
			c.pysimul_mutex_unlock(self._handle)

	def reset_series (self, key):
		if not self._expl_locked:
			c.pysimul_mutex_lock(self._handle)
		c.pysimul_reset_series(self._handle, key.encode('ascii'))
		if not self._expl_locked:
			c.pysimul_mutex_unlock(self._handle)

def cma (X, window):
	M = np.zeros(len(X))
	c.util_cma(ffi.cast("double*", X.ctypes.data), ffi.cast("double*", M.ctypes.data), len(X), window)
	return M
