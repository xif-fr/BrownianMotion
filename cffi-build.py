import os.path
from cffi import FFI
ffibuilder = FFI()

proto = open("cffi-proto.h", 'r').read()
ffibuilder.cdef(proto+'\n extern "Python" void pysimul_regular_callback_f(uint64_t, size_t, double);')
ffibuilder.set_source("pysimul_ffi", proto, extra_objects=[os.path.abspath("./libpysimul.so")])

if __name__ == "__main__":
	ffibuilder.compile(verbose=True)