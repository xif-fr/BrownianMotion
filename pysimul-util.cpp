#include <algorithm>

extern "C" {
	
	#include "cffi-proto.h"

	void util_cma (const double* X, double* M, size_t n, uint16_t window) {
		for (size_t i = 0; i < n; i++) {
			ssize_t a = i - window;
			if (a < 0) a = 0;
			size_t b = i + window;
			if (b > n-1) b = n-1;
			for (size_t j = a; j <= b; j++) {
				M[i] += X[j];
			}
			M[i] /= b-a+1;
		}
	}
	
}