#include <xifutils/maths/vecs2.hpp>
#include <cmath>
#include <cstdlib>

double vec2_t::r () const {
	return hypotf(x, y);
}

vecO_t::operator vec2_t() const {
	return vec2_t({
		.x = r*cos(θ),
		.y = r*sin(θ)
	});
}
