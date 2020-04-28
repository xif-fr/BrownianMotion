#include <cstddef>
#include <array>
#include <utility>
#include <vector>
#include <functional>

#ifndef __XIFUTILS_VEC2__
#define __XIFUTILS_VEC2__

struct vec2_t {
	double x, y;
	double& operator[] (size_t i)       { return ((double*)this)[i]; }
	void    operator/= (double λ)       { x /= λ; y /= λ; }
	void    operator*= (double λ)       { x *= λ; y *= λ; }
	void    operator+= (vec2_t o)       { x += o.x; y += o.y; }
	void    operator-= (vec2_t o)       { x -= o.x; y -= o.y; }
	vec2_t  operator*  (double λ) const { return vec2_t{ λ*x, λ*y }; }
	vec2_t  operator/  (double λ) const { return vec2_t{ x/λ, y/λ }; }
	vec2_t  operator+  (vec2_t o) const { return vec2_t{ x+o.x, y+o.y }; }
	vec2_t  operator-  ()         const { return vec2_t{ -x, -y }; }
	vec2_t  operator-  (vec2_t o) const { return vec2_t{ x-o.x, y-o.y }; }
	double  operator|  (vec2_t o) const { return x*o.x + y*o.y; }
	double  operator!  ()         const { return x*x + y*y; }
	double  r          ()         const;
};
inline vec2_t operator* (double λ, const vec2_t& v) { return v*λ; }
constexpr vec2_t O⃗ = {0.,0.};

struct pt2_t {
	double x, y;
	vec2_t operator- (pt2_t o) const { return vec2_t{ x-o.x, y-o.y }; }
	pt2_t operator+ (vec2_t v) const { return pt2_t{ x+v.x, y+v.y }; }
};

struct vecO_t {
	double r, θ;
	vecO_t operator* (double λ) const { vecO_t o; o.r *= λ; return o; }
	vecO_t operator/ (double λ) const { vecO_t o; o.r /= λ; return o; }
	void operator*= (double λ) { r *= λ; }
	void operator/= (double λ) { r /= λ; }
	operator vec2_t() const;
};

#endif
