
///////////////////////////////////////////////////////////////////////////////
// gfxnumeric.hpp
//
// Numeric operations and constants for graphics programming.
//
// Students: this header is complete as provided; you do not need to edit this
// file.
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <cassert>
#include <cmath>
#include <limits>

namespace gfx {

// Essential numerical constants, in double and float format.

const double DOUBLE_INFINITY = std::numeric_limits<double>::infinity(),
             DOUBLE_NEGATIVE_INFINITY = -DOUBLE_INFINITY,
             DOUBLE_NAN = std::numeric_limits<double>::quiet_NaN(),
             DOUBLE_PI = std::acos(-1.0);

const float FLOAT_INFINITY = std::numeric_limits<float>::infinity(),
            FLOAT_NEGATIVE_INFINITY = -FLOAT_INFINITY,
            FLOAT_NAN = std::numeric_limits<float>::quiet_NaN(),
            FLOAT_PI = std::acos(-1.0f);

// Convert degrees to radians.
double degrees_to_radians(int degrees) {
  return double(degrees) * DOUBLE_PI / 180.0;
}

// Approximate equality testing, within an additive difference (delta),
// robust to infinite and NaN values.
//
// number_type must be a floating point type, most likely double or float.
//
// delta must be finite and positive.
//
// Returns true only when:
//   - a is finite (note NaN, positive infinity, or negative infinity)
//   - b is finite
//   - |a-b| <= delta
template <typename number_type>
constexpr bool approx_equal(number_type a,
                            number_type b,
                            number_type delta) {
  static_assert(!std::numeric_limits<number_type>::is_integer,
                "approx_equal is only defined for floating point types");

  assert(std::isfinite(delta));
  assert(delta > 0);

  return (std::isfinite(a) &&
          std::isfinite(b) &&
          (std::abs(a - b) <= delta));
}

// Clamp: ensure that x is in the range [min, max].
//
// number_type must be a numeric type with < and <= operators.
//
// min must be less than or equal to max.
template <typename number_type>
constexpr number_type clamp(number_type min,
                            number_type x,
                            number_type max) {
  assert(min <= max);
  if (x < min) {
    return min;
  } else if (max < x) {
    return max;
  } else {
    return x;
  }
}

} // namespace gfx
