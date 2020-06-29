
///////////////////////////////////////////////////////////////////////////////
// gfximage.hpp
//
// Data structures for colors and raster images.
//
// This file builds upon gfxnumeric.hpp, so you may want to familiarize
// yourself with that header before diving into this one.
//
// Students: this header is complete as provided; you do not need to edit this
// file.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <vector>

#include "gfxnumeric.hpp" // for approx_equal

namespace gfx {

// Data type for one R, G, or B intensity value in a high dynamic range (HDR)
// color.
using hdr_intensity = float;

// Return true iff x is a valid intensity in the range [0.0, 1.0].
constexpr bool is_hdr_intensity_valid(hdr_intensity x) {
  return (x >= 0.0f) && (x <= 1.0f);
}

// Convert an 8bpp intensity in [0, 255] to an HDR intensity in [0.0, 1.0].
constexpr hdr_intensity byte_to_hdr(uint_least8_t byte) {
  assert(byte >= 0);
  assert(byte <= 255);
  hdr_intensity converted = float(byte) / 255.0f;
  assert(is_hdr_intensity_valid(converted));
  return converted;
}

// Convert an HDR intensity in [0.0, 1.0] to an 8bpp intensity in [0, 255].
constexpr uint_least8_t hdr_to_byte(hdr_intensity x) {
  assert(is_hdr_intensity_valid(x));
  unsigned converted = unsigned(x * 255.0f);
  assert(converted <= 255);
  return uint_least8_t(converted);
}

// Approximate-equality test of an hdr_intensity; this is just a special
// case of gfx::approx_equal declared in gfxnumeric.hpp.
constexpr bool hdr_intensity_approx_equal(hdr_intensity a,
                                          hdr_intensity b,
                                          hdr_intensity epsilon) {
  return approx_equal(a, b, epsilon);
}

// Class defining an HDR RGB color, comprised of a red, green, and blue
// intensity.
class hdr_rgb {
public:
  // Data type used internally to store the three intensities.
  using intensity_array = std::array<hdr_intensity, 3>;

  // Const iterator into an intensity_array.
  using iterator = intensity_array::const_iterator;

private:
  intensity_array intensities_;

public:

  // Construct from separate R, G, B values, each of which must be valid.
  constexpr hdr_rgb(hdr_intensity r, hdr_intensity g, hdr_intensity b)
  : intensities_{{r, g, b}} {
    assert(is_hdr_intensity_valid(r));
    assert(is_hdr_intensity_valid(g));
    assert(is_hdr_intensity_valid(b));
  }

  // Construct an all-zero (black) color.
  constexpr hdr_rgb()
  : hdr_rgb(0.0, 0.0, 0.0) { }

  // Accessors and mutators.
  constexpr hdr_intensity r() const { return intensities_[0]; }
  constexpr hdr_intensity g() const { return intensities_[1]; }
  constexpr hdr_intensity b() const { return intensities_[2]; }
  constexpr void r(hdr_intensity x) {
    assert(is_hdr_intensity_valid(x));
    intensities_[0] = x;
  }
  constexpr void g(hdr_intensity x) {
    assert(is_hdr_intensity_valid(x));
    intensities_[1] = x;
  }
  constexpr void b(hdr_intensity x) {
    assert(is_hdr_intensity_valid(x));
    intensities_[2] = x;
  }

  // Strict equality test. Two colors are identical when their corresponding
  // R, G, B values are equal.
  bool operator==(const hdr_rgb& rhs) const {
    return std::equal(intensities_.begin(), intensities_.end(), rhs.intensities_.begin());
  }

  // Approximate equality test.
  constexpr bool approx_equal(const hdr_rgb& other, hdr_intensity epsilon) const {
    return (hdr_intensity_approx_equal(r(), other.r(), epsilon) &&
            hdr_intensity_approx_equal(g(), other.g(), epsilon) &&
            hdr_intensity_approx_equal(b(), other.b(), epsilon));
  }

  // Reassign R, G, B to separate values all at once. Each intensity value must be valid.
  constexpr void assign(hdr_intensity new_r, hdr_intensity new_g, hdr_intensity new_b) {
    r(new_r);
    g(new_g);
    b(new_b);
  }

  // Const iterators through the R, G, B intensities.
  iterator begin() const { return intensities_.begin(); }
  iterator end()   const { return intensities_.end  (); }

  // Reassign R, G, B all to the same value, which must be valid.
  constexpr void fill(hdr_intensity x) {
    assert(is_hdr_intensity_valid(x));
    intensities_.fill(x);
  }

  // Create an hdr_rgb from 8bpp intensities, each of must be in the range
  // [0, 255].
  static constexpr hdr_rgb from_bytes(uint_least8_t r_byte,
                                      uint_least8_t g_byte,
                                      uint_least8_t b_byte) {
    return hdr_rgb(byte_to_hdr(r_byte),
                   byte_to_hdr(g_byte),
                   byte_to_hdr(b_byte));
  }

  // Create an hdr_rgb from a 24-bit HTML color code.
  static constexpr hdr_rgb from_hex(uint_least32_t web_hex_code) {
    assert(web_hex_code <= 0xFFFFFF);
    auto r_byte = (web_hex_code >> 16) & 0xFF,
         g_byte = (web_hex_code >>  8) & 0xFF,
         b_byte = (web_hex_code      ) & 0xFF;
    return hdr_rgb(byte_to_hdr(r_byte),
                   byte_to_hdr(g_byte),
                   byte_to_hdr(b_byte));
  }

  // Swap contents with another hdr_rgb.
  void swap(hdr_rgb& other) { intensities_.swap(other.intensities_); }
};

// The 16 named web colors, see
// https://en.wikipedia.org/wiki/Web_colors#HTML_color_names
const hdr_rgb WHITE  (hdr_rgb::from_hex(0xFFFFFF)),
              SILVER (hdr_rgb::from_hex(0xC0C0C0)),
              GRAY   (hdr_rgb::from_hex(0x808080)),
              BLACK  (hdr_rgb::from_hex(0x000000)),
              RED    (hdr_rgb::from_hex(0xFF0000)),
              MAROON (hdr_rgb::from_hex(0x800000)),
              YELLOW (hdr_rgb::from_hex(0xFFFF00)),
              OLIVE  (hdr_rgb::from_hex(0x808000)),
              LIME   (hdr_rgb::from_hex(0x00FF00)),
              GREEN  (hdr_rgb::from_hex(0x008000)),
              AQUA   (hdr_rgb::from_hex(0x00FFFF)),
              TEAL   (hdr_rgb::from_hex(0x008080)),
              BLUE   (hdr_rgb::from_hex(0x0000FF)),
              NAVY   (hdr_rgb::from_hex(0x000080)),
              FUSCHIA(hdr_rgb::from_hex(0xFF00FF)),
              PURPLE (hdr_rgb::from_hex(0x800080));

// A 2D raster image; a grid of pixels, where each pixel is an hdr_rgb.
//
// An hdr_image can be in either an empty state, containing no pixels, or
// in a nonempty state with positive width and positive height.
class hdr_image {
private:
    std::vector<std::vector<hdr_rgb>> rows_;

public:

  // Create an empty image.
  hdr_image() {
    assert(is_empty());
  }

  // Create an image with a given width, height, and color for all the pixels.
  // width and height must both be positive.
  hdr_image(size_t width,
            size_t height,
            const hdr_rgb& fill_color)
  : rows_(height, std::vector<hdr_rgb>(width, fill_color)) {
    assert(width > 0);
    assert(height > 0);
    assert(!is_empty());
  }

  // Copy constructor.
  hdr_image(const hdr_image&) = default;

  // Create an image with the same dimensions as other, but all pixels are
  // initialized to fill_color.
  hdr_image(const hdr_image& other,
            const hdr_rgb& fill_color)
  : hdr_image(other.width(), other.height(), fill_color) { }

  // Strict equality comparison. Two images are == when they have identical
  // dimensions, and every pair of corresponding pixels is ==. Two empty
  // images count as ==.
  bool operator==(const hdr_image& rhs) const {
    return (is_same_size(rhs) &&
            std::equal(rows_.begin(),
                       rows_.end(),
                       rhs.rows_.begin(),
                       [](auto& l, auto& r) {
                         return std::equal(l.begin(), l.end(), r.begin());
                       }));
  }

  // Approximate equality. To be approximately equal, both images must have
  // identical dimensions, and every pair of corresponding pixels must be
  // approximately equal according to hdr_rgb::approx_equal. Two empty
  // images count as approx_equal.
  bool approx_equal(const hdr_image& other, hdr_intensity epsilon) const {
    if (!is_same_size(other)) {
      return false;
    }

    for (size_t x = 0; x < width(); ++x) {
      for (size_t y = 0; y < height(); ++y) {
        if (!pixel(x, y).approx_equal(other.pixel(x, y), epsilon)) {
          return false;
        }
      }
    }

    return true;
  }

  // Make the image empty.
  void clear() {
    rows_.clear();
    assert(is_empty());
  }

  // Set every pixel to fill_color.
  void fill(const hdr_rgb& fill_color) {
    for (auto& row : rows_) {
      row.assign(width(), fill_color);
    }
  }

  // Return the height of the image. An empty image has height zero.
  size_t height() const {
    if (is_empty()) {
      return 0;
    } else {
      return rows_.size();
    }
  }

  // Return true iff the image is empty.
  bool is_empty() const { return rows_.empty(); }

  // Return true iff every pixel is == to color.
  bool is_every_pixel(const hdr_rgb& color) const {
    return std::all_of(rows_.begin(),
                       rows_.end(),
                       [&](auto& row) {
                         return std::all_of(row.begin(),
                                            row.end(),
                                            [&](auto& pixel) {
                                              return (pixel == color);
                                            });
                       });
  }

  // Return true when x or y is a valid coordinate for this image. When an
  // image is empty, no coordinate is valid.
  bool is_x(size_t x) const { return x < width();  }
  bool is_y(size_t y) const { return y < height(); }
  bool is_xy(size_t x, size_t y) const {
    return is_x(x) && is_y(y);
  }

  // Return true iff other has identical width and height to this image.
  // Two empty images count as having the same dimensions.
  bool is_same_size(const hdr_image& other) const {
    return (width() == other.width()) && (height() == other.height());
  }

  // Return the pixel color at (x, y).
  // x and y must both be valid coordinates according to is_xy.
  const hdr_rgb& pixel(size_t x, size_t y) const {
    assert(is_xy(x, y));
    return rows_[y][x];
  }

  // Assign the pixel at (x, y) to new_value.
  // x and y must both be valid coordinates according to is_xy.
  void pixel(size_t x, size_t y, const hdr_rgb& new_value) {
    assert(is_xy(x, y));
    rows_[y][x] = new_value;
  }

  // Change dimensions to new_width and new_height.
  //
  // Both new_width and new_height must be positive; if you want to make
  // the image empty, call clear() instead.
  //
  // If the resizing makes the image larger in one or both of the dimensions,
  // the newly-created pixels are initialized to fill_color.
  void resize(size_t new_width,
              size_t new_height,
              const hdr_rgb& fill_color = BLACK) {
    assert(new_width > 0);
    assert(new_height > 0);

    if ((width() == new_width) && (height() == new_height)) {
      return;
    }

    // fix number of rows; this may add empty rows if necessary
    rows_.resize(new_height);

    // fix number of columns
    for (auto& row : rows_) {
      row.resize(new_width, fill_color);
    }
  }

  // Swap contents with another image.
  void swap(hdr_image& other) { rows_.swap(other.rows_); }

  // Return the width of the image. An empty image has width zero.
  size_t width() const {
    if (is_empty()) {
      return 0;
    } else {
      return rows_.front().size();
    }
  }
};

} // namespace gfx
