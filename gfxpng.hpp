
///////////////////////////////////////////////////////////////////////////////
// gfxpng.hpp
//
// Read/write PNG images.
//
// This file builds upon gfximage.hpp, so you may want to familiarize
// yourself with that header before diving into this one.
//
// This code is a thin wrapper over the png++ library at
// https://www.nongnu.org/pngpp/ . In order to compile this code you must have
// png++ installed. The Makefile automatically installs png++ on Tuffix.
//
// Students: this header is complete as provided; you do not need to edit this
// file.
//
///////////////////////////////////////////////////////////////////////////////

#include <optional>
#include <string>

#include <png++/png.hpp>

#include "gfximage.hpp"

namespace gfx {

// Read a PNG file at the given path.
//
// On success, returns a non-empty optional<hdr_image> containing the hdr_image
// with the contents of the image file.
//
// On I/O error, returns an empty optional object.
//
std::optional<hdr_image> read_png(const std::string& path) {
  try {

    png::image<png::rgb_pixel> loaded(path);

    gfx::hdr_image result(loaded.get_width(), loaded.get_height(), BLACK);

    for (size_t y = 0; y < loaded.get_height(); ++y) {
      for (size_t x = 0; x < loaded.get_width(); ++x) {
        auto loaded_pixel = loaded.get_pixel(x, y);
        result.pixel(x, y,
                     hdr_rgb::from_bytes(loaded_pixel.red,
                                         loaded_pixel.green,
                                         loaded_pixel.blue));
      }
    }

    return result;

  } catch (const std::exception& error) {
    return std::optional<hdr_image>();
  }
}

// Write image to a PNG file at the given path.
//
// The given image must be non-empty.
//
// Returns true on success and false on I/O error.
bool write_png(const hdr_image& image, const std::string& path) {
  assert(!image.is_empty());

  try {

    png::image<png::rgb_pixel> truecolor(image.width(), image.height());

    for (size_t y = 0; y < image.height(); ++y) {
      for (size_t x = 0; x < image.width(); ++x) {
        auto& hdr_pixel = image.pixel(x, y);
        png::rgb_pixel byte_pixel(hdr_to_byte(hdr_pixel.r()),
                                  hdr_to_byte(hdr_pixel.g()),
                                  hdr_to_byte(hdr_pixel.b()));
        truecolor.set_pixel(x, y, byte_pixel);
      }
    }

    truecolor.write(path);

    return true;

  } catch (const std::exception& error) {
    return false;
  }
}

// Convenience function: returns true when the PNG images at path1 and path2
// are == . Returns false when either image cannot be loaded, or when the images
// are not ==.
bool png_equal(const std::string& path1, const std::string& path2) {
  auto png1 = read_png(path1),
       png2 = read_png(path2);
  return (png1 && png2 && (*png1 == *png2));
}

} // namsepace gfx
