
#include <iostream>
#include <string>

#include "gfxpng.hpp"
#include "gfxraytrace.hpp"

void print_usage() noexcept;

std::optional<gfx::scene> read_and_print_scene(const std::string& input_path) noexcept;

gfx::hdr_image render_and_print(const gfx::scene& s) noexcept;

bool write_and_print_png(const gfx::hdr_image& img,
                         const std::string& output_path) noexcept;

int main(int argc, char** argv) {

  if (argc != 3) {
    print_usage();
    return -1;
  }
  const std::string input_path = argv[1],
                    output_path = argv[2];

  auto maybe_scene = read_and_print_scene(input_path);
  if (!maybe_scene) {
    return 1;
  }

  auto image = render_and_print(*maybe_scene);

  if (!write_and_print_png(image, output_path)) {
    return 1;
  }

  return 0;
}

void print_usage() noexcept {
  std::cout << "mrraytracer usage:" << std::endl << std::endl
            << "./mrraytracer <INPUT-JSON-SCENE-PATH> <OUTPUT-PNG-PATH>"
            << std::endl << std::endl;
}

std::optional<gfx::scene>
read_and_print_scene(const std::string& input_path) noexcept {

  auto p = [](const gfx::hdr_rgb& c) {
    return (std::string("<") +
            std::to_string(c.r()) + ", " +
            std::to_string(c.g()) + ", " +
            std::to_string(c.b()) +
            ">");
  };

  try {

    gfx::scene s(gfx::scene::read_json(input_path));
    assert(s.complete());

    const std::string indent("    "),
                      indent2 = indent + indent;

    std::cout << "scene has:" << std::endl
              << indent << "camera"
                        << " eye=" << s.camera().eye()
                        << " u=" << s.camera().u()
                        << " v=" << s.camera().v()
                        << " w=" << s.camera().w()
                        << std::endl
              << indent << "viewport"
                        << " x_res=" << s.viewport().x_resolution()
                        << " y_res=" << s.viewport().y_resolution()
                        << " left=" << s.viewport().left()
                        << " top=" << s.viewport().top()
                        << " right=" << s.viewport().right()
                        << " bottom=" << s.viewport().bottom()
                        << std::endl;

    if (auto ptr = dynamic_cast<const gfx::orthographic_projection*>(&s.projection())) {
      std::cout << indent << "projection is ORTHOGRAPHIC" << std::endl;
    } else if (auto ptr = dynamic_cast<const gfx::perspective_projection*>(&s.projection())) {
      std::cout << indent << "projection is PERSPECTIVE with"
                << " focal_length=" << ptr->focal_length()
                << std::endl;
    } else {
      assert(false); // should never get here, projection is unknown type
    }

    if (auto ptr = dynamic_cast<const gfx::flat_shader*>(&s.shader())) {
      std::cout << indent << "shader is FLAT" << std::endl;
    } else if (auto ptr = dynamic_cast<const gfx::blinn_phong_shader*>(&s.shader())) {
      std::cout << indent << "shader is BLINN-PHONG with" << std::endl
                << indent2 << " ambient_coefficient=" << ptr->ambient_coefficient() << std::endl
                << indent2 << " ambient_color=" << p(ptr->ambient_color()) << std::endl
                << indent2 << " diffuse_coefficient=" << ptr->diffuse_coefficient() << std::endl
                << indent2 << " specular_coefficient=" << ptr->specular_coefficient() << std::endl;
    } else {
      assert(false); // should never get here, shader is unknown type
    }

    std::cout << indent << "background=" << p(s.background()) << std::endl;

    std::cout << indent << "count of lights=" << s.lights().size() << std::endl;
    for (auto& light : s.lights()) {
      std::cout << indent2
                << "location=" << light->location()
                << " color=" << p(light->color())
                << " intensity=" << light->intensity()
                << std::endl;
    }

    std::cout << indent << "count of objects=" << s.objects().size() << std::endl;

    return std::make_optional<gfx::scene>(std::move(s));

  } catch (gfx::scene_read_exception e) {
    std::cerr << "error loading scene \"" + input_path + "\": "
              << e.message() << std::endl;
    return std::nullopt;
  }
}

gfx::hdr_image render_and_print(const gfx::scene& s) noexcept {

  std::cout << "raytracing...";
  auto image = s.render();
  std::cout << "done" << std::endl;

  return image;
}

bool write_and_print_png(const gfx::hdr_image& img,
                         const std::string& output_path) noexcept {

  std::cout << "writing " << output_path << "...";
  if (!gfx::write_png(img, output_path)) {
    std::cerr << std::endl << "error: could not write PNG" << std::endl;
    return false;
  }
  std::cout << "done" << std::endl << std::endl;
  return true;
}
