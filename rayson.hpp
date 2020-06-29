
///////////////////////////////////////////////////////////////////////////////
// rayson.hpp
///////////////////////////////////////////////////////////////////////////////

#include <array>
#include <cmath>
#include <fstream>
#include <optional>
#include <memory>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include <iostream>

#include <nlohmann/json.hpp>

namespace rayson {

  class vector3 {
  private:
    double x_, y_, z_;

  public:

    constexpr vector3() noexcept
    : x_(0.0), y_(0.0), z_(0.0) { }

    constexpr vector3(double x, double y, double z) noexcept
    : x_(x), y_(y), z_(z) { }

    constexpr double x() const noexcept { return x_; }
    constexpr double y() const noexcept { return y_; }
    constexpr double z() const noexcept { return z_; }

    bool is_normalized() const noexcept {
      auto magnitude = std::sqrt(x_*x_ + y_*y_ + z_*z_),
           delta = magnitude - 1.0,
           distance = std::abs(delta);
      return distance <= .01;
    }

    constexpr bool operator==(const vector3& rhs) const noexcept {
      return (x_ == rhs.x_) && (y_ == rhs.y_) && (z_ == rhs.z_);
    }
    constexpr bool operator!=(const vector3& rhs) const noexcept {
      return !(*this == rhs);
    }
  };

  class color {
  private:
    double r_, g_, b_;

  public:

    constexpr color() noexcept
    : r_(0.0), g_(0.0), b_(0.0) { }

    constexpr color(double r, double g, double b) noexcept
    : r_(r), g_(g), b_(b) {
      assert((r >= 0.0) && (r <= 1.0));
      assert((g >= 0.0) && (g <= 1.0));
      assert((b >= 0.0) && (b <= 1.0));
    }

    constexpr double r() const noexcept { return r_; }
    constexpr double g() const noexcept { return g_; }
    constexpr double b() const noexcept { return b_; }
  };

  class camera {
  private:
    vector3 eye_, up_, view_;

  public:

    constexpr camera(
      const vector3& eye,
      const vector3& up,
      const vector3& view
    ) noexcept
    : eye_(eye), up_(up), view_(view) { }

    constexpr const vector3& eye () const noexcept { return eye_   ; }
    constexpr const vector3& up  () const noexcept { return up_    ; }
    constexpr const vector3& view() const noexcept { return view_  ; }
  };

  class viewport {
  private:
    unsigned x_resolution_, y_resolution_;
    double left_, top_, right_, bottom_;

  public:

    // left: negative
    // top: positive
    // right: positive
    // bottom: negative
    constexpr viewport(
      unsigned x_resolution,
      unsigned y_resolution,
      double left,
      double top,
      double right,
      double bottom) noexcept
    : x_resolution_(x_resolution),
      y_resolution_(y_resolution),
      left_(left),
      top_(top),
      right_(right),
      bottom_(bottom) {
      assert(x_resolution > 0);
      assert(y_resolution > 0);
      assert(left < 0.0);
      assert(right > 0.0);
      assert(top > 0.0);
      assert(bottom < 0.0);
    }

    constexpr unsigned x_resolution() const noexcept { return x_resolution_; }
    constexpr unsigned y_resolution() const noexcept { return y_resolution_; }

    constexpr double left  () const noexcept { return left_  ; }
    constexpr double top   () const noexcept { return top_   ; }
    constexpr double right () const noexcept { return right_ ; }
    constexpr double bottom() const noexcept { return bottom_; }
  };

  class ortho_projection {

  };

  class persp_projection {
  private:
    double focal_length_;

  public:

    constexpr persp_projection(double focal_length) noexcept
    : focal_length_(focal_length) {
      assert(focal_length > 0.0);
    }

    constexpr double focal_length() const noexcept { return focal_length_; }
  };

  using projection = std::variant<ortho_projection, persp_projection>;

  class flat_shader {

  };

  class phong_shader {
  private:
    double ambient_coeff_, diffuse_coeff_, specular_coeff_;
    color ambient_color_;

  public:
    phong_shader(
      double ambient_coeff,
      double diffuse_coeff,
      double specular_coeff,
      const color& ambient_color) noexcept
    : ambient_coeff_(ambient_coeff),
      diffuse_coeff_(diffuse_coeff),
      specular_coeff_(specular_coeff),
      ambient_color_(ambient_color) {

      assert(ambient_coeff >= 0.0);
      assert(diffuse_coeff >= 0.0);
      assert(specular_coeff >= 0.0);
    }

    constexpr double ambient_coeff () const noexcept { return ambient_coeff_ ; }
    constexpr double diffuse_coeff () const noexcept { return diffuse_coeff_ ; }
    constexpr double specular_coeff() const noexcept { return specular_coeff_; }

    constexpr const color& ambient_color() const noexcept { return ambient_color_; }
  };

  using shader = std::variant<flat_shader, phong_shader>;

  class material {
  private:
    std::string name_;
    double shininess_;
    color color_;

  public:

    material(
      const std::string& name,
      double shininess,
      const color& color
      ) noexcept
    : name_(name), shininess_(shininess), color_(color) {
      assert(shininess > 0.0);
    }

    material(
      const std::string&& name,
      double shininess,
      const color& color
      ) noexcept
    : name_(name), shininess_(shininess), color_(color) {
      assert(shininess > 0.0);
    }

    constexpr const std::string& name() const noexcept { return name_; }
    constexpr double shininess() const noexcept { return shininess_; }
    constexpr const color& color() const noexcept { return color_; }
  };

  class point_light {
  private:
    vector3 location_;
    color color_;
    double intensity_;

  public:

    constexpr point_light(
      const vector3& location,
      const color& color,
      double intensity
      ) noexcept
    : location_(location), color_(color), intensity_(intensity) {
      assert(intensity > 0.0);
    }

    constexpr const vector3& location() const noexcept { return location_; }
    constexpr const color& color() const noexcept { return color_; }
    constexpr double intensity() const noexcept { return intensity_; }
  };

  class sphere {
  private:
    const material* material_;
    vector3 center_;
    double radius_;

  public:

    sphere(
      const material* material,
      const vector3& center,
      double radius
      ) noexcept
    : material_(material), center_(center), radius_(radius) {
      assert(material != nullptr);
      assert(radius > 0.0);
    }

    constexpr const material& material() const noexcept { return *material_; }
    constexpr const vector3& center() const noexcept { return center_; }
    constexpr double radius() const noexcept { return radius_; }
  };

  class triangle {
  private:
    const material* material_;
    vector3 a_, b_, c_;

  public:
    triangle(
      const material* material,
      const vector3& a,
      const vector3& b,
      const vector3& c)
    : material_(material), a_(a), b_(b), c_(c) {
      assert(material != nullptr);
      assert(a != b);
      assert(a != c);
      assert(b != c);
    }

    constexpr const material& material() const noexcept { return *material_; }
    constexpr const vector3& a() const noexcept { return a_; }
    constexpr const vector3& b() const noexcept { return b_; }
    constexpr const vector3& c() const noexcept { return c_; }
  };

  class scene {
  public:

    using material_container = std::vector<material>;
    using point_light_container = std::vector<point_light>;
    using sphere_container = std::vector<sphere>;
    using triangle_container = std::vector<triangle>;

  private:
    camera camera_;
    viewport viewport_;
    projection projection_;
    shader shader_;
    color background_;
    material_container materials_;
    point_light_container point_lights_;
    sphere_container spheres_;
    triangle_container triangles_;

  public:

    scene(
      camera&& camera,
      viewport&& viewport,
      projection&& projection,
      shader&& shader,
      const color& background
    ) noexcept
    : camera_(camera),
      viewport_(viewport),
      projection_(projection),
      shader_(shader),
      background_(background) { }

    constexpr const camera&     camera    () const noexcept { return camera_;     }
    constexpr const viewport&   viewport  () const noexcept { return viewport_;   }
    constexpr const projection& projection() const noexcept { return projection_; }
    constexpr const shader&     shader    () const noexcept { return shader_;     }
    constexpr const color&      background() const noexcept { return background_; }

    constexpr const point_light_container& point_lights() const noexcept { return point_lights_; }
    constexpr const material_container&    materials   () const noexcept { return materials_;    }
    constexpr const sphere_container&      spheres     () const noexcept { return spheres_;      }
    constexpr const triangle_container&    triangles   () const noexcept { return triangles_;    }

    void emplace_point_light (point_light&& x) noexcept { point_lights_.emplace_back(x); }
    void emplace_material    (material&&    x) noexcept { materials_   .emplace_back(x); }
    void emplace_sphere      (sphere&&      x) noexcept { spheres_     .emplace_back(x); }
    void emplace_triangle    (triangle&&    x) noexcept { triangles_   .emplace_back(x); }
  };

  // An error encountered while trying to read and parse a scene file.
  class read_exception {
  private:
    std::string message_;

  public:

    read_exception(const std::string& message)
    : message_(message) { }

    read_exception(std::string&& message)
    : message_(message) { }

    constexpr const std::string& message() const noexcept { return message_; }
  };

  scene read_json(const nlohmann::json& j) {

    if (!j.is_object()) {
      throw read_exception("rayson must be comprised of one JSON object");
    }

    auto has = [](auto& j_obj, const std::string& key) noexcept {
      return j_obj.find(key) != j_obj.end();
    };

    auto check_has_key = [&](auto& j_obj, const std::string& key) {
      if (!has(j_obj, key)) {
        throw read_exception("missing key \"" + key + "\"");
      }
    };

    auto get_bool = [&](auto& j_obj, const std::string& key) {
      check_has_key(j_obj, key);
      if (!j_obj[key].is_boolean()) {
        throw read_exception("key \"" + key + "\" must be a bool");
      }
      return j_obj[key].template get<bool>();
    };

    auto get_double = [&](auto& j_obj, const std::string& key) {
      check_has_key(j_obj, key);
      if (!j_obj[key].is_number_float()) {
        throw read_exception("key \"" + key + "\" must be a float");
      }
      return j_obj[key].template get<double>();
    };

    auto get_string = [&](auto& j_obj, const std::string& key) {
      check_has_key(j_obj, key);
      if (!j_obj[key].is_string()) {
        throw read_exception("key \"" + key + "\" must be a string");
      }
      return j_obj[key].template get<std::string>();
    };

    auto get_positive_unsigned = [&](auto& j_obj, const std::string& key) {
      check_has_key(j_obj, key);
      if (!j_obj[key].is_number_integer()) {
        throw read_exception("key \"" + key + "\" must be an integer");
      }
      int x = j_obj[key].template get<int>();
      if (x <= 0) {
        throw read_exception("key \"" + key + "\" must be positive");
      }
      return static_cast<unsigned>(x);
    };

    auto get_vector3 = [&](auto& j_obj, const std::string& key) {
      check_has_key(j_obj, key);
      auto& a = j_obj.at(key);
      if (!a.is_array()) {
        throw read_exception("expected " + key + " to be an array");
      }
      if (a.size() != 3) {
        throw read_exception("expected array " + key + " to have 3 elements, but found " +
                             std::to_string(j_obj.size()) + " elements");
      }
      try {
        return vector3(a.at(0).template get<double>(),
                       a.at(1).template get<double>(),
                       a.at(2).template get<double>());
      } catch (std::exception& e) {
        throw read_exception("vector3 must contain numbers");
      }
    };

    auto get_positive_double = [&](auto& j_obj, const std::string& key) {
      double x = get_double(j_obj, key);
      if (x <= 0.0) {
        throw read_exception(std::string("key \"") + key + "\" must be positive");
      }
      return x;
    };

    auto get_negative_double = [&](auto& j_obj, const std::string& key) {
      double x = get_double(j_obj, key);
      if (x >= 0.0) {
        throw read_exception(std::string("key \"") + key + "\" must be negative");
      }
      return x;
    };

    auto get_nonnegative_double = [&](auto& j_obj, const std::string& key) {
      double x = get_double(j_obj, key);
      if (x < 0.0) {
        throw read_exception(std::string("key \"") + key + "\" must be non-negative");
      }
      return x;
    };

    auto get_color = [&](auto& j_obj, const std::string& key) {
      auto vect = get_vector3(j_obj, key);
      if ((vect.x() < 0.0) || (vect.x() > 1.0)) {
        throw read_exception("color has r component outside the range [0, 1]");
      }
      if ((vect.y() < 0.0) || (vect.y() > 1.0)) {
        throw read_exception("color has g component outside the range [0, 1]");
      }
      if ((vect.z() < 0.0) || (vect.z() > 1.0)) {
        throw read_exception("color has b component outside the range [0, 1]");
      }
      return color(vect.x(), vect.y(), vect.z());
    };

    camera camera(get_vector3(j, "camera_eye"),
                  get_vector3(j, "camera_up"),
                  get_vector3(j, "camera_view"));

    viewport viewport(get_positive_unsigned(j, "x_resolution"),
                      get_positive_unsigned(j, "y_resolution"),
                      get_negative_double(j, "viewport_left"),
                      get_positive_double(j, "viewport_top"),
                      get_positive_double(j, "viewport_right"),
                      get_negative_double(j, "viewport_bottom"));

    projection projection;
    {
      const std::string ortho_key = "ortho_projection",
                        persp_key = "persp_focal_length";
      bool has_orth = has(j, ortho_key), has_persp = has(j, persp_key);
      if (has_orth && has_persp) {
        throw read_exception("cannot have both " + ortho_key + " and " + persp_key);
      }  else if (!has_orth && !has_persp) {
        throw read_exception("must have " + ortho_key + " or " + persp_key);
      } else if (has_orth) {
        bool b = get_bool(j, ortho_key);
        if (b == false) {
          throw read_exception(ortho_key + ", if present, must be true");
        }
        projection = ortho_projection();
      } else {
        assert(has_persp);
        projection = persp_projection(get_positive_double(j, persp_key));
      }
    }

    shader shader;
    if (has(j, "flat_shader")) {
      bool b = get_bool(j, "flat_shader");
      if (b == false) {
        throw read_exception("flat_shader, if present, must be true");
      }
      shader = flat_shader();
    } else if (has(j, "phong_shader")) {
      auto child = j["phong_shader"];
      shader = phong_shader(get_nonnegative_double(child, "ambient_coeff"),
                            get_nonnegative_double(child, "diffuse_coeff"),
                            get_nonnegative_double(child, "specular_coeff"),
                            get_color(child, "ambient_color"));
    } else {
      throw read_exception("must have flat_shader or phong_shader");
    }

    color background = get_color(j, "background");

    auto result = scene(std::move(camera),
                        std::move(viewport),
                        std::move(projection),
                        std::move(shader),
                        background);

    if (has(j, "point_lights")) {
      auto& child = j["point_lights"];
      if (!child.is_array()) {
        throw read_exception("expected point_lights to be an array");
      }
      for (auto& it : child) {
        result.emplace_point_light(point_light(get_vector3(it, "location"),
                                               get_color(it, "color"),
                                               get_positive_double(it, "intensity")));
      }
    }

    if (!has(j, "materials")) {
      throw read_exception("no materials list");
    } else {
      auto& child = j["materials"];
      for (auto& it : child) {
        result.emplace_material(material(get_string(it, "name"),
                                         get_positive_double(it, "shininess"),
                                         get_color(it, "color")));
      }
    }

    std::unordered_map<std::string, const material*> material_map;
    for (auto& m : result.materials()) {
      auto& key = m.name();
      if (material_map.count(key) > 0) {
        throw read_exception("duplicate material name \"" + key + "\"");
      } else {
        material_map[key] = &m;
      }
    }

    if (has(j, "spheres")) {
      for (auto& i : j["spheres"]) {
        auto material_name = get_string(i, "material");
        auto material_found = material_map.find(material_name);
        if (material_found == material_map.end()) {
          throw read_exception("sphere references undefined material \"" + material_name + "\"");
        }
        result.emplace_sphere(sphere(material_found->second,
                                     get_vector3(i, "center"),
                                     get_positive_double(i, "radius")));
      }
    }

    if (has(j, "triangles")) {
      for (auto& i : j["triangles"]) {

        auto material_name = get_string(i, "material");
        auto a = get_vector3(i, "a");
        auto b = get_vector3(i, "b");
        auto c = get_vector3(i, "c");

        if ((a == b) || (a == c) || (b == c)) {
          throw read_exception("triangle is degenerate due to duplicated vertices");
        }

        auto material_found = material_map.find(material_name);
        if (material_found == material_map.end()) {
          throw read_exception("triangle references undefined material \"" + material_name + "\"");
        }

        result.emplace_triangle(triangle(material_found->second, a, b, c));
      }
    }

    return result;
  }

  scene read_file(const std::string& path) {

    std::ifstream f(path);
    if (!f) {
      throw read_exception("could not open \"" + path + "\"");
    }

    nlohmann::json j;
    try {
      f >> j;
    } catch (std::invalid_argument e) {
      f.close();
      throw read_exception("JSON parse error reading \"" + path + "\"");
    }

    f.close();

    return read_json(j);
  }
}
