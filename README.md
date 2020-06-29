# Raytracer
Minimalist Raytracer in C++ that supports the fundamental raytracing algorithm, orthographic and perspective projection, flat and Blinn-Phong shading, point lights lighting, and scene objects comprised of sphere and triangle primitives. 

### Note: 
Due to the group project nature of this project, **this README file highlights only my own contributions to the project**. Skeleton code was provided by Dr. Kevin Wortman during a Principles of Computer Graphics class and the project was implemented using a Linux environment. 

## About
This repository contains the implementation files of a minimalist raytracer built using C++. 

**This raytracer supports:** fundamental raytracing algorithm, orthographic and perspective projections, flat and Blinn-Phong shaders, point lights only lighting, and the rendering of basic scene objects comprised of spheres and triangles (including a tea pot scene). 

**Some future work for this raytracer includes:** adding shadowing features, mirrored surfaces, other scene object shapes (besides triangles and spheres), and additional lighting features. 

## Linear Algebra Contributions
This project builds on a previous linear algebra project I built in which I implemented vector and matrix algebra operations in C++. The following linear algebra contributions were all performed by me using skeleton code provided by Dr. Kevin Wortman. Code for all math implementations can be found in the gfxalgebra.hpp file, which contains my implementation of basic addition, subtraction, multiplication by scalar, and division by scalar operations for both vector and matrix classes as well as:

*Additional vector functions:* Dot and cross products, vector equality check functions, resizing operations (grow/shrink/subvector conversion), unit vector check and normalization, vector conversion to matrix row/column functions.

*Additional matrix functions:* matrix by matrix multiplication, matrix by scalar multiplication, matrix equality checks, matrix conversion functions (column n as a width-1 matrix, column n to vector, row n as a height-1 matrix, row n to vector, and submatrix).

*Additional linear algebra functions:* Finding the determinant of a matrix, solving a linear system, finding identity matrix, and matrix transpose.

## The Raytracer
All of my contributions for this raytracer can be found in the gfxraytrace.hpp file in this repository, and can be broken down as follows:

### The Basic Raytracer Algorithm:
The structure of this basic ray tracing program follows the basic raytracing approach:

```
for each pixel in the scene we must:
  compute the viewing ray
  find: first object hit by the ray and the surface normal
  set the pixel color to the value computed using the hit point, light and surface normal
```

For this project, I implemented this idea as follows: (see: scene::render() in gfxraytrace.hpp)

```
// For each pixel:
for (size_t y = 0; y < h; ++y) {
    for (size_t x = 0; x < w; ++x) {

      // Compute the viewing ray
      vector2<double> uv = viewport_->uv(x,y);
      view_ray viewing_ray = projection_->compute_view_ray(camera(), uv[0], uv[1]);
      
      // Find first object hit by the ray and the surface normal
      std::optional<intersection> intersection_check = intersect(viewing_ray);
      
      // Color the pixel
      if (intersection_check == std::nullopt) {
        result.pixel(x, y, background_);
      }
      else {
        hdr_rgb pixel_color = shader_->shade(*this, camera(), *intersection_check);
        result.pixel(x, y, pixel_color);
      }

    }
  }
```

### Ray-Object Intersection:
In order for the raytracer to find the ray-object intersection required to calculate what color each pixel must be shaded as, the program follows the approach:

#### For Ray-Spehere Intersections:
Intersection points occur when points on the raysatisfy the implicit equation which can be solved using a quadratic equation approach.

[]

where:

    *t*, the intersection

      d, the view direction
      
      e, the origin
      
      c, the center of the sphere
      
      R, the radius of the sphere

This solution can be simplified by calculating the discriminant (B^2 - 4AC) where: 

    *A* = d^2 

    *B* = d * (e-c)

    *C* = (e-c)^2 - R^2


Note that by calculating the discriminant we can find the following:
- if the discriminant is negative: the line and the sphere do not intersect
- if the discriminant is positive: there are two solutions (one where the ray enters the sphere and one where it leaves)
- if the discriminant is zero: the ray grazes the sphere (touches it at exactly one point)

For the c++ implementation see: scene_sphere::intersect()

#### For Ray-Triangle Intersections:
We use a barycentric coordinates approach and matrix algebra to solve a linear system that solves the equation:

[]

where:

      d, the view direction

      e, the origin
      
      a,b,c, the triangle vertices
      
      β, γ, *t*, unknowns we solve for to calculate intersection. 
      
The intersection is inside the triangle if and only if β>0, γ>0, and β+γ<1.

For the c++ implementation see: scene_triangle::intersect()

### Perspective:
This raytracer supports the rendering of scenes using orthographic and perspective projection.

#### Ortographic Projection
Ortographic view in this raytracer follows the approach:

**Set the ray direction following:**

the view direction -w

**Set the ray origin following:**

*uv* coordinated of viewport (pixel's positions on the image plane, measured respect to origin e and basis {u,v}
Find the ray origin using the equation e + *u*u + *v*v


I implemented this idea as follows: (see: orthographic_projection::compute_view_ray() in gfxraytrace.hpp)

```
  // Find ray origin
  vector3<double> origin = c.eye() + (c.u() * u) + (c.v() * v);
  // Set ray origin and direction (using the opposite of the view direction as the direction)
  return view_ray(origin, -c.w());
```

#### Perspective Projection
Perspective view in this raytracer follows the approach:

**Set the ray direction following:**

*uv* coordinated of viewport (pixel's positions on the image plane), w the image plane normal, image plane positioned at some distance *d* in front of the origin (image plane distance or focal length), so the direction is defined by the viewpoint and the position of the pixel in the image plane.

Find the ray direction using the equation -*d*w + *u*u + *v*v

**Set the ray origin following:**

For perspective view all rays have origin at the viewpoint e. 


I implemented this idea as follows: (see: perspective_projection::compute_view_ray() in gfxraytrace.hpp)

```
  // Find ray direction
  vector3<double> direction = -c.w() + (c.u() * u) + (c.v() * v);
  // Set the ray origin and direction
  return view_ray(c.eye(), direction);
```

### Shading
This raytracer supports the rendering of scenes using flat and Blinn-Phong shading.

#### Flat Shading
Flat shading is implemented as the simplest form of shading, simply returning the color the pixel should be shaded as at the intersection. 

See flat_shader::shade() for c++ implementation.

#### Blinn-Phong Shading
Follows the simple model for specular highlights that uses the Phong reflection model as modified by Jim Blinn:

[]
where:

        *L*, Pixel color

       *kd*, diffuse coefficient (surface color)
       
       *I*, intensity of the light source
       
       n and l, unit vectors
       
       h, half vector
       
       p, Phong exponent that controls the apparent shininess of the surface
       
       *ks*, specular coefficient (specular color of the surface)
       
See blinn_phong_shader::shade() for c++ implementation.

### Additional Contributions
Some of my additional contributions include:
- implementation of function that calculates u, v viewport coordinates (See  viewport::uv())
- implementation of camera constructor that computes the basis in terms of a given view-direction and up vector (See camera::camera())
- implementation of helper function that traces a ray to find the closest intersecting scene object (See scene::intersect())

## Final Results
Scenes with ortographic projection and flat shading:

Scenes with ortographic projection and Blinn-Phong shading:

Scenes with perspective projection and flat shading:

Scenes with perspective projection and Blinn-Phong shading:

Teapot scene:

## Usage
- Run `make` to raytrace the fast scenes (all except *teatime.png*)
- Run `make all` to raytrace all the scenes (including *teatime.png*)
