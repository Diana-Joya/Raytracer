
///////////////////////////////////////////////////////////////////////////////
// gfxalgebra.hpp
//
// Linear algebra for graphics.
//
// This header defines a gfx::vector class representing a linear algebra
// vector, and a gfx::matrix class representing a linear algebra matrix.
//
// This file builds upon gfxnumeric.hpp, so you may want to familiarize
// yourself with that header before diving into this one.
//
// Students: all of your work should go in this file, and the only files that
// you need to modify in project 2 are this file, and README.md.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include <functional>
#include <numeric>

#include "gfxnumeric.hpp"

namespace gfx {

// Forward-declare the matrix type, because a few vector operations involve
// matrices.
template <typename scalar_type,
          size_t HEIGHT,
          size_t WIDTH>
class matrix;

// A mathematical vector, in the spirit of linear algebra.
//
// This is very different from a general-purpose self-resizing array data
// structure such as std::vector. gfx::vector has a fixed dimension (size)
// and only supports mathematical operations.
//
// scalar_type is the type of each element, which must be a numeric type
// such as double, float, or int. At a minimum, it must be possible to
// assign a scalar_type to 0.
//
// DIMENSION is the size of the vector. Dimension should be positive;
// zero-dimension vectors are technically supported, but seem pointless.
template <typename scalar_type,
          size_t DIMENSION>
class vector {
public:

  // Type aliases.
  using same_type = gfx::vector<scalar_type, DIMENSION>;
  using storage_type = std::array<scalar_type, DIMENSION>;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;

private:

  storage_type elements_;

public:

  ////////
  // Constructors and destructor.
  ////////

  // Default constructor. Every element is initialized to zero.
  constexpr vector() noexcept
  : vector(0) { }

  // Copy and move contructors.
  constexpr vector(const same_type&) noexcept = default;
  constexpr vector(same_type&&) noexcept = default;

  // Fill constructor. Every element is initialized to default_value.
  constexpr vector(scalar_type default_value) noexcept { fill(default_value); }

  // Iterator constructor. If the iterator range has fewer than DIMENSION
  // elements, the unspecified elements default to 0. If the iterator range
  // has extra elements, the extras are ignored.
  template <typename input_iterator>
  constexpr vector(input_iterator first, input_iterator last) noexcept {
    auto iter = first;
    for (size_t i = 0; i < DIMENSION; ++i) {
      elements_[i] = (iter == last) ? scalar_type(0) : *iter++;
    }
  }

  // Initializer list constructor. If the list has fewer than DIMENSION
  // elements, the unspecified elements default to 0. If the list
  // has extra elements, the extras are ignored.
  constexpr vector(std::initializer_list<scalar_type> il) noexcept
  : vector(il.begin(), il.end()) { }

  // Destructor.
  ~vector() = default;

  ////////
  // Operator overloads.
  ////////

  constexpr same_type& operator= (const same_type&) noexcept = default;

  constexpr bool operator== (const same_type& rhs) const noexcept {
    return std::equal(elements_.begin(), elements_.end(), rhs.elements_.begin());
  }

  constexpr bool operator!= (const same_type& rhs) const noexcept {
    return !(*this == rhs);
  }

  constexpr const scalar_type& operator[](size_t i) const noexcept {
    assert(is_index(i));
    return elements_[i];
  }

  constexpr scalar_type& operator[](size_t i) noexcept {
    assert(is_index(i));
    return elements_[i];
  }

  constexpr same_type operator+(const same_type& rhs) const noexcept {
    same_type res;

    for (size_t i = 0; i < DIMENSION; i++) {
      res[i] = elements_[i] + rhs[i];
    }

    return res;
  }

  constexpr same_type operator-() const noexcept {
    same_type res;

    for (size_t i = 0; i < DIMENSION; i++) {
      res[i] = elements_[i] * -1;
    }

    return res;
  }

  constexpr same_type operator-(const same_type& rhs) const noexcept {
    same_type res;

    for (size_t i = 0; i < DIMENSION; i++) {
      res[i] = elements_[i] - rhs[i];
    }

    return res;
  }

  // Vector-scalar product.
  constexpr same_type operator*(scalar_type rhs) const noexcept {
    same_type res;

    for (size_t i = 0; i < DIMENSION; i++) {
      res[i] = elements_[i] * rhs;
    }

    return res;
  }

  // Vector-vector product (dot product).
  constexpr scalar_type operator*(const same_type& rhs) const noexcept {
    scalar_type res = 0;

    for (size_t i = 0; i < DIMENSION; i++) {
      res += elements_[i] * rhs[i];
    }

    return res;
  }

  // Vector divided by scalar.
  constexpr same_type operator/(scalar_type rhs) const noexcept {
    same_type res;

    for (size_t i = 0; i < DIMENSION; i++) {
      res[i] = elements_[i] / rhs;
    }

    return res;
  }

  // Stream insertion operator, for printing.
  friend std::ostream& operator<<(std::ostream& stream, const same_type& rhs) {
    stream << '<';
    if constexpr (DIMENSION > 0) {
      stream << rhs.elements_[0];
    }
    for (size_t i = 1; i < DIMENSION; ++i) {
      stream << ", " << rhs.elements_[i];
    }
    stream << '>';
    return stream;
  }

  ////////
  // Approximate equality.
  ////////

  // Return true iff each element of this vector is approximately equal to
  // the corresponding element of other, using delta, as determined by
  // gfx::approx_equal.
  constexpr bool approx_equal(const same_type& other, scalar_type delta) const noexcept {
    for (size_t i = 0; i < DIMENSION; i++) {
      if (gfx::approx_equal(elements_[i], other[i], delta) == false) return false;
    }

    return true;
  }

  ////////
  // Iterators.
  ////////

  constexpr const_iterator begin() const noexcept { return elements_.cbegin(); }
  constexpr const_iterator end  () const noexcept { return elements_.cend  (); }

  constexpr iterator begin() noexcept { return elements_.begin(); }
  constexpr iterator end  () noexcept { return elements_.end  (); }

  ////////
  // Size and indices.
  ////////

  constexpr size_t dimension() const noexcept { return DIMENSION; }

  constexpr bool is_index(size_t i) const noexcept { return (i < DIMENSION); }

  ////////
  // Converting to other types.
  ////////

  // Return a vector of size NEW_DIMENSION, based on this vector.
  // NEW_DIMENSION must be greater than DIMENSION.
  // The first DIMENSION elements are copied from this vector.
  // The remaining, newly-created elements are all initialized to default_value.
  template <size_t NEW_DIMENSION>
  vector<scalar_type, NEW_DIMENSION>
  grow(scalar_type default_value = 0) const noexcept {

    static_assert(NEW_DIMENSION > DIMENSION,
                  "new dimension must be larger than old dimension");

    vector<scalar_type, NEW_DIMENSION> res;
    size_t count = 0;
    while (count < DIMENSION) {
      res[count] = elements_[count];
      count++;
    }

    while (count < NEW_DIMENSION) {
      res[count] = default_value;
      count++;
    }

    return res;
  }

  // Return a vector of size NEW_DIMENSION, based on this vector.
  // NEW_DIMENSION must be less than DIMENSION.
  // The returned vector contains the first NEW_DIMENSION elements of this vector.
  template <size_t NEW_DIMENSION>
  vector<scalar_type, NEW_DIMENSION>
  shrink() const noexcept {

    static_assert(NEW_DIMENSION < DIMENSION,
                  "new dimension must be smaller than old dimension");

    vector<scalar_type, NEW_DIMENSION> res;

    for (size_t i = 0; i < NEW_DIMENSION; i++) {
      res[i] = elements_[i];
    }

    return res;
  }

  // Return a vector of size NEW_DIMENSION, based on this vector.
  // Copies NEW_DIMENSION elements, starting at index start.
  // The specified range of indices must all fit within this vector.
  template <size_t NEW_DIMENSION>
  vector<scalar_type, NEW_DIMENSION>
  subvector(size_t start = 0) const noexcept {

    static_assert(NEW_DIMENSION <= DIMENSION,
                  "new dimension cannot be larger than old dimension");
    assert((start + NEW_DIMENSION) <= DIMENSION);

    vector<scalar_type, NEW_DIMENSION> res;
    size_t count = start;

    for (size_t i = 0; i < NEW_DIMENSION; i++) {
      res[i] = elements_[count];
      count++;
    }

    return res;
  }

  // Convert this vector to a column matrix, i.e. a matrix of height
  // DIMENSION and width 1.
  // The definition of this function needs to come after gfx::matrix, and is
  // near the bottom of this source file.
  constexpr matrix<scalar_type, DIMENSION, 1> to_column_matrix() const noexcept;

  // Convert this vector to a row matrix, i.e. a matrix of height 1
  // and width DIMENSION.
  // The definition of this function needs to come after gfx::matrix, and is
  // near the bottom of this source file.
  constexpr matrix<scalar_type, 1, DIMENSION> to_row_matrix() const noexcept;

  ////////
  // Miscellaneous operations.
  ////////

  // Cross product. This function is only defined on vectors of DIMENSION 3.
  constexpr same_type cross(const same_type& rhs) const noexcept {

    static_assert(3 == DIMENSION,
                  "cross product is only defined for 3D vectors");

    scalar_type x = elements_[1] * rhs[2] - elements_[2] * rhs[1];
    scalar_type y = elements_[2] * rhs[0] - elements_[0] * rhs[2];
    scalar_type z = elements_[0] * rhs[1] - elements_[1] * rhs[0];

    same_type res{x, y, z};

    return res;
  }

  // Fill; assign every element the value x.
  constexpr void fill(scalar_type x) noexcept { elements_.fill(x); }

  // Return true iff this is a unit vector, i.e. the length of this vector is
  // approximately equal to 1, within delta.
  bool is_unit(scalar_type delta) const noexcept {
    return gfx::approx_equal(magnitude(), (scalar_type)1, delta);
  }

  // Return true iff every element of this vector is == 0.
  constexpr bool is_zero() const noexcept {
    using namespace std::placeholders;
    return std::all_of(elements_.begin(), elements_.end(),
                       std::bind(std::equal_to{}, _1, 0));
  }

  // Return the magnitude of this vector, squared (raised to the second power).
  // Computing this quantity does not require a square root, so it is faster
  // than magnitude().
  constexpr scalar_type magnitude_squared() const noexcept {
    scalar_type res = 0;
    for (size_t i = 0; i < DIMENSION; i++) {
      res += elements_[i] * elements_[i];
    }

    return res;
  }

  // Return the magnitude (length) of this vector.
  scalar_type magnitude() const noexcept {
    return sqrt(magnitude_squared());
  }

  // Return a version of this vector that has been normalized.
  // In other words, return a vector with magnitude 1.0, and the same direction
  // as this vector.
  same_type normalized() const noexcept {
    return *this / magnitude();
  }
};

// Type aliases for 2, 3, and 4-dimensional vectors.
template <typename scalar_type> using vector2 = vector<scalar_type, 2>;
template <typename scalar_type> using vector3 = vector<scalar_type, 3>;
template <typename scalar_type> using vector4 = vector<scalar_type, 4>;

// A mathematical matrix.
//
// Like gfx::vector, this is intended for linear algebra operations, and is not
// a general purpose multidemensional array data structure.
//
// scalar_type is the type of each element, which must be a number type that
// can be initialized to 0.
// HEIGHT is the number of rows in the matrix.
// WIDTH is the number of columns in the matrix.
//
// This class is intended for graphics applications where matrices are rarely
// larger than 4x4.
template <typename scalar_type,
          size_t HEIGHT,
          size_t WIDTH>
class matrix {
public:

  // Type aliases.
  using same_type = matrix<scalar_type, HEIGHT, WIDTH>;
  using row_type = vector<scalar_type, WIDTH>;
  using row_container_type = std::array<row_type, HEIGHT>;
  using row_iterator = typename row_container_type::iterator;
  using const_row_iterator = typename row_container_type::const_iterator;

private:

  // The rows (elements).
  row_container_type rows_;

public:

  ////////
  // Constructors and destructor.
  ////////

  // Default constructor, which initializes each element to 0.
  constexpr matrix() noexcept
  : matrix(scalar_type(0)) {}

  // Copy constructor.
  constexpr matrix(const same_type&) noexcept = default;

  // Move constructor.
  constexpr matrix(same_type&&) noexcept = default;

  // Fill constructor, each element is initialized to default_value.
  constexpr matrix(scalar_type default_value) noexcept { fill(default_value); }

  // Iterator constructor. Elements are filled in row-major order, i.e.
  // the first row left-to-right, then the second row left-to-right, and so on.
  // If the iterator range contains fewer elements than this matrix, the
  // unspecified elements are initialized to 0.
  // If the iterator range contains more elements than this matrix, the extras
  // are ignored.
  template <typename input_iterator>
  constexpr matrix(input_iterator first, input_iterator last) noexcept {
    input_iterator iter = first;
    for (size_t r = 0; r < HEIGHT; ++r) {
      for (size_t c = 0; c < WIDTH; ++c) {
        rows_[r][c] = (iter == last) ? scalar_type(0) : *iter++;
      }
    }
  }

  // Initializer list constructor. Elements are filled in row-major order, i.e.
  // the first row left-to-right, then the second row left-to-right, and so on.
  // If the list contains fewer elements than this matrix, the
  // unspecified elements are initialized to 0.
  // If the list contains more elements than this matrix, the extras
  // are ignored.
  constexpr matrix(std::initializer_list<scalar_type> il) noexcept
  : matrix(il.begin(), il.end()) { }

  ////////
  // Operator overloads.
  ////////

  constexpr same_type& operator= (const same_type&) noexcept = default;

  constexpr bool operator== (const same_type& rhs) const noexcept {
    return std::equal(rows_.begin(), rows_.end(), rhs.rows_.begin());
  }

  constexpr bool operator!= (const same_type& rhs) const noexcept {
    return !(*this == rhs);
  }

  constexpr const row_type& operator[] (size_t row) const noexcept {
    assert(is_row(row));
    return rows_[row];
  }

  constexpr row_type& operator[] (size_t row) noexcept {
    assert(is_row(row));
    return rows_[row];
  }

  constexpr same_type operator+ (const same_type& rhs) const noexcept {
    same_type res;
    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < WIDTH; c++) {
        res[r][c] = rows_[r][c] + rhs[r][c];
      }
    }

    return res;
  }

  // Negation.
  constexpr same_type operator- () const noexcept {
    same_type res;
    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < WIDTH; c++) {
        res[r][c] = rows_[r][c] * -1;
      }
    }

    return res;
  }

  // Matrix-matrix subtraction.
  constexpr same_type operator- (const same_type& rhs) const noexcept {
    same_type res;
    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < WIDTH; c++) {
        res[r][c] = rows_[r][c] - rhs[r][c];
      }
    }

    return res;
  }

  // Matrix-scalar multiplication.
  constexpr same_type operator* (const scalar_type rhs) const noexcept {
    same_type res;
    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < WIDTH; c++) {
        res[r][c] = rows_[r][c] * rhs;
      }
    }

    return res;
  }

  // Matrix-matrix multiplication.
  // The rhs matrix's height must be the same as this matrix' WIDTH.
  // The result vector's dimension is the same as the rhs matrix's width.
  // These dimensions are specified as template parameters, so trying to
  // multiply with invalid dimensions will cause a compile error.
  template <size_t RESULT_WIDTH>
  constexpr
  matrix<scalar_type, HEIGHT, RESULT_WIDTH>
  operator* (const matrix<scalar_type, WIDTH, RESULT_WIDTH>& rhs) const noexcept {
    matrix<scalar_type, HEIGHT, RESULT_WIDTH> res;
    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < RESULT_WIDTH; c++) {
        res[r][c] = 0;
        for (size_t k = 0; k < WIDTH; k++) {
          res[r][c] += rows_[r][k] * rhs[k][c];
        }
      }
    }
    return res;
  }

  // Division by a scalar.
  constexpr same_type operator/ (scalar_type rhs) const noexcept {
    same_type res;
    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < WIDTH; c++) {
        res[r][c] = rows_[r][c] / rhs;
      }
    }

    return res;
  }

  // Stream insertion operator, for printing.
  friend std::ostream& operator<<(std::ostream& stream, const same_type& rhs) {
    for (auto& row : rhs.rows_) {
      stream << '[';
      if constexpr (WIDTH > 0) {
        stream << row[0];
      }
      for (size_t i = 1; i < WIDTH; ++i) {
        stream << ' ' << row[i];
      }
      stream << ']' << std::endl;
    }
    return stream;
  }

  ////////
  // Approximate equality.
  ////////

  // Return true iff each element of this matrix is approximately equal to
  // the corresponding element of other, using delta, as determined by
  // gfx::approx_equal.
  constexpr bool approx_equal(const same_type& other, scalar_type delta) const noexcept {
    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < WIDTH; c++) {
        if (gfx::approx_equal(rows_[r][c], other[r][c], delta) == false) return false;
      }
    }

    return true;
  }

  ////////
  // Iterators.
  ////////

  constexpr const_row_iterator begin() const noexcept { return rows_.begin(); }
  constexpr       row_iterator begin()       noexcept { return rows_.begin(); }

  constexpr const_row_iterator end() const noexcept { return rows_.end(); }
  constexpr       row_iterator end()       noexcept { return rows_.end(); }

  ////////
  // Size and indices.
  ////////

  static constexpr size_t height() noexcept { return HEIGHT; }

  constexpr size_t width() const noexcept { return WIDTH; }

  constexpr bool is_column(size_t column) const noexcept { return (column < WIDTH); }

  constexpr bool is_row(size_t row) const noexcept { return (row < HEIGHT); }

  constexpr bool is_row_column(size_t row, size_t column) const noexcept {
    return is_row(row) && is_column(column);
  }

  template <typename OTHER_SCALAR_TYPE,
            size_t OTHER_WIDTH,
            size_t OTHER_HEIGHT>
  constexpr bool
  is_same_size(const matrix<OTHER_SCALAR_TYPE, OTHER_HEIGHT, OTHER_WIDTH>& other)
  const noexcept {
    return (width() == other.width()) && (height() == other.height());
  }

  static constexpr bool is_square() noexcept { return (WIDTH == HEIGHT); }

  // Return true iff every element of this matrix is == 0.
  constexpr bool is_zero() const noexcept {
    return std::all_of(rows_.begin(), rows_.end(),
                       [](auto& row) { return row.is_zero(); });
  }

  ////////
  // Converting to other types.
  ////////

  // Return column number c of this matrix as a width-1 matrix.
  // c must be a valid column number.
  constexpr matrix<scalar_type, HEIGHT, 1> column_matrix(size_t c) const noexcept {

    assert(is_column(c));

    matrix<scalar_type, HEIGHT, 1> res;
    for (size_t r = 0; r < HEIGHT; r++) {
        res[r][0] = rows_[r][c];
    }

    return res;
  }

  // Return column number c of this matrix as a vector.
  // c must be a valid column number.
  constexpr vector<scalar_type, HEIGHT> column_vector(size_t c) const noexcept {

    assert(is_column(c));

    vector<scalar_type, HEIGHT> res;

    for (size_t r = 0; r < HEIGHT; r++) {
        res[r] = rows_[r][c];
    }

    return res;
  }

  // Return row number r of this matrix as a height-1 matrix.
  // r must be a valid row number.
  constexpr const gfx::matrix<scalar_type, 1, WIDTH> row_matrix(size_t r) const noexcept {

    assert(is_row(r));

    matrix<scalar_type, 1, WIDTH> res;

    for (size_t c = 0; c < WIDTH; c++) {
        res[0][c] = rows_[r][c];
    }

    return res;
  }

  // Return row number r of this matrix as a vector.
  // r must be a valid row number.
  constexpr const row_type& row_vector(size_t r) const noexcept {
    assert(is_row(r));
    return rows_[r];
  }

  // Return a portion of this matrix as a new matrix object.
  // The dimensions of the sub-matrix are template parameters NEW_HEIGHT
  // and NEW_WIDTH.
  // The location of the sub-matrix is specified by top_row and left_column.
  // The specified range of locations must all fit within this matrix.
  template <size_t NEW_HEIGHT,
            size_t NEW_WIDTH>
  constexpr matrix<scalar_type, NEW_HEIGHT, NEW_WIDTH>
  submatrix(size_t top_row = 0, size_t left_column = 0) const noexcept {

    assert(is_row_column(top_row, left_column));
    assert((top_row + NEW_HEIGHT) <= HEIGHT);
    assert((left_column + NEW_WIDTH) <= WIDTH);

    matrix<scalar_type, NEW_HEIGHT, NEW_WIDTH> res;
    size_t i = top_row;
    for (size_t r = 0; r < NEW_HEIGHT; r++) {
      size_t j = left_column;
      for (size_t c = 0; c < NEW_WIDTH; c++) {
          res[r][c] = rows_[i][j];
          j++;
        }
      i++;
    }

    return res;
  }

  ////////
  // Determinant.
  ////////

  // Return the determinant of this matrix.
  // This function is only defined on square 2x2 or 3x3 matrices.
  scalar_type determinant() const noexcept {

    static_assert(is_square(),
                  "determinant is only defined for square matrices");
    static_assert((WIDTH == 2) || (WIDTH == 3),
	                "determinant only implemented for 2x2 and 3x3 matrices");

    scalar_type res = 0;
    if (WIDTH == 2) {
      res = determinant_helper(3);
    }
    else {
      scalar_type a = rows_[0][0] * determinant_helper(0);
      scalar_type b = rows_[0][1] * determinant_helper(1);
      scalar_type c = rows_[0][2] * determinant_helper(2);
      res = a - b + c;
    }

    return res;
  }

  scalar_type determinant_helper(size_t it = 0) const noexcept {
    scalar_type res = 0;
    if (it == 0 || it == 2) {
      res = (rows_[1][(it + 1)%3] * rows_[2][(it + 2)%3] - rows_[1][(it + 2)%3] * rows_[2][(it + 1)%3]);
    }
    else if (it == 1) {
      res = (rows_[1][(it + 2)%3] * rows_[2][(it + 1)%3] - rows_[1][(it + 1)%3] * rows_[2][(it + 2)%3]);
    }
    else {
      res = rows_[0][0] * rows_[1][1] - rows_[1][0] * rows_[0][1];
    }
    return res;
  }

  ////////
  // Solving linear systems.
  ////////

  // Solve a linear system.
  // This matrix is considered to be the M matrix of coefficients.
  // b is the vector containing scalars on the right-hand-side of the equations.
  // Returns the x vector of values for each variable in the system.
  // This function is only defined for square 2x2 or 3x3 matrices.
  vector<scalar_type, HEIGHT> solve(const vector<scalar_type, HEIGHT>& b) const noexcept {

    static_assert(is_square(),
                  "only square linear systems can be solved");
    static_assert((WIDTH == 2) || (WIDTH == 3),
                  "solve is only implemented for 2x2 and 3x3 matrices");

    // This solve function uses Cramer's rule to solve the linear system
    vector<scalar_type, HEIGHT> res;
    matrix<scalar_type, HEIGHT, 1> ansCol = b.to_column_matrix();
    scalar_type d = determinant();

    // Replace x/y/z column with answer column to get Dx, Dy, Dz
    for (size_t c = 0; c < WIDTH; c++) {
        matrix<scalar_type, WIDTH, HEIGHT> curr = replace_column(ansCol, c);
        res[c] = curr.determinant();
    }

    // Solve for x, y, z
    for (size_t i = 0; i < res.dimension(); i++) {
      res[i] = res[i] / d;
    }

    return res;
  }

  // Replace column with answer column
  // needed to calculate determinants for Cramer's rule
  matrix<scalar_type, WIDTH, HEIGHT> replace_column (const matrix<scalar_type, WIDTH, 1>& b, size_t pos) const noexcept {
    matrix<scalar_type, WIDTH, HEIGHT> res = *this;

    for (size_t i = 0; i < HEIGHT; i++) {
      res[i][pos] = b[i][0];
    }

    return res;
  }

  ////////
  // Miscellaneous operations.
  ////////

  // Fill; assign each element to x.
  constexpr void fill(scalar_type x) noexcept {
    std::for_each(rows_.begin(), rows_.end(),
                  [&](auto& row) { row.fill(x); });
  }

  // Create and return an identity matrix with the same dimensions as this
  // matrix.
  // This function is only defined for square matrices.
  static constexpr same_type identity() noexcept {
    same_type res;

    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < WIDTH; c++) {
        if (r == c) res[r][c] = 1;
        else res[r][c] = 0;
      }
    }

    return res;
  }

  // Return the transpose of this matrix, i.e. exchange rows with columns.
  constexpr matrix<scalar_type, WIDTH, HEIGHT> transpose() const noexcept {
    matrix<scalar_type, WIDTH, HEIGHT> res;

    for (size_t r = 0; r < HEIGHT; r++) {
      for (size_t c = 0; c < WIDTH; c++) {
        res[c][r] = rows_[r][c];
      }
    }

    return res;
  }
};

// Type aliases for 2x2, 3x3, and 4x4 square matrices.
template <typename scalar_type> using matrix2x2 = matrix<scalar_type, 2, 2>;
template <typename scalar_type> using matrix3x3 = matrix<scalar_type, 3, 3>;
template <typename scalar_type> using matrix4x4 = matrix<scalar_type, 4, 4>;

// Now we can finally define vector::to_column_matrix and vector::to_row_matrix.

template <typename scalar_type,
          size_t DIMENSION>
constexpr matrix<scalar_type, DIMENSION, 1>
vector<scalar_type, DIMENSION>::to_column_matrix() const noexcept {
  matrix<scalar_type, DIMENSION, 1> res;

  for (size_t r = 0; r < DIMENSION; r++) {
    res[r][0] = elements_[r];
  }

  return res;
}

template <typename scalar_type,
          size_t DIMENSION>
constexpr matrix<scalar_type, 1, DIMENSION>
vector<scalar_type, DIMENSION>::to_row_matrix() const noexcept {
  matrix<scalar_type, 1, DIMENSION> res;

  for (size_t c = 0; c < DIMENSION; c++) {
    res[0][c] = elements_[c];
  }

  return res;
}

} // namespace gfx
