/**
 * @file unit-capi-dense_nullable.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2020 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * Tests nullable, dense arrays.
 */

#include "catch.hpp"
#include "test/src/helpers.h"
#ifdef _WIN32
#include "tiledb/sm/filesystem/win.h"
#else
#include "tiledb/sm/filesystem/posix.h"
#endif
#include "tiledb/sm/c_api/tiledb.h"
#include "tiledb/sm/c_api/tiledb_serialization.h"
#include "tiledb/sm/misc/utils.h"
#include "tiledb/sm/serialization/query.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace tiledb::sm;
using namespace tiledb::test;

struct test_dim_t;
struct test_attr_t;
struct test_query_buffer_t;

struct DenseNullableArrayFx {
#ifdef _WIN32
  const string FILE_URI_PREFIX = "";
  const string FILE_TEMP_DIR =
      tiledb::sm::Win::current_dir() + "\\tiledb_test\\";
#else
  const string FILE_URI_PREFIX = "file://";
  const string FILE_TEMP_DIR =
      tiledb::sm::Posix::current_dir() + "/tiledb_test/";
#endif

  DenseNullableArrayFx();
  ~DenseNullableArrayFx();

  tiledb_ctx_t* ctx_;
  tiledb_vfs_t* vfs_;

  void create_dir(const string& path);
  void remove_dir(const string& path);

  void create_array(
      const string& array_name,
      const vector<test_dim_t>& test_dims,
      const vector<test_attr_t>& test_attrs,
      const tiledb_layout_t cell_order,
      const tiledb_layout_t tile_order);

  void write(
      const string& array_name,
      const vector<test_query_buffer_t>& test_query_buffers,
      const tiledb_layout_t* layout = nullptr);

  void read(
      const string& array_name,
      const vector<test_query_buffer_t>& test_query_buffers,
      const void* subarray,
      const tiledb_layout_t* layout = nullptr);
};

DenseNullableArrayFx::DenseNullableArrayFx() {
  // Create a config.
  tiledb_config_t* config = nullptr;
  tiledb_error_t* error = nullptr;
  REQUIRE(tiledb_config_alloc(&config, &error) == TILEDB_OK);
  REQUIRE(error == nullptr);

  // Create the context.
  REQUIRE(tiledb_ctx_alloc(config, &ctx_) == TILEDB_OK);
  REQUIRE(error == nullptr);

  // Create the VFS.
  REQUIRE(tiledb_vfs_alloc(ctx_, config, &vfs_) == TILEDB_OK);

  tiledb_config_free(&config);
}

DenseNullableArrayFx::~DenseNullableArrayFx() {
}

void DenseNullableArrayFx::create_dir(const std::string& path) {
  REQUIRE(tiledb_vfs_create_dir(ctx_, vfs_, path.c_str()) == TILEDB_OK);
}

void DenseNullableArrayFx::remove_dir(const string& path) {
  int is_dir = 0;
  REQUIRE(tiledb_vfs_is_dir(ctx_, vfs_, path.c_str(), &is_dir) == TILEDB_OK);
  if (is_dir)
    REQUIRE(tiledb_vfs_remove_dir(ctx_, vfs_, path.c_str()) == TILEDB_OK);
}

struct test_dim_t {
  test_dim_t(
      const string& name,
      const tiledb_datatype_t type,
      const void* const domain,
      const uint64_t tile_extent)
      : name_(name)
      , type_(type)
      , domain_(domain)
      , tile_extent_(tile_extent) {
  }

  string name_;
  tiledb_datatype_t type_;
  const void* domain_;
  uint64_t tile_extent_;
};

struct test_attr_t {
  test_attr_t(
      const string& name, const tiledb_datatype_t type, const bool nullable)
      : name_(name)
      , type_(type)
      , nullable_(nullable) {
  }

  string name_;
  tiledb_datatype_t type_;
  bool nullable_;
};

void DenseNullableArrayFx::create_array(
    const string& array_name,
    const vector<test_dim_t>& test_dims,
    const vector<test_attr_t>& test_attrs,
    const tiledb_layout_t cell_order,
    const tiledb_layout_t tile_order) {
  remove_dir(FILE_TEMP_DIR);
  create_dir(FILE_TEMP_DIR);

  // Create the dimensions.
  std::vector<tiledb_dimension_t*> dims;
  dims.reserve(test_dims.size());
  for (const auto& test_dim : test_dims) {
    tiledb_dimension_t* dim;
    const int rc = tiledb_dimension_alloc(
        ctx_,
        test_dim.name_.c_str(),
        test_dim.type_,
        test_dim.domain_,
        &test_dim.tile_extent_,
        &dim);
    REQUIRE(rc == TILEDB_OK);

    dims.emplace_back(dim);
  }

  // Create the domain.
  tiledb_domain_t* domain;
  int rc = tiledb_domain_alloc(ctx_, &domain);
  REQUIRE(rc == TILEDB_OK);
  for (const auto& dim : dims) {
    rc = tiledb_domain_add_dimension(ctx_, domain, dim);
    REQUIRE(rc == TILEDB_OK);
  }

  // Create attributes
  std::vector<tiledb_attribute_t*> attrs;
  attrs.reserve(test_attrs.size());
  for (const auto& test_attr : test_attrs) {
    tiledb_attribute_t* attr;
    rc = tiledb_attribute_alloc(
        ctx_, test_attr.name_.c_str(), test_attr.type_, &attr);
    REQUIRE(rc == TILEDB_OK);

    if (test_attr.nullable_) {
      rc = tiledb_attribute_set_nullable(ctx_, attr, 1);
      REQUIRE(rc == TILEDB_OK);
    }

    attrs.emplace_back(attr);
  }

  // Create array schema
  tiledb_array_schema_t* array_schema;
  rc = tiledb_array_schema_alloc(ctx_, TILEDB_DENSE, &array_schema);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_cell_order(ctx_, array_schema, cell_order);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_tile_order(ctx_, array_schema, tile_order);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_domain(ctx_, array_schema, domain);
  REQUIRE(rc == TILEDB_OK);
  for (const auto& attr : attrs) {
    rc = tiledb_array_schema_add_attribute(ctx_, array_schema, attr);
    REQUIRE(rc == TILEDB_OK);
  }

  // Check array schema
  rc = tiledb_array_schema_check(ctx_, array_schema);
  REQUIRE(rc == TILEDB_OK);

  // Create array
  rc = tiledb_array_create(
      ctx_, (FILE_TEMP_DIR + array_name).c_str(), array_schema);
  REQUIRE(rc == TILEDB_OK);

  // Free attributes.
  for (auto& attr : attrs) {
    tiledb_attribute_free(&attr);
  }

  // Free dimensions.
  for (auto& dim : dims) {
    tiledb_dimension_free(&dim);
  }

  // Free the domain.
  tiledb_domain_free(&domain);

  // Free the array schema.
  tiledb_array_schema_free(&array_schema);
}

struct test_query_buffer_t {
  test_query_buffer_t(
      const string& name, void* const buffer, uint64_t* const buffer_size)
      : name_(name)
      , buffer_(buffer)
      , buffer_size_(buffer_size)
      , buffer_validity_(nullptr)
      , buffer_validity_size_(nullptr) {
  }

  test_query_buffer_t(
      const string& name,
      void* const buffer,
      uint64_t* const buffer_size,
      uint8_t* const buffer_validity,
      uint64_t* const buffer_validity_size)
      : name_(name)
      , buffer_(buffer)
      , buffer_size_(buffer_size)
      , buffer_validity_(buffer_validity)
      , buffer_validity_size_(buffer_validity_size) {
  }

  string name_;
  void* buffer_;
  uint64_t* buffer_size_;
  uint8_t* buffer_validity_;
  uint64_t* buffer_validity_size_;
};

void DenseNullableArrayFx::write(
    const string& array_name,
    const vector<test_query_buffer_t>& test_query_buffers,
    const tiledb_layout_t* const layout) {
  // Open the array for writing.
  tiledb_array_t* array;
  int rc =
      tiledb_array_alloc(ctx_, (FILE_TEMP_DIR + array_name).c_str(), &array);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_WRITE);
  REQUIRE(rc == TILEDB_OK);

  // Create the write query.
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_WRITE, &query);
  REQUIRE(rc == TILEDB_OK);

  // Set the query layout, if requested.
  if (layout != nullptr) {
    rc = tiledb_query_set_layout(ctx_, query, *layout);
    REQUIRE(rc == TILEDB_OK);
  }

  // Create stack storage for validity vector objects to
  // outlive the query.
  vector<tiledb_validity_vec_t*> validity_vecs;

  // Set the query buffers.
  for (const auto& test_query_buffer : test_query_buffers) {
    if (test_query_buffer.buffer_validity_size_ == nullptr) {
      rc = tiledb_query_set_buffer(
          ctx_,
          query,
          test_query_buffer.name_.c_str(),
          test_query_buffer.buffer_,
          test_query_buffer.buffer_size_);
      REQUIRE(rc == TILEDB_OK);
    } else {
      tiledb_validity_vec_t* validity_vec = nullptr;
      rc = tiledb_validity_vec_alloc(
          ctx_, *test_query_buffer.buffer_validity_size_, &validity_vec);
      REQUIRE(rc == TILEDB_OK);
      rc = tiledb_validity_vec_set(
          ctx_, validity_vec, test_query_buffer.buffer_validity_);
      REQUIRE(rc == TILEDB_OK);

      rc = tiledb_query_set_buffer_nullable(
          ctx_,
          query,
          test_query_buffer.name_.c_str(),
          test_query_buffer.buffer_,
          test_query_buffer.buffer_size_,
          validity_vec);
      REQUIRE(rc == TILEDB_OK);

      validity_vecs.emplace_back(validity_vec);
    }
  }

  // Submit the query.
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Finalize the query, a no-op for non-global writes.
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Clean up
  rc = tiledb_array_close(ctx_, array);
  REQUIRE(rc == TILEDB_OK);
  tiledb_array_free(&array);
  tiledb_query_free(&query);
  for (auto& validity_vec : validity_vecs)
    tiledb_validity_vec_free(&validity_vec);
}

void DenseNullableArrayFx::read(
    const string& array_name,
    const vector<test_query_buffer_t>& test_query_buffers,
    const void* const subarray,
    const tiledb_layout_t* const layout) {
  // Open the array for reading.
  tiledb_array_t* array;
  int rc =
      tiledb_array_alloc(ctx_, (FILE_TEMP_DIR + array_name).c_str(), &array);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  REQUIRE(rc == TILEDB_OK);

  // Create the read query.
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);

  // Set the query layout, if requested.
  if (layout != nullptr) {
    rc = tiledb_query_set_layout(ctx_, query, *layout);
    REQUIRE(rc == TILEDB_OK);
  }

  // Create a map from the index of each element in `test_query_buffers`
  // to a validity object that we read into.
  unordered_map<size_t, tiledb_validity_vec_t*> validity_vec_map;

  // Set the query buffers.
  for (size_t i = 0; i < test_query_buffers.size(); ++i) {
    const test_query_buffer_t& test_query_buffer = test_query_buffers[i];
    if (test_query_buffer.buffer_validity_size_ == nullptr) {
      rc = tiledb_query_set_buffer(
          ctx_,
          query,
          test_query_buffer.name_.c_str(),
          test_query_buffer.buffer_,
          test_query_buffer.buffer_size_);
      REQUIRE(rc == TILEDB_OK);
    } else {
      tiledb_validity_vec_t* validity_vec = nullptr;
      rc = tiledb_validity_vec_alloc(
          ctx_, *test_query_buffer.buffer_validity_size_, &validity_vec);
      REQUIRE(rc == TILEDB_OK);

      rc = tiledb_query_set_buffer_nullable(
          ctx_,
          query,
          test_query_buffer.name_.c_str(),
          test_query_buffer.buffer_,
          test_query_buffer.buffer_size_,
          validity_vec);
      REQUIRE(rc == TILEDB_OK);

      validity_vec_map[i] = validity_vec;
    }
  }

  // Set the subarray to read.
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);

  // Submit the query.
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Finalize the query, a no-op for non-global writes.
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Update the `query_buffer_t` objects that have a validity vector.
  for (const auto& kv : validity_vec_map) {
    const size_t i = kv.first;
    tiledb_validity_vec_t* const validity_vec = kv.second;

    // Update the query buffer.
    const test_query_buffer_t& test_query_buffer = test_query_buffers[i];
    rc = tiledb_validity_vec_size(
        ctx_, validity_vec, test_query_buffer.buffer_validity_size_);
    REQUIRE(rc == TILEDB_OK);
    rc = tiledb_validity_vec_values(
        ctx_,
        validity_vec,
        test_query_buffer.buffer_validity_,
        test_query_buffer.buffer_validity_size_);
    CHECK(rc == TILEDB_OK);
  }

  // Clean up
  rc = tiledb_array_close(ctx_, array);
  REQUIRE(rc == TILEDB_OK);
  tiledb_array_free(&array);
  tiledb_query_free(&query);
  for (auto& kv : validity_vec_map) {
    tiledb_validity_vec_t* validity_vec = kv.second;
    tiledb_validity_vec_free(&validity_vec);
  }
}

#if 0
TEST_CASE_METHOD(
    DenseNullableArrayFx,
    "C API: Test 1D dense array, full write, nullable, fixed attributes",
    "[1d][dense][full_write][nullable][fixed]") {
  const string array_name = "dense_nullable_array";

  // Define the dimensions.
  vector<test_dim_t> test_dims;
  const uint64_t d1_domain[] = {1, 4};
  const uint64_t d1_tile_extent = 2;
  test_dims.emplace_back("d1", TILEDB_UINT64, d1_domain, d1_tile_extent);

  // Define the attributes.
  vector<test_attr_t> test_attrs;
  test_attrs.emplace_back("a1", TILEDB_INT32, true);

  // Create the array.
  create_array(
      array_name, test_dims, test_attrs, TILEDB_ROW_MAJOR, TILEDB_ROW_MAJOR);

  // Define the write query buffers.
  vector<test_query_buffer_t> write_query_buffers;
  int a1_write_buffer[] = {1, 2, 3, 4};
  uint64_t a1_write_buffer_size = sizeof(a1_write_buffer);
  uint8_t a1_write_buffer_validity[] = {1, 1, 0, 1};
  uint64_t a1_write_buffer_validity_size = sizeof(a1_write_buffer_validity);
  write_query_buffers.emplace_back(
      "a1",
      a1_write_buffer,
      &a1_write_buffer_size,
      a1_write_buffer_validity,
      &a1_write_buffer_validity_size);

  // Execute the write query.
  write(array_name, write_query_buffers);

  // Define the read query buffers to read the entire domain.
  vector<test_query_buffer_t> read_query_buffers;
  int a1_read_buffer[4] = {0};
  uint64_t a1_read_buffer_size = sizeof(a1_read_buffer);
  uint8_t a1_read_buffer_validity[4] = {0};
  uint64_t a1_read_buffer_validity_size = sizeof(a1_read_buffer_validity);
  read_query_buffers.emplace_back(
      "a1",
      a1_read_buffer,
      &a1_read_buffer_size,
      a1_read_buffer_validity,
      &a1_read_buffer_validity_size);

  // Execute a read query over the entire domain.
  const uint64_t subarray_full[] = {1, 4};
  read(array_name, read_query_buffers, subarray_full);

  // Compare the fixed-value buffers.
  CHECK(a1_read_buffer_size == a1_write_buffer_size);
  CHECK(!memcmp(a1_read_buffer, a1_write_buffer, a1_read_buffer_size));

  // Compare the validity buffers.
  CHECK(a1_read_buffer_validity_size == a1_write_buffer_validity_size);
  CHECK(!memcmp(
      a1_read_buffer_validity,
      a1_write_buffer_validity,
      a1_read_buffer_validity_size));

  // Execute a read query over a partial domain.
  const uint64_t subarray_partial[] = {2, 3};
  read(array_name, read_query_buffers, subarray_partial);

  // Compare the fixed-value buffers.
  const int expected_a1_read_buffer[] = {2, 3};
  const uint64_t expected_a1_read_buffer_size = sizeof(expected_a1_read_buffer);
  CHECK(a1_read_buffer_size == expected_a1_read_buffer_size);
  CHECK(!memcmp(a1_read_buffer, expected_a1_read_buffer, a1_read_buffer_size));

  // Compare the validity buffers.
  const uint8_t expected_a1_read_buffer_validity[] = {1, 0};
  const uint64_t expected_a1_read_buffer_validity_size = sizeof(expected_a1_read_buffer_validity);
  CHECK(a1_read_buffer_validity_size == expected_a1_read_buffer_validity_size);
  CHECK(!memcmp(
      a1_read_buffer_validity,
      expected_a1_read_buffer_validity,
      a1_read_buffer_validity_size));
}
#endif

TEST_CASE_METHOD(
    DenseNullableArrayFx,
    "C API: Test 2D dense array, full write, nullable, fixed attributes",
    "[2d][dense][full_write][nullable][fixed]") {
  const string array_name = "dense_nullable_array";

  // Define the dimensions.
  vector<test_dim_t> test_dims;
  const uint64_t d1_domain[] = {1, 4};
  const uint64_t d1_tile_extent = 2;
  test_dims.emplace_back("d1", TILEDB_UINT64, d1_domain, d1_tile_extent);
  const uint64_t d2_domain[] = {1, 4};
  const uint64_t d2_tile_extent = 2;
  test_dims.emplace_back("d2", TILEDB_UINT64, d2_domain, d2_tile_extent);

  // Define the attributes.
  vector<test_attr_t> test_attrs;
  SECTION("1 attribute(s)") {
    test_attrs.emplace_back("a1", TILEDB_INT32, true);
  }
  SECTION("2 attribute(s)") {
    test_attrs.emplace_back("a1", TILEDB_INT32, true);
    test_attrs.emplace_back("a2", TILEDB_INT32, true);
  }

  // Create the array.
  create_array(
      array_name, test_dims, test_attrs, TILEDB_ROW_MAJOR, TILEDB_ROW_MAJOR);

  // Define the write query buffers for "a1".
  vector<test_query_buffer_t> write_query_buffers;
  int a1_write_buffer[] = {
      1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  uint64_t a1_write_buffer_size = sizeof(a1_write_buffer);
  uint8_t a1_write_buffer_validity[] = {
      1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
  uint64_t a1_write_buffer_validity_size = sizeof(a1_write_buffer_validity);
  write_query_buffers.emplace_back(
      "a1",
      a1_write_buffer,
      &a1_write_buffer_size,
      a1_write_buffer_validity,
      &a1_write_buffer_validity_size);

  // Define the write query buffers for "a2".
  int a2_write_buffer[16];
  uint64_t a2_write_buffer_size = sizeof(a2_write_buffer);
  for (uint64_t i = 0; i < 16; ++i)
    a2_write_buffer[i] = a1_write_buffer[16 - 1 - i];
  uint8_t a2_write_buffer_validity[16];
  uint64_t a2_write_buffer_validity_size = sizeof(a2_write_buffer_validity);
  for (uint64_t i = 0; i < 16; ++i)
    a2_write_buffer_validity[i] = a1_write_buffer_validity[16 - 1 - i];
  if (test_attrs.size() == 2) {
    write_query_buffers.emplace_back(
        "a2",
        a2_write_buffer,
        &a2_write_buffer_size,
        a2_write_buffer_validity,
        &a2_write_buffer_validity_size);
  }

  // Execute the write query.
  write(array_name, write_query_buffers);

  // Define the read query buffers for "a1".
  vector<test_query_buffer_t> read_query_buffers;
  int a1_read_buffer[16] = {0};
  uint64_t a1_read_buffer_size = sizeof(a1_read_buffer);
  uint8_t a1_read_buffer_validity[16] = {0};
  uint64_t a1_read_buffer_validity_size = sizeof(a1_read_buffer_validity);
  read_query_buffers.emplace_back(
      "a1",
      a1_read_buffer,
      &a1_read_buffer_size,
      a1_read_buffer_validity,
      &a1_read_buffer_validity_size);

  // Define the read query buffers for "a2".
  int a2_read_buffer[16] = {0};
  uint64_t a2_read_buffer_size = sizeof(a2_read_buffer);
  uint8_t a2_read_buffer_validity[16] = {0};
  uint64_t a2_read_buffer_validity_size = sizeof(a2_read_buffer_validity);
  if (test_attrs.size() == 2) {
    read_query_buffers.emplace_back(
        "a2",
        a2_read_buffer,
        &a2_read_buffer_size,
        a2_read_buffer_validity,
        &a2_read_buffer_validity_size);
  }

  // Execute a read query over the entire domain.
  const uint64_t subarray_full[] = {1, 4, 1, 4};
  read(array_name, read_query_buffers, subarray_full);

  // Compare the fixed-value buffers.
  CHECK(a1_read_buffer_size == a1_write_buffer_size);
  CHECK(!memcmp(a1_read_buffer, a1_write_buffer, a1_read_buffer_size));
  if (test_attrs.size() == 2) {
    CHECK(a2_read_buffer_size == a2_write_buffer_size);
    CHECK(!memcmp(a2_read_buffer, a2_write_buffer, a2_read_buffer_size));
  }

  // Compare the validity buffers.
  CHECK(a1_read_buffer_validity_size == a1_write_buffer_validity_size);
  CHECK(!memcmp(
      a1_read_buffer_validity,
      a1_write_buffer_validity,
      a1_read_buffer_validity_size));
  if (test_attrs.size() == 2) {
    CHECK(a2_read_buffer_validity_size == a2_write_buffer_validity_size);
    CHECK(!memcmp(
        a2_read_buffer_validity,
        a2_write_buffer_validity,
        a2_read_buffer_validity_size));
  }

  // Execute a read query over a partial domain.
  const uint64_t subarray_partial[] = {2, 3, 2, 3};
  read(array_name, read_query_buffers, subarray_partial);

  // Compare the fixed-value buffers.
  const int expected_a1_read_buffer[] = {6, 7, 10, 11};
  const uint64_t expected_a1_read_buffer_size = sizeof(expected_a1_read_buffer);
  CHECK(a1_read_buffer_size == expected_a1_read_buffer_size);
  CHECK(!memcmp(a1_read_buffer, expected_a1_read_buffer, a1_read_buffer_size));
  if (test_attrs.size() == 2) {
    const int expected_a2_read_buffer[] = {11, 10, 7, 6};
    const uint64_t expected_a2_read_buffer_size =
        sizeof(expected_a2_read_buffer);
    CHECK(a2_read_buffer_size == expected_a2_read_buffer_size);
    CHECK(
        !memcmp(a2_read_buffer, expected_a2_read_buffer, a2_read_buffer_size));
  }

  // Compare the validity buffers.
  const uint8_t expected_a1_read_buffer_validity[] = {1, 1, 0, 1};
  const uint64_t expected_a1_read_buffer_validity_size =
      sizeof(expected_a1_read_buffer_validity);
  CHECK(a1_read_buffer_validity_size == expected_a1_read_buffer_validity_size);
  CHECK(!memcmp(
      a1_read_buffer_validity,
      expected_a1_read_buffer_validity,
      a1_read_buffer_validity_size));
  if (test_attrs.size() == 2) {
    const uint8_t expected_a2_read_buffer_validity[] = {1, 0, 1, 1};
    const uint64_t expected_a2_read_buffer_validity_size =
        sizeof(expected_a2_read_buffer_validity);
    CHECK(
        a2_read_buffer_validity_size == expected_a2_read_buffer_validity_size);
    CHECK(!memcmp(
        a2_read_buffer_validity,
        expected_a2_read_buffer_validity,
        a2_read_buffer_validity_size));
  }
}

#if 0
TEST_CASE_METHOD(
    DenseNullableArrayFx,
    "C API: Test 3D dense array, full write, nullable, fixed attributes",
    "[3d][dense][full_write][nullable][fixed]") {
  const string array_name = "dense_nullable_array";

  // Define the dimensions.
  vector<test_dim_t> test_dims;
  const uint64_t d1_domain[] = {1, 4};
  const uint64_t d1_tile_extent = 2;
  test_dims.emplace_back("d1", TILEDB_UINT64, d1_domain, d1_tile_extent);
  const uint64_t d2_domain[] = {1, 4};
  const uint64_t d2_tile_extent = 2;
  test_dims.emplace_back("d2", TILEDB_UINT64, d2_domain, d2_tile_extent);
  const uint64_t d3_domain[] = {1, 4};
  const uint64_t d3_tile_extent = 2;
  test_dims.emplace_back("d3", TILEDB_UINT64, d3_domain, d3_tile_extent);

  // Define the attributes.
  vector<test_attr_t> test_attrs;
  test_attrs.emplace_back("a1", TILEDB_INT32, true);

  // Create the array.
  create_array(
      array_name, test_dims, test_attrs, TILEDB_ROW_MAJOR, TILEDB_ROW_MAJOR);

  // Define the write query buffers.
  vector<test_query_buffer_t> write_query_buffers;
  int a1_write_buffer[] = {
      1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
      17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
      33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
      49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64
  };
  uint64_t a1_write_buffer_size = sizeof(a1_write_buffer);
  uint8_t a1_write_buffer_validity[] = {
      1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
      1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
      1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
      1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0
  };
  uint64_t a1_write_buffer_validity_size = sizeof(a1_write_buffer_validity);
  write_query_buffers.emplace_back(
      "a1",
      a1_write_buffer,
      &a1_write_buffer_size,
      a1_write_buffer_validity,
      &a1_write_buffer_validity_size);

  // Execute the write query.
  write(array_name, write_query_buffers);

  // Define the read query buffers to read the entire domain.
  vector<test_query_buffer_t> read_query_buffers;
  int a1_read_buffer[64] = {0};
  uint64_t a1_read_buffer_size = sizeof(a1_read_buffer);
  uint8_t a1_read_buffer_validity[64] = {0};
  uint64_t a1_read_buffer_validity_size = sizeof(a1_read_buffer_validity);
  read_query_buffers.emplace_back(
      "a1",
      a1_read_buffer,
      &a1_read_buffer_size,
      a1_read_buffer_validity,
      &a1_read_buffer_validity_size);

  // Execute a read query over the entire domain.
  const uint64_t subarray_full[] = {1, 4, 1, 4, 1, 4};
  read(array_name, read_query_buffers, subarray_full);

  // Compare the fixed-value buffers.
  CHECK(a1_read_buffer_size == a1_write_buffer_size);
  CHECK(!memcmp(a1_read_buffer, a1_write_buffer, a1_read_buffer_size));

  // Compare the validity buffers.
  CHECK(a1_read_buffer_validity_size == a1_write_buffer_validity_size);
  CHECK(!memcmp(
      a1_read_buffer_validity,
      a1_write_buffer_validity,
      a1_read_buffer_validity_size));

  // Execute a read query over a partial domain.
  const uint64_t subarray_partial[] = {2, 3, 2, 3, 2, 3};
  read(array_name, read_query_buffers, subarray_partial);

  // Compare the fixed-value buffers.
  const int expected_a1_read_buffer[] = {22, 23, 26, 27, 38, 39, 42, 43};
  const uint64_t expected_a1_read_buffer_size = sizeof(expected_a1_read_buffer);
  CHECK(a1_read_buffer_size == expected_a1_read_buffer_size);
  CHECK(!memcmp(a1_read_buffer, expected_a1_read_buffer, a1_read_buffer_size));

  // Compare the validity buffers.
  const uint8_t expected_a1_read_buffer_validity[] = {1, 0, 0, 1, 1, 0, 0, 1};
  const uint64_t expected_a1_read_buffer_validity_size = sizeof(expected_a1_read_buffer_validity);
  CHECK(a1_read_buffer_validity_size == expected_a1_read_buffer_validity_size);
  CHECK(!memcmp(
      a1_read_buffer_validity,
      expected_a1_read_buffer_validity,
      a1_read_buffer_validity_size));
}
#endif

// 1d, 3d]
// var attr
// diff layouts