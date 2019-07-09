/**
 * @file unit-ResultCellSlabIter.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2019 TileDB, Inc.
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
 * Tests the `ResultCellSlabIter` class.
 */

#include "test/src/helpers.h"
#include "tiledb/sm/c_api/tiledb_struct_def.h"
#include "tiledb/sm/query/result_cell_slab_iter.h"

#ifdef _WIN32
#include "tiledb/sm/filesystem/win.h"
#else
#include "tiledb/sm/filesystem/posix.h"
#endif

#include <catch.hpp>
#include <iostream>

using namespace tiledb::sm;

/* ********************************* */
/*         STRUCT DEFINITION         */
/* ********************************* */

struct ResultCellSlabIterFx {
  tiledb_ctx_t* ctx_;
  tiledb_vfs_t* vfs_;
  bool s3_supported_, hdfs_supported_;
  std::string temp_dir_;
  const std::string s3_bucket_name_ =
      "s3://" + random_bucket_name("tiledb") + "/";
  std::string array_name_;
  const char* ARRAY_NAME = "result_cell_slab_iter";
  tiledb_array_t* array_ = nullptr;

  ResultCellSlabIterFx();
  ~ResultCellSlabIterFx();

  template <class T>
  void check_iter(
      ResultCellSlabIter<T>* iter,
      const std::vector<std::vector<uint64_t>>& c_result_cell_slabs);
  template <class T>
  void create_result_space_tiles(
      unsigned dim_num,
      Layout layout,
      const T* domain,
      const T* tile_extents,
      const std::vector<std::vector<T>>& domain_slices,
      const std::vector<std::vector<uint8_t>>& tile_coords,
      std::map<const T*, ResultSpaceTile<T>>* result_space_tiles);
};

ResultCellSlabIterFx::ResultCellSlabIterFx() {
  ctx_ = nullptr;
  vfs_ = nullptr;
  hdfs_supported_ = false;
  s3_supported_ = false;

  get_supported_fs(&s3_supported_, &hdfs_supported_);
  create_ctx_and_vfs(s3_supported_, &ctx_, &vfs_);
  create_s3_bucket(s3_bucket_name_, s3_supported_, ctx_, vfs_);

// Create temporary directory based on the supported filesystem
#ifdef _WIN32
  temp_dir_ = tiledb::sm::Win::current_dir() + "\\tiledb_test\\";
#else
  temp_dir_ = "file://" + tiledb::sm::Posix::current_dir() + "/tiledb_test/";
#endif
  if (s3_supported_)
    temp_dir_ = s3_bucket_name_ + "tiledb/test/";
  // Removing this for now as these tests take too long on HDFS
  // if (hdfs_supported_)
  //   temp_dir_ = "hdfs:///tiledb_test/";
  create_dir(temp_dir_, ctx_, vfs_);

  array_name_ = temp_dir_ + ARRAY_NAME;
  int rc = tiledb_array_alloc(ctx_, array_name_.c_str(), &array_);
  CHECK(rc == TILEDB_OK);
}

ResultCellSlabIterFx::~ResultCellSlabIterFx() {
  tiledb_array_free(&array_);
  remove_dir(temp_dir_, ctx_, vfs_);
  tiledb_ctx_free(&ctx_);
  tiledb_vfs_free(&vfs_);
}

template <class T>
void ResultCellSlabIterFx::check_iter(
    ResultCellSlabIter<T>* iter,
    const std::vector<std::vector<uint64_t>>& c_result_cell_slabs) {
  CHECK(iter->end());
  CHECK(iter->begin().ok());
  for (const auto& rcs : c_result_cell_slabs) {
    auto result_cell_slab = iter->result_cell_slab();
    // result_cell_slab.print();

    if (rcs[0] == UINT64_MAX) {
      CHECK(result_cell_slab.tile_ == nullptr);
    } else {
      CHECK(result_cell_slab.tile_ != nullptr);
      CHECK(result_cell_slab.tile_->frag_idx_ == rcs[0]);
      CHECK(result_cell_slab.tile_->tile_idx_ == rcs[1]);
    }
    CHECK(result_cell_slab.start_ == rcs[2]);
    CHECK(result_cell_slab.length_ == rcs[3]);
    CHECK(!iter->end());
    ++(*iter);
  }

  CHECK(iter->end());
}

template <class T>
void ResultCellSlabIterFx::create_result_space_tiles(
    unsigned dim_num,
    Layout layout,
    const T* domain,
    const T* tile_extents,
    const std::vector<std::vector<T>>& domain_slices,
    const std::vector<std::vector<uint8_t>>& tile_coords,
    std::map<const T*, ResultSpaceTile<T>>* result_space_tiles) {
  std::vector<TileDomain<T>> frag_tile_domains;
  for (size_t i = 0; i < domain_slices.size(); ++i) {
    frag_tile_domains.emplace_back(
        (unsigned)(domain_slices.size() - i),
        dim_num,
        domain,
        &(domain_slices[i][0]),
        tile_extents,
        layout);
  }
  TileDomain<T> array_tile_domain(
      UINT32_MAX, dim_num, domain, domain, tile_extents, layout);
  Reader::compute_result_space_tiles<T>(
      tile_coords, array_tile_domain, frag_tile_domains, result_space_tiles);
}

/* ********************************* */
/*                TESTS              */
/* ********************************* */

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Empty iterator",
    "[ResultCellSlabIter][empty]") {
  Subarray* subarray = nullptr;
  std::map<const int32_t*, ResultSpaceTile<int32_t>> result_space_tiles;
  std::vector<ResultCoords<int32_t>> result_coords;
  ResultCellSlabIter<int32_t> iter(
      subarray, &result_space_tiles, &result_coords);
  CHECK(iter.end());
  CHECK(iter.begin().ok());
  CHECK(iter.end());
}

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Test 1D slabs, 1 fragment, full overlap",
    "[ResultCellSlabIter][slabs][1d][1f][full_overlap]") {
  // Create array
  uint64_t domain[] = {1, 100};
  uint64_t tile_extent = 10;
  create_array(
      ctx_,
      array_name_,
      TILEDB_DENSE,
      {"d"},
      {TILEDB_UINT64},
      {domain},
      {&tile_extent},
      {"a", "b"},
      {TILEDB_INT32, TILEDB_INT32},
      {1, TILEDB_VAR_NUM},
      {::Compressor(TILEDB_FILTER_LZ4, -1),
       ::Compressor(TILEDB_FILTER_LZ4, -1)},
      TILEDB_ROW_MAJOR,
      TILEDB_ROW_MAJOR,
      2);

  // Create subarray
  open_array(ctx_, array_, TILEDB_READ);
  Subarray subarray;
  SubarrayRanges<uint64_t> ranges = {{5, 15}};
  Layout subarray_layout = Layout::ROW_MAJOR;
  create_subarray(array_->array_, ranges, subarray_layout, &subarray);
  subarray.compute_tile_coords<uint64_t>();

  // Create result space tiles
  std::vector<std::vector<uint64_t>> domain_slices = {{1, 100}};
  const auto& tile_coords = subarray.tile_coords();
  std::map<const uint64_t*, ResultSpaceTile<uint64_t>> result_space_tiles;
  create_result_space_tiles(
      1,
      subarray_layout,
      domain,
      &tile_extent,
      domain_slices,
      tile_coords,
      &result_space_tiles);

  // Check iterator
  std::vector<ResultCoords<uint64_t>> result_coords;
  ResultCellSlabIter<uint64_t> iter(
      &subarray, &result_space_tiles, &result_coords);
  std::vector<std::vector<uint64_t>> c_result_cell_slabs = {
      {1, 0, 4, 6},
      {1, 1, 0, 5},
  };
  check_iter<uint64_t>(&iter, c_result_cell_slabs);

  close_array(ctx_, array_);
}

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Test 1D slabs, 1 fragment, no overlap",
    "[ResultCellSlabIter][slabs][1d][1f][no_overlap]") {
  // Create array
  uint64_t domain[] = {1, 100};
  uint64_t tile_extent = 10;
  create_array(
      ctx_,
      array_name_,
      TILEDB_DENSE,
      {"d"},
      {TILEDB_UINT64},
      {domain},
      {&tile_extent},
      {"a", "b"},
      {TILEDB_INT32, TILEDB_INT32},
      {1, TILEDB_VAR_NUM},
      {::Compressor(TILEDB_FILTER_LZ4, -1),
       ::Compressor(TILEDB_FILTER_LZ4, -1)},
      TILEDB_ROW_MAJOR,
      TILEDB_ROW_MAJOR,
      2);

  // Create subarray
  open_array(ctx_, array_, TILEDB_READ);
  Subarray subarray;
  SubarrayRanges<uint64_t> ranges = {{5, 15}};
  Layout subarray_layout = Layout::ROW_MAJOR;
  create_subarray(array_->array_, ranges, subarray_layout, &subarray);
  subarray.compute_tile_coords<uint64_t>();

  // Create result space tiles
  std::vector<std::vector<uint64_t>> domain_slices = {{20, 30}};
  const auto& tile_coords = subarray.tile_coords();
  std::map<const uint64_t*, ResultSpaceTile<uint64_t>> result_space_tiles;
  create_result_space_tiles(
      1,
      subarray_layout,
      domain,
      &tile_extent,
      domain_slices,
      tile_coords,
      &result_space_tiles);

  // Check iterator
  std::vector<ResultCoords<uint64_t>> result_coords;
  ResultCellSlabIter<uint64_t> iter(
      &subarray, &result_space_tiles, &result_coords);
  std::vector<std::vector<uint64_t>> c_result_cell_slabs = {
      {UINT64_MAX, 0, 4, 6},
      {UINT64_MAX, 1, 0, 5},
  };
  check_iter<uint64_t>(&iter, c_result_cell_slabs);

  close_array(ctx_, array_);
}

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Test 1D slabs, 2 fragments",
    "[ResultCellSlabIter][slabs][1d][2f]") {
  // Create array
  uint64_t domain[] = {1, 100};
  uint64_t tile_extent = 10;
  create_array(
      ctx_,
      array_name_,
      TILEDB_DENSE,
      {"d"},
      {TILEDB_UINT64},
      {domain},
      {&tile_extent},
      {"a", "b"},
      {TILEDB_INT32, TILEDB_INT32},
      {1, TILEDB_VAR_NUM},
      {::Compressor(TILEDB_FILTER_LZ4, -1),
       ::Compressor(TILEDB_FILTER_LZ4, -1)},
      TILEDB_ROW_MAJOR,
      TILEDB_ROW_MAJOR,
      2);

  // Create subarray
  open_array(ctx_, array_, TILEDB_READ);
  Subarray subarray;
  SubarrayRanges<uint64_t> ranges = {{5, 15, 3, 5, 11, 14}};
  Layout subarray_layout = Layout::ROW_MAJOR;
  create_subarray(array_->array_, ranges, subarray_layout, &subarray);
  subarray.compute_tile_coords<uint64_t>();

  // Create result space tiles
  std::vector<std::vector<uint64_t>> domain_slices = {{5, 12}, {4, 15}};
  const auto& tile_coords = subarray.tile_coords();
  std::map<const uint64_t*, ResultSpaceTile<uint64_t>> result_space_tiles;
  create_result_space_tiles(
      1,
      subarray_layout,
      domain,
      &tile_extent,
      domain_slices,
      tile_coords,
      &result_space_tiles);

  // Check iterator
  std::vector<ResultCoords<uint64_t>> result_coords;
  ResultCellSlabIter<uint64_t> iter(
      &subarray, &result_space_tiles, &result_coords);
  std::vector<std::vector<uint64_t>> c_result_cell_slabs = {
      {2, 0, 4, 6},
      {2, 1, 0, 2},
      {1, 1, 2, 3},
      {UINT64_MAX, 0, 2, 1},
      {1, 0, 3, 1},
      {2, 0, 4, 1},
      {2, 1, 0, 2},
      {1, 1, 2, 2},
  };
  check_iter<uint64_t>(&iter, c_result_cell_slabs);

  close_array(ctx_, array_);
}

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Test 1D slabs, 1 dense fragment, 2 sparse fragments",
    "[ResultCellSlabIter][slabs][1d][1df][2sf]") {
  // Create array
  uint64_t domain[] = {1, 100};
  uint64_t tile_extent = 10;
  create_array(
      ctx_,
      array_name_,
      TILEDB_DENSE,
      {"d"},
      {TILEDB_UINT64},
      {domain},
      {&tile_extent},
      {"a", "b"},
      {TILEDB_INT32, TILEDB_INT32},
      {1, TILEDB_VAR_NUM},
      {::Compressor(TILEDB_FILTER_LZ4, -1),
       ::Compressor(TILEDB_FILTER_LZ4, -1)},
      TILEDB_ROW_MAJOR,
      TILEDB_ROW_MAJOR,
      2);

  // Create subarray
  open_array(ctx_, array_, TILEDB_READ);
  Subarray subarray;
  SubarrayRanges<uint64_t> ranges = {{3, 15, 18, 20}};
  Layout subarray_layout = Layout::ROW_MAJOR;
  create_subarray(array_->array_, ranges, subarray_layout, &subarray);
  subarray.compute_tile_coords<uint64_t>();

  // Create result space tiles
  std::vector<std::vector<uint64_t>> domain_slices = {{3, 12}};
  const auto& tile_coords = subarray.tile_coords();
  std::map<const uint64_t*, ResultSpaceTile<uint64_t>> result_space_tiles;
  create_result_space_tiles(
      1,
      subarray_layout,
      domain,
      &tile_extent,
      domain_slices,
      tile_coords,
      &result_space_tiles);

  // Create result coordinates
  std::vector<ResultCoords<uint64_t>> result_coords;
  ResultTile result_tile_2_0(2, 0);
  ResultTile result_tile_3_0(3, 0);
  ResultTile result_tile_3_1(3, 1);
  uint64_t coords_2_3 = 3;
  uint64_t coords_2_5 = 5;
  uint64_t coords_3_8 = 8;
  uint64_t coords_3_9 = 9;
  uint64_t coords_3_12 = 12;
  uint64_t coords_3_19 = 19;
  result_coords.emplace_back(&result_tile_2_0, &coords_2_3, 1);
  result_coords.emplace_back(&result_tile_2_0, &coords_2_5, 3);
  result_coords.emplace_back(&result_tile_3_0, &coords_3_8, 2);
  result_coords.emplace_back(&result_tile_3_0, &coords_3_9, 3);
  result_coords.back().invalidate();
  result_coords.emplace_back(&result_tile_3_1, &coords_3_12, 1);
  result_coords.emplace_back(&result_tile_3_1, &coords_3_19, 2);

  // Check iterator
  ResultCellSlabIter<uint64_t> iter(
      &subarray, &result_space_tiles, &result_coords);
  std::vector<std::vector<uint64_t>> c_result_cell_slabs = {
      {2, 0, 1, 1},
      {1, 0, 3, 1},
      {2, 0, 3, 1},
      {1, 0, 5, 2},
      {3, 0, 2, 1},
      {1, 0, 8, 2},
      {1, 1, 0, 1},
      {3, 1, 1, 1},
      {UINT64_MAX, 1, 2, 3},
      {UINT64_MAX, 1, 7, 1},
      {3, 1, 2, 1},
      {UINT64_MAX, 1, 9, 1},
  };
  check_iter<uint64_t>(&iter, c_result_cell_slabs);

  close_array(ctx_, array_);
}

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Test 2D slabs, 1 range, 1 dense fragment, full "
    "overlap",
    "[ResultCellSlabIter][slabs][2d][1r][1f][full_overlap]") {
  Layout subarray_layout = Layout::ROW_MAJOR;
  Layout tile_domain_layout = Layout::ROW_MAJOR;
  tiledb_layout_t tile_order = TILEDB_ROW_MAJOR;
  tiledb_layout_t cell_order = TILEDB_ROW_MAJOR;
  std::vector<std::vector<uint64_t>> c_result_cell_slabs = {};

  SECTION("- tile: row, cell: row, subarray: row") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {1, 0, 4, 2},
        {1, 1, 3, 3},
        {1, 0, 7, 2},
        {1, 1, 6, 3},
    };
  }

  SECTION("- tile: row, cell: col, subarray: row") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {1, 0, 4, 2},
        {1, 1, 1, 3},
        {1, 0, 5, 2},
        {1, 1, 2, 3},
    };
  }

  SECTION("- tile: col, cell: row, subarray: row") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {1, 0, 4, 2},
        {1, 2, 3, 3},
        {1, 0, 7, 2},
        {1, 2, 6, 3},
    };
  }

  SECTION("- tile: col, cell: col, subarray: row") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {1, 0, 4, 2},
        {1, 2, 1, 3},
        {1, 0, 5, 2},
        {1, 2, 2, 3},
    };
  }

  SECTION("- tile: row, cell: row, subarray: col") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {1, 0, 4, 2},
        {1, 0, 5, 2},
        {1, 1, 3, 2},
        {1, 1, 4, 2},
        {1, 1, 5, 2},
    };
  }

  SECTION("- tile: row, cell: col, subarray: col") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {1, 0, 4, 2},
        {1, 0, 7, 2},
        {1, 1, 1, 2},
        {1, 1, 4, 2},
        {1, 1, 7, 2},
    };
  }

  SECTION("- tile: col, cell: row, subarray: col") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {1, 0, 4, 2},
        {1, 0, 5, 2},
        {1, 2, 3, 2},
        {1, 2, 4, 2},
        {1, 2, 5, 2},
    };
  }

  SECTION("- tile: col, cell: col, subarray: col") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {1, 0, 4, 2},
        {1, 0, 7, 2},
        {1, 2, 1, 2},
        {1, 2, 4, 2},
        {1, 2, 7, 2},
    };
  }

  // Create array
  uint64_t domain[] = {1, 6, 1, 6};
  uint64_t tile_extents[] = {3, 3};
  create_array(
      ctx_,
      array_name_,
      TILEDB_DENSE,
      {"d1", "d2"},
      {TILEDB_UINT64, TILEDB_UINT64},
      {&domain[0], &domain[2]},
      {&tile_extents[0], &tile_extents[1]},
      {"a", "b"},
      {TILEDB_INT32, TILEDB_INT32},
      {1, TILEDB_VAR_NUM},
      {::Compressor(TILEDB_FILTER_LZ4, -1),
       ::Compressor(TILEDB_FILTER_LZ4, -1)},
      tile_order,
      cell_order,
      2);

  // Create subarray
  open_array(ctx_, array_, TILEDB_READ);
  Subarray subarray;
  SubarrayRanges<uint64_t> ranges = {{2, 3}, {2, 6}};
  create_subarray(array_->array_, ranges, subarray_layout, &subarray);
  subarray.compute_tile_coords<uint64_t>();

  // Create result space tiles
  std::vector<std::vector<uint64_t>> domain_slices = {{1, 6, 1, 6}};
  const auto& tile_coords = subarray.tile_coords();
  std::map<const uint64_t*, ResultSpaceTile<uint64_t>> result_space_tiles;
  create_result_space_tiles(
      2,
      tile_domain_layout,
      domain,
      tile_extents,
      domain_slices,
      tile_coords,
      &result_space_tiles);

  // Create result coordinates
  std::vector<ResultCoords<uint64_t>> result_coords;

  // Check iterator
  ResultCellSlabIter<uint64_t> iter(
      &subarray, &result_space_tiles, &result_coords);
  check_iter<uint64_t>(&iter, c_result_cell_slabs);

  close_array(ctx_, array_);
}

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Test 2D slabs, 1 range, 1 dense fragment, no overlap",
    "[ResultCellSlabIter][slabs][2d][1r][1f][no_overlap]") {
  Layout subarray_layout = Layout::ROW_MAJOR;
  Layout tile_domain_layout = Layout::ROW_MAJOR;
  tiledb_layout_t tile_order = TILEDB_ROW_MAJOR;
  tiledb_layout_t cell_order = TILEDB_ROW_MAJOR;
  std::vector<std::vector<uint64_t>> c_result_cell_slabs = {};

  SECTION("- tile: row, cell: row, subarray: row") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 1, 3, 3},
        {UINT64_MAX, 0, 7, 2},
        {UINT64_MAX, 1, 6, 3},
    };
  }

  SECTION("- tile: row, cell: col, subarray: row") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 1, 1, 3},
        {UINT64_MAX, 0, 5, 2},
        {UINT64_MAX, 1, 2, 3},
    };
  }

  SECTION("- tile: col, cell: row, subarray: row") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 2, 3, 3},
        {UINT64_MAX, 0, 7, 2},
        {UINT64_MAX, 2, 6, 3},
    };
  }

  SECTION("- tile: col, cell: col, subarray: row") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 2, 1, 3},
        {UINT64_MAX, 0, 5, 2},
        {UINT64_MAX, 2, 2, 3},
    };
  }

  SECTION("- tile: row, cell: row, subarray: col") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 0, 5, 2},
        {UINT64_MAX, 1, 3, 2},
        {UINT64_MAX, 1, 4, 2},
        {UINT64_MAX, 1, 5, 2},
    };
  }

  SECTION("- tile: row, cell: col, subarray: col") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 0, 7, 2},
        {UINT64_MAX, 1, 1, 2},
        {UINT64_MAX, 1, 4, 2},
        {UINT64_MAX, 1, 7, 2},
    };
  }

  SECTION("- tile: col, cell: row, subarray: col") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 0, 5, 2},
        {UINT64_MAX, 2, 3, 2},
        {UINT64_MAX, 2, 4, 2},
        {UINT64_MAX, 2, 5, 2},
    };
  }

  SECTION("- tile: col, cell: col, subarray: col") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 0, 7, 2},
        {UINT64_MAX, 2, 1, 2},
        {UINT64_MAX, 2, 4, 2},
        {UINT64_MAX, 2, 7, 2},
    };
  }

  // Create array
  uint64_t domain[] = {1, 6, 1, 6};
  uint64_t tile_extents[] = {3, 3};
  create_array(
      ctx_,
      array_name_,
      TILEDB_DENSE,
      {"d1", "d2"},
      {TILEDB_UINT64, TILEDB_UINT64},
      {&domain[0], &domain[2]},
      {&tile_extents[0], &tile_extents[1]},
      {"a", "b"},
      {TILEDB_INT32, TILEDB_INT32},
      {1, TILEDB_VAR_NUM},
      {::Compressor(TILEDB_FILTER_LZ4, -1),
       ::Compressor(TILEDB_FILTER_LZ4, -1)},
      tile_order,
      cell_order,
      2);

  // Create subarray
  open_array(ctx_, array_, TILEDB_READ);
  Subarray subarray;
  SubarrayRanges<uint64_t> ranges = {{2, 3}, {2, 6}};
  create_subarray(array_->array_, ranges, subarray_layout, &subarray);
  subarray.compute_tile_coords<uint64_t>();

  // Create result space tiles
  std::vector<std::vector<uint64_t>> domain_slices = {{6, 6, 6, 6}};
  const auto& tile_coords = subarray.tile_coords();
  std::map<const uint64_t*, ResultSpaceTile<uint64_t>> result_space_tiles;
  create_result_space_tiles(
      2,
      tile_domain_layout,
      domain,
      tile_extents,
      domain_slices,
      tile_coords,
      &result_space_tiles);

  // Create result coordinates
  std::vector<ResultCoords<uint64_t>> result_coords;

  // Check iterator
  ResultCellSlabIter<uint64_t> iter(
      &subarray, &result_space_tiles, &result_coords);
  check_iter<uint64_t>(&iter, c_result_cell_slabs);

  close_array(ctx_, array_);
}

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Test 2D slabs, 1 range, 1 dense fragment, partial "
    "overlap",
    "[ResultCellSlabIter][slabs][2d][1r][1f][partial_overlap]") {
  Layout subarray_layout = Layout::ROW_MAJOR;
  Layout tile_domain_layout = Layout::ROW_MAJOR;
  tiledb_layout_t tile_order = TILEDB_ROW_MAJOR;
  tiledb_layout_t cell_order = TILEDB_ROW_MAJOR;
  std::vector<std::vector<uint64_t>> c_result_cell_slabs = {};

  SECTION("- tile: row, cell: row, subarray: row") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 1, 3, 3},
        {UINT64_MAX, 0, 7, 2},
        {UINT64_MAX, 1, 6, 1},
        {1, 0, 7, 2},
    };
  }

  SECTION("- tile: row, cell: col, subarray: row") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 1, 1, 3},
        {UINT64_MAX, 0, 5, 2},
        {UINT64_MAX, 1, 2, 1},
        {1, 0, 5, 2},
    };
  }

  SECTION("- tile: col, cell: row, subarray: row") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 2, 3, 3},
        {UINT64_MAX, 0, 7, 2},
        {UINT64_MAX, 2, 6, 1},
        {1, 0, 7, 2},
    };
  }

  SECTION("- tile: col, cell: col, subarray: row") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 2, 1, 3},
        {UINT64_MAX, 0, 5, 2},
        {UINT64_MAX, 2, 2, 1},
        {1, 0, 5, 2},
    };
  }

  SECTION("- tile: row, cell: row, subarray: col") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 0, 5, 2},
        {UINT64_MAX, 1, 3, 2},
        {UINT64_MAX, 1, 4, 1},
        {1, 0, 7, 1},
        {UINT64_MAX, 1, 5, 1},
        {1, 0, 8, 1},
    };
  }

  SECTION("- tile: row, cell: col, subarray: col") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 0, 7, 2},
        {UINT64_MAX, 1, 1, 2},
        {UINT64_MAX, 1, 4, 1},
        {1, 0, 5, 1},
        {UINT64_MAX, 1, 7, 1},
        {1, 0, 8, 1},
    };
  }

  SECTION("- tile: col, cell: row, subarray: col") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 0, 5, 2},
        {UINT64_MAX, 2, 3, 2},
        {UINT64_MAX, 2, 4, 1},
        {1, 0, 7, 1},
        {UINT64_MAX, 2, 5, 1},
        {1, 0, 8, 1},
    };
  }

  SECTION("- tile: col, cell: col, subarray: col") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {UINT64_MAX, 0, 4, 2},
        {UINT64_MAX, 0, 7, 2},
        {UINT64_MAX, 2, 1, 2},
        {UINT64_MAX, 2, 4, 1},
        {1, 0, 5, 1},
        {UINT64_MAX, 2, 7, 1},
        {1, 0, 8, 1},
    };
  }

  // Create array
  uint64_t domain[] = {1, 6, 1, 6};
  uint64_t tile_extents[] = {3, 3};
  create_array(
      ctx_,
      array_name_,
      TILEDB_DENSE,
      {"d1", "d2"},
      {TILEDB_UINT64, TILEDB_UINT64},
      {&domain[0], &domain[2]},
      {&tile_extents[0], &tile_extents[1]},
      {"a", "b"},
      {TILEDB_INT32, TILEDB_INT32},
      {1, TILEDB_VAR_NUM},
      {::Compressor(TILEDB_FILTER_LZ4, -1),
       ::Compressor(TILEDB_FILTER_LZ4, -1)},
      tile_order,
      cell_order,
      2);

  // Create subarray
  open_array(ctx_, array_, TILEDB_READ);
  Subarray subarray;
  SubarrayRanges<uint64_t> ranges = {{2, 3}, {2, 6}};
  create_subarray(array_->array_, ranges, subarray_layout, &subarray);
  subarray.compute_tile_coords<uint64_t>();

  // Create result space tiles
  std::vector<std::vector<uint64_t>> domain_slices = {{3, 6, 5, 6}};
  const auto& tile_coords = subarray.tile_coords();
  std::map<const uint64_t*, ResultSpaceTile<uint64_t>> result_space_tiles;
  create_result_space_tiles(
      2,
      tile_domain_layout,
      domain,
      tile_extents,
      domain_slices,
      tile_coords,
      &result_space_tiles);

  // Create result coordinates
  std::vector<ResultCoords<uint64_t>> result_coords;

  // Check iterator
  ResultCellSlabIter<uint64_t> iter(
      &subarray, &result_space_tiles, &result_coords);
  check_iter<uint64_t>(&iter, c_result_cell_slabs);

  close_array(ctx_, array_);
}

TEST_CASE_METHOD(
    ResultCellSlabIterFx,
    "ResultCellSlabIter: Test 2D slabs, multiple ranges, 2 denses fragments, 1 "
    "sparse "
    "overlap",
    "[ResultCellSlabIter][slabs][2d][mr][2df1sf]") {
  Layout subarray_layout = Layout::ROW_MAJOR;
  Layout tile_domain_layout = Layout::ROW_MAJOR;
  tiledb_layout_t tile_order = TILEDB_ROW_MAJOR;
  tiledb_layout_t cell_order = TILEDB_ROW_MAJOR;
  std::vector<std::vector<uint64_t>> c_result_cell_slabs = {};

  SECTION("- tile: row, cell: row, subarray: row") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {2, 0, 7, 1},
        {3, 0, 1, 1},
        {2, 1, 6, 1},
        {1, 1, 7, 2},
        {2, 2, 1, 2},
        {2, 3, 0, 1},
        {UINT64_MAX, 3, 1, 2},
        {2, 2, 4, 2},
        {2, 3, 3, 1},
        {3, 1, 0, 1},
        {3, 1, 2, 1},
    };
  }

  SECTION("- tile: row, cell: col, subarray: row") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {2, 0, 5, 1},
        {3, 0, 1, 1},
        {2, 1, 2, 1},
        {1, 1, 5, 2},
        {2, 2, 3, 2},
        {2, 3, 0, 1},
        {UINT64_MAX, 3, 3, 2},
        {2, 2, 4, 2},
        {2, 3, 1, 1},
        {3, 1, 0, 1},
        {3, 1, 2, 1},
    };
  }

  SECTION("- tile: col, cell: row, subarray: row") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {2, 0, 7, 1},
        {3, 0, 1, 1},
        {2, 2, 6, 1},
        {1, 1, 7, 2},
        {2, 1, 1, 2},
        {2, 3, 0, 1},
        {UINT64_MAX, 3, 1, 2},
        {2, 1, 4, 2},
        {2, 3, 3, 1},
        {3, 1, 0, 1},
        {3, 1, 2, 1},
    };
  }

  SECTION("- tile: col, cell: col, subarray: row") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::ROW_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {2, 0, 5, 1},
        {3, 0, 1, 1},
        {2, 2, 2, 1},
        {1, 1, 5, 2},
        {2, 1, 3, 2},
        {2, 3, 0, 1},
        {UINT64_MAX, 3, 3, 2},
        {2, 1, 4, 2},
        {2, 3, 1, 1},
        {3, 1, 0, 1},
        {3, 1, 2, 1},
    };
  }

  SECTION("- tile: row, cell: row, subarray: col") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {2, 0, 7, 1},
        {2, 2, 1, 2},
        {3, 0, 1, 1},
        {2, 2, 2, 2},
        {2, 1, 6, 1},
        {2, 3, 0, 2},
        {1, 1, 7, 1},
        {UINT64_MAX, 3, 1, 1},
        {3, 1, 0, 1},
        {1, 1, 8, 1},
        {UINT64_MAX, 3, 2, 1},
        {3, 1, 2, 1},
    };
  }

  SECTION("- tile: row, cell: col, subarray: col") {
    tile_order = TILEDB_ROW_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::ROW_MAJOR;
    c_result_cell_slabs = {
        {2, 0, 5, 1},
        {2, 2, 3, 2},
        {3, 0, 1, 1},
        {2, 2, 6, 2},
        {2, 1, 2, 1},
        {2, 3, 0, 2},
        {1, 1, 5, 1},
        {UINT64_MAX, 3, 3, 1},
        {3, 1, 0, 1},
        {1, 1, 8, 1},
        {UINT64_MAX, 3, 6, 1},
        {3, 1, 2, 1},
    };
  }

  SECTION("- tile: col, cell: row, subarray: col") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_ROW_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {2, 0, 7, 1},
        {2, 1, 1, 2},
        {3, 0, 1, 1},
        {2, 1, 2, 2},
        {2, 2, 6, 1},
        {2, 3, 0, 2},
        {1, 1, 7, 1},
        {UINT64_MAX, 3, 1, 1},
        {3, 1, 0, 1},
        {1, 1, 8, 1},
        {UINT64_MAX, 3, 2, 1},
        {3, 1, 2, 1},
    };
  }

  SECTION("- tile: col, cell: col, subarray: col") {
    tile_order = TILEDB_COL_MAJOR;
    cell_order = TILEDB_COL_MAJOR;
    subarray_layout = Layout::COL_MAJOR;
    tile_domain_layout = Layout::COL_MAJOR;
    c_result_cell_slabs = {
        {2, 0, 5, 1},
        {2, 1, 3, 2},
        {3, 0, 1, 1},
        {2, 1, 6, 2},
        {2, 2, 2, 1},
        {2, 3, 0, 2},
        {1, 1, 5, 1},
        {UINT64_MAX, 3, 3, 1},
        {3, 1, 0, 1},
        {1, 1, 8, 1},
        {UINT64_MAX, 3, 6, 1},
        {3, 1, 2, 1},
    };
  }

  // Create array
  uint64_t domain[] = {1, 6, 1, 6};
  uint64_t tile_extents[] = {3, 3};
  create_array(
      ctx_,
      array_name_,
      TILEDB_DENSE,
      {"d1", "d2"},
      {TILEDB_UINT64, TILEDB_UINT64},
      {&domain[0], &domain[2]},
      {&tile_extents[0], &tile_extents[1]},
      {"a", "b"},
      {TILEDB_INT32, TILEDB_INT32},
      {1, TILEDB_VAR_NUM},
      {::Compressor(TILEDB_FILTER_LZ4, -1),
       ::Compressor(TILEDB_FILTER_LZ4, -1)},
      tile_order,
      cell_order,
      2);

  // Create subarray
  open_array(ctx_, array_, TILEDB_READ);
  Subarray subarray;
  SubarrayRanges<uint64_t> ranges = {{3, 5}, {2, 4, 5, 6}};
  create_subarray(array_->array_, ranges, subarray_layout, &subarray);
  subarray.compute_tile_coords<uint64_t>();

  // Create result space tiles
  std::vector<std::vector<uint64_t>> domain_slices = {{3, 5, 2, 4},
                                                      {2, 3, 1, 6}};
  const auto& tile_coords = subarray.tile_coords();
  std::map<const uint64_t*, ResultSpaceTile<uint64_t>> result_space_tiles;
  create_result_space_tiles(
      2,
      tile_domain_layout,
      domain,
      tile_extents,
      domain_slices,
      tile_coords,
      &result_space_tiles);

  // Create result coordinates
  std::vector<ResultCoords<uint64_t>> result_coords;
  ResultTile result_tile_3_0(3, 0);
  ResultTile result_tile_3_1(3, 1);
  uint64_t coords_3_3_3[] = {3, 3};
  uint64_t coords_3_5_5[] = {5, 5};
  uint64_t coords_3_5_6[] = {5, 6};
  result_coords.emplace_back(&result_tile_3_0, coords_3_3_3, 1);
  result_coords.emplace_back(&result_tile_3_1, coords_3_5_5, 0);
  result_coords.emplace_back(&result_tile_3_1, coords_3_5_6, 2);

  // Check iterator
  ResultCellSlabIter<uint64_t> iter(
      &subarray, &result_space_tiles, &result_coords);
  check_iter<uint64_t>(&iter, c_result_cell_slabs);

  close_array(ctx_, array_);
}