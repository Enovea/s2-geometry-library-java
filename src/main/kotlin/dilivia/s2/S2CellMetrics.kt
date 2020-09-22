/**
 * This project is a kotlin port of the Google s2 geometry library (Copyright 2005 Google Inc. All Rights Reserved.):
 *                                 https://github.com/google/s2geometry.git
 *
 * Copyright © 2020 Dilivia (contact@dilivia.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package dilivia.s2

import com.google.common.geometry.S2.M_SQRT2
import dilivia.s2.Assertions.assert
import kotlin.math.max
import kotlin.math.min

/**
// The following are various constants that describe the shapes and sizes of
// S2Cells (see s2coords.h and s2cell_id.h).  They are useful for deciding
// which cell level to use in order to satisfy a given condition (e.g. that
// cell vertices must be no further than "x" apart).  All of the raw constants
// are differential quantities; you can use the GetValue(level) method to
// compute the corresponding length or area on the unit sphere for cells at a
// given level.  The minimum and maximum bounds are valid for cells at all
// levels, but they may be somewhat conservative for very large cells
// (e.g. face cells).
 *
 * Defines a cell metric of the given dimension (1 == length, 2 == area).
 *
 * @property dim The metric dimension (1 == length, 2 == area).
 * @property deriv The "deriv" value of a metric is a derivative, and must be multiplied by a length or area in
 * (s,t)-space to get a useful value.
*/
open class S2CellMetric(val dim: Int, val deriv: Double) {

  // Return the value of a metric for cells at the given level. The value is
  // either a length or an area on the unit sphere, depending on the
  // particular metric.
  fun getValue(level: Int): Double = StrictMath.scalb(deriv, - dim * level)

  // Return the level at which the metric has approximately the given
  // value.  For example, S2::kAvgEdge.GetClosestLevel(0.1) returns the
  // level at which the average cell edge length is approximately 0.1.
  // The return value is always a valid level.
  fun getClosestLevel(value: Double): Int = getLevelForMaxValue((if (dim == 1) M_SQRT2 else 2.0) * value)

  // Return the minimum level such that the metric is at most the given
  // value, or S2CellId::kMaxLevel if there is no such level.  For example,
  // S2::kMaxDiag.GetLevelForMaxValue(0.1) returns the minimum level such
  // that all cell diagonal lengths are 0.1 or smaller.  The return value
  // is always a valid level.
  fun getLevelForMaxValue(value: Double): Int {
      if (value <= 0) return S2Coords.kMaxCellLevel

      // This code is equivalent to computing a floating-point "level"
      // value and rounding up.  ilogb() returns the exponent corresponding to a
      // fraction in the range [1,2).
      val exponent = StrictMath.getExponent(value / deriv)
      val level = max(0, min(S2Coords.kMaxCellLevel, -(exponent shr (dim - 1))))
      assert({ level == S2Coords.kMaxCellLevel || getValue(level) <= value }, { "getLevelForMaxValue($value): deriv = $deriv => level = $level, getValue(level) = ${getValue(level)}" })
      assert { level == 0 || getValue(level - 1) > value }
      return level;
  }

  // Return the maximum level such that the metric is at least the given
  // value, or zero if there is no such level.  For example,
  // S2::kMinWidth.GetLevelForMinValue(0.1) returns the maximum level such
  // that all cells have a minimum width of 0.1 or larger.  The return value
  // is always a valid level.
  fun getLevelForMinValue(value: Double): Int {
      if (value <= 0) return S2Coords.kMaxCellLevel

      // This code is equivalent to computing a floating-point "level"
      // value and rounding down.
      val exponent = StrictMath.getExponent(deriv / value)
      val level = max(0, min(S2Coords.kMaxCellLevel, exponent shr (dim - 1)))
      assert { level == 0 || getValue(level) >= value }
      assert { level == S2Coords.kMaxCellLevel || getValue(level + 1) < value }
      return level;
  }

}

class LengthMetric(deriv: Double): S2CellMetric(1, deriv)
class AreaMetric(deriv: Double): S2CellMetric(2, deriv)

object S2CellMetrics {
    // Each cell is bounded by four planes passing through its four edges and
    // the center of the sphere.  These metrics relate to the angle between each
    // pair of opposite bounding planes, or equivalently, between the planes
    // corresponding to two different s-values or two different t-values.  For
    // example, the maximum angle between opposite bounding planes for a cell at
    // level k is kMaxAngleSpan.GetValue(k), and the average angle span for all
    // cells at level k is approximately kAvgAngleSpan.GetValue(k).
    val kMinAngleSpan: LengthMetric = S2Coords.projection.kMinAngleSpan
    val kMaxAngleSpan: LengthMetric = S2Coords.projection.kMaxAngleSpan
    val kAvgAngleSpan: LengthMetric = S2Coords.projection.kAvgAngleSpan

    // The width of geometric figure is defined as the distance between two
    // parallel bounding lines in a given direction.  For cells, the minimum
    // width is always attained between two opposite edges, and the maximum
    // width is attained between two opposite vertices.  However, for our
    // purposes we redefine the width of a cell as the perpendicular distance
    // between a pair of opposite edges.  A cell therefore has two widths, one
    // in each direction.  The minimum width according to this definition agrees
    // with the classic geometric one, but the maximum width is different.  (The
    // maximum geometric width corresponds to kMaxDiag defined below.)
    //
    // For a cell at level k, the distance between opposite edges is at least
    // kMinWidth.GetValue(k) and at most kMaxWidth.GetValue(k).  The average
    // width in both directions for all cells at level k is approximately
    // kAvgWidth.GetValue(k).
    //
    // The width is useful for bounding the minimum or maximum distance from a
    // point on one edge of a cell to the closest point on the opposite edge.
    // For example, this is useful when "growing" regions by a fixed distance.
    //
    // Note that because S2Cells are not usually rectangles, the minimum width of
    // a cell is generally smaller than its minimum edge length.  (The interior
    // angles of an S2Cell range from 60 to 120 degrees.)
    val kMinWidth: LengthMetric = S2Coords.projection.kMinWidth
    val kMaxWidth: LengthMetric = S2Coords.projection.kMaxWidth
    val kAvgWidth: LengthMetric = S2Coords.projection.kAvgWidth

    // The minimum edge length of any cell at level k is at least
    // kMinEdge.GetValue(k), and the maximum is at most kMaxEdge.GetValue(k).
    // The average edge length is approximately kAvgEdge.GetValue(k).
    //
    // The edge length metrics can also be used to bound the minimum, maximum,
    // or average distance from the center of one cell to the center of one of
    // its edge neighbors.  In particular, it can be used to bound the distance
    // between adjacent cell centers along the space-filling Hilbert curve for
    // cells at any given level.
    val kMinEdge:  LengthMetric = S2Coords.projection.kMinEdge
    val kMaxEdge:  LengthMetric = S2Coords.projection.kMaxEdge
    val kAvgEdge:  LengthMetric = S2Coords.projection.kAvgEdge

    // The minimum diagonal length of any cell at level k is at least
    // kMinDiag.GetValue(k), and the maximum is at most kMaxDiag.GetValue(k).
    // The average diagonal length is approximately kAvgDiag.GetValue(k).
    //
    // The maximum diagonal also happens to be the maximum diameter of any cell,
    // and also the maximum geometric width (see the discussion above).  So for
    // example, the distance from an arbitrary point to the closest cell center
    // at a given level is at most half the maximum diagonal length.
    val kMinDiag: LengthMetric = S2Coords.projection.kMinDiag
    val kMaxDiag: LengthMetric = S2Coords.projection.kMaxDiag
    val kAvgDiag: LengthMetric = S2Coords.projection.kAvgDiag

    // The minimum area of any cell at level k is at least kMinArea.GetValue(k),
    // and the maximum is at most kMaxArea.GetValue(k).  The average area of all
    // cells at level k is exactly kAvgArea.GetValue(k).
    val kMinArea:  AreaMetric = S2Coords.projection.kMinArea
    val kMaxArea:  AreaMetric = S2Coords.projection.kMaxArea
    val kAvgArea:  AreaMetric = S2Coords.projection.kAvgArea

    // This is the maximum edge aspect ratio over all cells at any level, where
    // the edge aspect ratio of a cell is defined as the ratio of its longest
    // edge length to its shortest edge length.
    val kMaxEdgeAspect: Double = S2Coords.projection.kMaxEdgeAspect

    // This is the maximum diagonal aspect ratio over all cells at any level,
    // where the diagonal aspect ratio of a cell is defined as the ratio of its
    // longest diagonal length to its shortest diagonal length.
    val kMaxDiagAspect: Double = S2Coords.projection.kMaxDiagAspect

}


