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

import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow


// Note: obviously, I could have defined a bundle of metrics like this in the
// S2 class itself rather than just for testing.  However, it's not clear that
// this is useful other than for testing purposes, and I find
//  S2CellMetrics.kMinWidth.GetLevelForMinValue(width) to be slightly more readable than
// than  S2CellMetrics.kWidth.min().GetLevelForMinValue(width).  Also, there is no
// fundamental reason that we need to analyze the minimum, maximum, and average
// values of every metric; it would be perfectly reasonable to just define
// one of these.
data class MetricBundle(val dim: Int, val min: S2CellMetric, val max: S2CellMetric, val avg: S2CellMetric)

class S2CellMetricsTest : S2GeometryTestCase() {

    fun checkMinMaxAvg(bundle: MetricBundle) {
        assertTrue(bundle.min.deriv <= bundle.avg.deriv)
        assertTrue(bundle.avg.deriv <= bundle.max.deriv)
    }

    fun checkLessOrEqual(a: MetricBundle, b: MetricBundle) {
        assertTrue(a.min.deriv <= b.min.deriv)
        assertTrue(a.max.deriv <= b.max.deriv)
        assertTrue(a.avg.deriv <= b.avg.deriv)
    }

    fun test() {
        val angleSpanBundle = MetricBundle(1, S2CellMetrics.kMinAngleSpan, S2CellMetrics.kMaxAngleSpan, S2CellMetrics.kAvgAngleSpan)
        val widthBundle = MetricBundle(1, S2CellMetrics.kMinWidth, S2CellMetrics.kMaxWidth, S2CellMetrics.kAvgWidth)
        val edgeBundle = MetricBundle(1, S2CellMetrics.kMinEdge, S2CellMetrics.kMaxEdge, S2CellMetrics.kAvgEdge)
        val diagBundle = MetricBundle(1, S2CellMetrics.kMinDiag, S2CellMetrics.kMaxDiag, S2CellMetrics.kAvgDiag)
        val areaBundle = MetricBundle(2, S2CellMetrics.kMinArea, S2CellMetrics.kMaxArea, S2CellMetrics.kAvgArea)

        // First, check that min <= avg <= max for each metric.
        checkMinMaxAvg(angleSpanBundle)
        checkMinMaxAvg(widthBundle)
        checkMinMaxAvg(edgeBundle)
        checkMinMaxAvg(diagBundle)
        checkMinMaxAvg(areaBundle)

        // Check that the maximum aspect ratio of an individual cell is consistent
        // with the global minimums and maximums.
        assertTrue(S2CellMetrics.kMaxEdgeAspect >= 1)
        assertTrue(S2CellMetrics.kMaxEdgeAspect <= S2CellMetrics.kMaxEdge.deriv / S2CellMetrics.kMinEdge.deriv)
        assertTrue(S2CellMetrics.kMaxDiagAspect >= 1)
        assertTrue(S2CellMetrics.kMaxDiagAspect <= S2CellMetrics.kMaxDiag.deriv / S2CellMetrics.kMinDiag.deriv)

        // Check various conditions that are provable mathematically.
        checkLessOrEqual(widthBundle, angleSpanBundle)
        checkLessOrEqual(widthBundle, edgeBundle)
        checkLessOrEqual(edgeBundle, diagBundle)

        assertTrue(S2CellMetrics.kMinArea.deriv >= S2CellMetrics.kMinWidth.deriv * S2CellMetrics.kMinEdge.deriv - 1e-15)
        assertTrue(S2CellMetrics.kMaxArea.deriv <= S2CellMetrics.kMaxWidth.deriv * S2CellMetrics.kMaxEdge.deriv + 1e-15)

        // GetLevelForMaxValue() and friends have built-in assertions, we just need
        // to call these functions to test them.
        //
        // We don't actually check that the metrics are correct here, e.g. that
        // GetMinWidth(10) is a lower bound on the width of cells at level 10.
        // It is easier to check these properties in s2cell_test, since
        // S2Cell has methods to compute the cell vertices, etc.

        for (level in -2..(S2Coords.kMaxCellLevel + 3)) {
            var width = S2CellMetrics.kMinWidth.deriv * 2.0.pow((-level).toDouble())
            if (level >= S2Coords.kMaxCellLevel + 3) width = 0.0
            // Check boundary cases (exactly equal to a threshold value).
            val expectedLevel = max(0, min(S2Coords.kMaxCellLevel, level))
            assertEquals(S2CellMetrics.kMinWidth.getLevelForMaxValue(width), expectedLevel)
            assertEquals(S2CellMetrics.kMinWidth.getLevelForMinValue(width), expectedLevel)
            assertEquals(S2CellMetrics.kMinWidth.getClosestLevel(width), expectedLevel)

            // Also check non-boundary cases.
            assertEquals(S2CellMetrics.kMinWidth.getLevelForMaxValue(1.2 * width), expectedLevel)
            assertEquals(S2CellMetrics.kMinWidth.getLevelForMinValue(0.8 * width), expectedLevel)
            assertEquals(S2CellMetrics.kMinWidth.getClosestLevel(1.2 * width), expectedLevel)
            assertEquals(S2CellMetrics.kMinWidth.getClosestLevel(0.8 * width), expectedLevel)

            // Same thing for area.
            var area = S2CellMetrics.kMinArea.deriv * 4.0.pow((-level).toDouble())
            if (level <= -3) area = 0.0
            assertEquals(S2CellMetrics.kMinArea.getLevelForMaxValue(area), expectedLevel)
            assertEquals(S2CellMetrics.kMinArea.getLevelForMinValue(area), expectedLevel)
            assertEquals(S2CellMetrics.kMinArea.getClosestLevel(area), expectedLevel)
            assertEquals(S2CellMetrics.kMinArea.getLevelForMaxValue(1.2 * area), expectedLevel)
            assertEquals(S2CellMetrics.kMinArea.getLevelForMinValue(0.8 * area), expectedLevel)
            assertEquals(S2CellMetrics.kMinArea.getClosestLevel(1.2 * area), expectedLevel)
            assertEquals(S2CellMetrics.kMinArea.getClosestLevel(0.8 * area), expectedLevel)
        }
    }
}
