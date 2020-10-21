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
package dilivia.s2.index

import com.google.common.collect.ComparisonChain
import dilivia.s2.S2.M_PI
import dilivia.s2.S2.M_SQRT1_2
import dilivia.s2.Assertions
import dilivia.s2.S1Angle
import dilivia.s2.S2CellId
import dilivia.s2.S2CellMetrics
import dilivia.s2.S2EdgeDistances
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2LatLng
import dilivia.s2.S2Measures
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.builder.snap.IntLatLngSnapFunction
import dilivia.s2.builder.snap.S2CellIdSnapFunction
import dilivia.s2.builder.snap.SnapFunction
import dilivia.s2.region.S2Cell
import java.math.BigDecimal
import java.math.MathContext
import java.math.RoundingMode
import kotlin.math.IEEErem
import kotlin.math.abs
import kotlin.math.min
import kotlin.math.roundToLong


// A scaled S2LatLng with integer coordinates, similar to E7 coordinates,
// except that the scale is variable (see LatLngConfig below).
typealias IntLatLng = Pair<Long, Long>

operator fun IntLatLng.get(idx: Int): Long = when (idx) {
    0 -> this.first
    1 -> this.second
    else -> throw ArrayIndexOutOfBoundsException(idx)
}

operator fun IntLatLng.plus(other: IntLatLng) = IntLatLng(this.first + other.first, this.second + other.second)

//
// The bulk of this file consists of "tests" that attempt to construct worst
// cases for the various constants used in S2CellIdSnapFunction and
// IntLatLngSnapFunction implementations.  For all of these constants I have
// done hand analysis of the planar configurations, but sometimes the
// spherical case is slightly better or worse because of the spherical
// distortion.

class S2SnapFunctionsTest : S2GeometryTestCase() {

    fun testS2CellIdSnapFunctionLevelToFromSnapRadius() {
        for (level in 0..S2CellId.kMaxLevel) {
            val radius = S2CellIdSnapFunction.minSnapRadiusForLevel(level)
            assertEquals(level, S2CellIdSnapFunction.levelForMaxSnapRadius(radius))
            assertEquals(min(level + 1, S2CellId.kMaxLevel), S2CellIdSnapFunction.levelForMaxSnapRadius(radius * 0.999))
        }
        assertEquals(0, S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.radians(5)))
        assertEquals(S2CellId.kMaxLevel, S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.radians(1e-30)))
    }

    fun testS2CellIdSnapFunctionSnapPoint() {
        repeat(1000) {
            for (level in 0..S2CellId.kMaxLevel) {
                // This checks that points are snapped to the correct level, since
                // S2CellId centers at different levels are always different.
                val f = S2CellIdSnapFunction(level)
                val p = S2Random.randomCellId(level).toPoint()
                assertEquals(p, f.snapPoint(p))
            }
        }
    }

    fun testIntLatLngSnapFunctionExponentToFromSnapRadius() {
        for (exponent in IntLatLngSnapFunction.kMinExponent..IntLatLngSnapFunction.kMaxExponent) {
            val radius = IntLatLngSnapFunction.minSnapRadiusForExponent(exponent)
            assertEquals(exponent, IntLatLngSnapFunction.exponentForMaxSnapRadius(radius))
            assertEquals(min(exponent + 1, IntLatLngSnapFunction.kMaxExponent), IntLatLngSnapFunction.exponentForMaxSnapRadius(radius * 0.999))
        }
        assertEquals(IntLatLngSnapFunction.kMinExponent, IntLatLngSnapFunction.exponentForMaxSnapRadius(S1Angle.radians(5)))
        assertEquals(IntLatLngSnapFunction.kMaxExponent, IntLatLngSnapFunction.exponentForMaxSnapRadius(S1Angle.radians(1e-30)))
    }

    fun testIntLatLngSnapFunctionSnapPoint() {
        repeat(1000) {
            // Test that IntLatLngSnapFunction does not modify points that were
            // generated using the S2LatLng::From{E5,E6,E7} methods.  This ensures
            // that both functions are using bitwise-compatible conversion methods.
            val p = S2Random.randomPoint()
            val ll = S2LatLng.fromPoint(p)
            val p5 = S2LatLng.fromE5(ll.lat().e5(), ll.lng().e5()).toPoint()
            assertEquals(p5, IntLatLngSnapFunction(5).snapPoint(p5))
            val p6 = S2LatLng.fromE6(ll.lat().e6(), ll.lng().e6()).toPoint()
            assertEquals(p6, IntLatLngSnapFunction(6).snapPoint(p6))
            val p7 = S2LatLng.fromE7(ll.lat().e7(), ll.lng().e7()).toPoint()
            assertEquals(p7, IntLatLngSnapFunction(7).snapPoint(p7))

            // Make sure that we're not snapping using some lower exponent.
            val p7not6 = S2LatLng.fromE7(10 * ll.lat().e6() + 1, 10 * ll.lng().e6() + 1).toPoint()
            assertTrue(p7not6 != IntLatLngSnapFunction(6).snapPoint(p7not6))
        }
    }


    fun testS2CellIdSnapFunctionMinVertexSeparationSnapRadiusRatio() {
        // The purpose of this "test" is to compute a lower bound to the fraction
        // (min_vertex_separation() / snap_radius()).  Essentially this involves
        // searching for two adjacent cells A and B such when one of the corner
        // vertices of B is snapped to the center of B, the distance to the center
        // of A decreases as much as possible.  In other words, we want the ratio
        //
        //   distance(center(A), center(B)) / distance(center(A), vertex(B))
        //
        // to be as small as possible.  We do this by considering one cell level at
        // a time, and remembering the cells that had the lowest ratios.  When we
        // proceed from one level to the next, we consider all the children of those
        // cells and keep the best ones.
        //
        // The reason we can restrict the search to children of cells at the
        // previous level is that the ratio above is essentially a function of the
        // local distortions created by projecting the S2 cube space onto the
        // sphere.  These distortions change smoothly over the sphere, so by keeping
        // a fairly large number of candidates ("num_to_keep"), we are essentially
        // keeping all the neighbors of the optimal cell as well.
        var best_score = 1e10
        val best_cells = mutableSetOf<S2CellId>()
        for (level in 0..S2CellId.kMaxLevel) {
            val score = getS2CellIdMinVertexSeparation(level, best_cells).radians
            best_score = min(best_score, score)
        }
        println("min_vertex_sep / snap_radius ratio: %.15f\n".format(best_score))
    }


    fun testS2CellIdSnapFunctionMinEdgeVertexSeparationForLevel() {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for any
        // snap radius at each level.
        val score = getS2CellIdMinEdgeSeparation("min_sep_for_level", object : S2CellIdMinEdgeSeparationFunction {
            override fun apply(level: Int, edge_sep: S1Angle, min_snap_radius: S1Angle, max_snap_radius: S1Angle) = edge_sep.radians / S2CellMetrics.kMinDiag.getValue(level)
        })
        println("min_edge_vertex_sep / kMinDiag ratio: %.15f".format(score))
    }

    fun testS2CellIdSnapFunctionMinEdgeVertexSeparationAtMinSnapRadius() {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for the
        // special case where the minimum snap radius is being used.
        val score = getS2CellIdMinEdgeSeparation("min_sep_at_min_radius", object : S2CellIdMinEdgeSeparationFunction {
            override fun apply(level: Int, edge_sep: S1Angle, min_snap_radius: S1Angle, max_snap_radius: S1Angle): Double {
                val min_radius_at_level = S2CellMetrics.kMaxDiag.getValue(level) / 2
                return if (min_snap_radius.radians <= (1 + 1e-10) * min_radius_at_level)
                    (edge_sep.radians / S2CellMetrics.kMinDiag.getValue(level))
                else 100.0
            }
        })
        println("min_edge_vertex_sep / kMinDiag at MinSnapRadiusForLevel: %.15f".format(score))
    }

    fun testS2CellIdSnapFunctionMinEdgeVertexSeparationSnapRadiusRatio() {
        // Computes the minimum edge separation expressed as a fraction of the
        // maximum snap radius that could yield that edge separation.
        val score = getS2CellIdMinEdgeSeparation("min_sep_snap_radius_ratio", object : S2CellIdMinEdgeSeparationFunction {
            override fun apply(level: Int, edge_sep: S1Angle, min_snap_radius: S1Angle, max_snap_radius: S1Angle): Double {
                return edge_sep.radians / max_snap_radius.radians
            }
        })
        println("min_edge_vertex_sep / snap_radius ratio: %.15f".format(score))
    }

    fun testIntLatLngSnapFunctionMinVertexSeparationSnapRadiusRatio() {
        var best_score = 1e10
        val best_configs = mutableSetOf<IntLatLng>()
        var scale = 18L
        for (lat0 in 0L..9L) {
            best_configs.add(IntLatLng(lat0, 0));
        }
        for (exp in 0..10) {
            val score = getLatLngMinVertexSeparation(scale, 10 * scale, best_configs)
            best_score = min(best_score, score);
            scale *= 10
        }
        println("min_vertex_sep / snap_radius ratio: %.15f".format(best_score))
    }

    fun testIntLatLngSnapFunctionMinEdgeVertexSeparationForLevel() {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for any
        // snap radius at each level.
        val score = getLatLngMinEdgeSeparation("min_sep_for_level", object : LatLngMinEdgeSeparationFunction {
            override fun apply(scale: Long, edge_sep: S1Angle, max_snap_radius: S1Angle): Double {
                val e_unit = M_PI / scale
                return edge_sep.radians / e_unit
            }
        })
        println("min_edge_vertex_sep / e_unit ratio: %.15f".format(score))
    }

    fun testIntLatLngSnapFunctionMinEdgeVertexSeparationSnapRadiusRatio() {
        // Computes the minimum edge separation expressed as a fraction of the
        // maximum snap radius that could yield that edge separation.
        val score = getLatLngMinEdgeSeparation ("min_sep_snap_radius_ratio", object : LatLngMinEdgeSeparationFunction {
            override fun apply(scale: Long, edge_sep: S1Angle, max_snap_radius: S1Angle): Double {
                return edge_sep.radians / max_snap_radius.radians
            }
        })
        println("min_edge_vertex_sep / snap_radius ratio: %.15f".format(score))
    }

    companion object {

        val kSearchRootId = S2CellId.fromFace(0);
        val kSearchFocusId = S2CellId.fromFace(0).child(3)

        fun getMaxVertexDistance(p: S2Point, id: S2CellId): S1Angle {
            val cell = S2Cell(id)
            return maxOf(
                    maxOf(S1Angle(p, cell.getVertex(0)), S1Angle(p, cell.getVertex(1))),
                    maxOf(S1Angle(p, cell.getVertex(2)), S1Angle(p, cell.getVertex(3)))
            )
        }

        // Helper function that computes the vertex separation between "id0" and its
        // neighbors.
        fun updateS2CellIdMinVertexSeparation(id0: S2CellId, scores: MutableList<Pair<S1Angle, S2CellId>>) {
            val site0 = id0.toPoint()
            val nbrs = mutableListOf<S2CellId>()
            id0.appendAllNeighbors(id0.level(), nbrs)
            for (id1 in nbrs) {
                val site1 = id1.toPoint()
                val vertexSep = S1Angle(site0, site1)
                val maxSnapRadius = getMaxVertexDistance(site0, id1)
                assertTrue(maxSnapRadius >= S2CellIdSnapFunction.minSnapRadiusForLevel(id0.level()))
                val r = vertexSep / maxSnapRadius;
                scores.add(Pair(r, id0))
            }
        }

        fun getS2CellIdMinVertexSeparation(level: Int, best_cells: MutableSet<S2CellId>): S1Angle {
            // The worst-case separation ratios always occur when the snap_radius is not
            // much larger than the minimum, since this allows the site spacing to be
            // reduced by as large a fraction as possible.
            //
            // For the minimum vertex separation ratio, we choose a site and one of its
            // 8-way neighbors, then look at the ratio of the distance to the center of
            // that neighbor to the distance to the furthest corner of that neighbor
            // (which is the largest possible snap radius for this configuration).
            var scores = mutableListOf<Pair<S1Angle, S2CellId>>()
            if (level == 0) {
                updateS2CellIdMinVertexSeparation(kSearchRootId, scores)
            } else {
                for (parent in best_cells) {
                    var id0 = parent.childBegin()
                    while (id0 != parent.childEnd()) {
                        updateS2CellIdMinVertexSeparation(id0, scores)
                        id0 = id0.next()
                    }
                }
            }
            // Now sort the entries, print out the "num_to_print" best ones, and keep
            // the best "num_to_keep" of them to seed the next round.
            scores = scores.distinct().toMutableList()
            scores.sortWith { e1, e2 -> ComparisonChain.start().compare(e1.first, e2.first).compare(e1.second, e2.second).result() }
            best_cells.clear()
            var num_to_keep = 300
            var num_to_print = 1
            for (entry in scores) {
                val id = entry.second
                if (--num_to_print >= 0) {
                    val uv = id.getCenterUV()
                    println("Level %2d: min_vertex_sep_ratio = %.15f u=%.6f v=%.6f %s".format(level, entry.first.radians, uv[0], uv[1], id.toToken()))
                }
                if (kSearchFocusId.contains(id) || id.contains(kSearchFocusId)) {
                    if (best_cells.add(id) && --num_to_keep <= 0) break;
                }
            }
            return scores[0].first
        }

        fun getCircumRadius(a: S2Point, b: S2Point, c: S2Point): S1Angle {
            // We return this value is the circumradius is very large.
            val kTooBig = S1Angle.radians(M_PI)
            val turnAngle = S2Measures.turnAngle(a, b, c)
            if (abs(turnAngle.IEEErem(M_PI)) < 1e-2) return kTooBig

            val a2 = (b - c).norm2().toBigDecimal(MathContext.UNLIMITED)
            a2.setScale(2, RoundingMode.HALF_EVEN)
            val b2 = (c - a).norm2().toBigDecimal(MathContext.UNLIMITED)
            b2.setScale(2, RoundingMode.HALF_EVEN)
            val c2 = (a - b).norm2().toBigDecimal(MathContext.UNLIMITED)
            c2.setScale(2, RoundingMode.HALF_EVEN)
            val two = 2.toBigDecimal(MathContext.UNLIMITED)
            two.setScale(2, RoundingMode.HALF_EVEN)
            if (a2 > two || b2 > two || c2 > two) return kTooBig
            val ma = a2 * (b2 + c2 - a2)
            val mb = b2 * (c2 + a2 - b2)
            val mc = c2 * (a2 + b2 - c2)
            val sum = ma + mb + mc

            if (sum.compareTo(BigDecimal.ZERO) == 0) {
                return S1Angle.zero
            }

            val p = (a * (ma / sum).toDouble() + b * (mb / sum).toDouble() + c * (mc / sum).toDouble())
            return S1Angle(p, a)
        }

        fun getNeighbors(id: S2CellId): List<S2CellId> {
            val kNumLayers = 2;
            var nbrs = mutableListOf<S2CellId>()
            nbrs.add(id)
            for (layer in 0 until kNumLayers) {
                val new_nbrs = mutableListOf<S2CellId>()
                for (nbr in nbrs) {
                    nbr.appendAllNeighbors(id.level(), new_nbrs)
                }
                nbrs.addAll(new_nbrs)
                nbrs = nbrs.distinct().toMutableList()
                nbrs.sort()
            }
            return nbrs;
        }

        // S2CellIdMinEdgeSeparationFunction defines an objective function that will
        // be optimized by GetS2CellIdMinEdgeSeparation() by finding worst-case
        // configurations of S2CellIds.  We use this to find the worst cases under
        // various conditions (e.g., when the minimum snap radius at a given level is
        // being used).  The objective function is called for a specific configuration
        // of vertices that are snapped at the given S2CellId level.  "edge_sep" is
        // the edge-vertex distance that is achieved by this configuration, and
        // "min_snap_radius" and "max_snap_radius" are the minimum and maximum snap
        // radii for which this configuration is valid (i.e., where the desired
        // snapping will take place).
        @FunctionalInterface
        interface S2CellIdMinEdgeSeparationFunction {
            fun apply(level: Int, edge_sep: S1Angle, min_snap_radius: S1Angle, max_snap_radius: S1Angle): Double
        }

        // Returns the minimum value of the given objective function over sets of
        // nearby vertices that are designed to minimize the edge-vertex separation
        // when an edge is snapped.
        fun getS2CellIdMinEdgeSeparation(label: String, objective: S2CellIdMinEdgeSeparationFunction, level: Int, best_cells: MutableSet<S2CellId>): Double {
            // To find minimum edge separations, we choose a cell ("id0") and two nearby
            // cells ("id1" and "id2"), where "nearby" is defined by GetNeighbors().
            // Let "site0", "site1", and "site2" be the centers of these cells.  The
            // idea is to consider an input edge E that intersects the Voronoi regions
            // of "site1" and "site2" (and therefore snaps to an edge E' between these
            // sites) but does not not intersect the Voronoi region of "site0" (and
            // therefore can't be snapped to site0).  The goal is to search for snapped
            // edges E' that approach site0 as closely as possible.
            //
            // To do this, we first compute the circumradius of the three cell centers
            // ("site0", "site1", and "site2"); this is the minimum snap radius in order
            // for it to be possible to construct an edge E that snaps to "site1" and
            // "site2" but not to "site0".  We also compute the distance from "site0" to
            // the snapped edge.  Next we find the corner vertex of "id1" and "id2" that
            // is furthest from "site0"; the smaller of these two distances is the
            // maximum snap radius such that "site1" and "site2" can be chosen as
            // sites after choosing "site0".  If the maximum is less than the minimum,
            // then this configuration is rejected; otherwise we evaluate the given
            // objective function and keep the configurations that result in the
            // smallest values.
            //
            // The optimization process works by keeping track of the set of S2CellIds
            // that yielded the best results at the previous level, and exploring all
            // the nearby neighbor combinations of the children of those cells at the
            // next level.  In order to get better coverage, we keep track of the best
            // score and configuration (i.e. the two neighboring cells "id1" and "id2")
            // for each initial cell "id0".
            val best_scores = mutableMapOf<S2CellId, Double>()
            val best_configs = mutableMapOf<S2CellId, Pair<S2CellId, S2CellId>>()
            for (parent in best_cells) {
                var id0 = parent.childBegin(level)
                while (id0 != parent.childEnd(level)) {
                    val site0 = id0.toPoint()
                    val nbrs = getNeighbors(id0)
                    for (id1 in nbrs) {
                        val site1 = id1.toPoint()
                        val max_v1 = getMaxVertexDistance(site0, id1)
                        for (id2 in nbrs) {
                            if (id2 <= id1) continue;
                            val site2 = id2.toPoint()
                            val min_snap_radius = getCircumRadius(site0, site1, site2)
                            if (min_snap_radius > SnapFunction.kMaxSnapRadius) {
                                continue
                            }
                            // Note that it is only the original points *before* snapping that
                            // need to be at least "snap_radius" away from "site0".  The points
                            // after snapping ("site1" and "site2") may be closer.
                            val max_v2 = getMaxVertexDistance(site0, id2)
                            val max_snap_radius = minOf(max_v1, max_v2)
                            if (min_snap_radius > max_snap_radius) continue
                            val minSnapRadiusForLevel = S2CellIdSnapFunction.minSnapRadiusForLevel(level)
                            assertTrue("max_snap_radius = $max_snap_radius <= $minSnapRadiusForLevel", max_snap_radius >= minSnapRadiusForLevel)

                            // This is a valid configuration, so evaluate it.
                            val edge_sep = S2EdgeDistances.getDistance(site0, site1, site2)
                            val score = objective.apply(level, edge_sep, min_snap_radius, max_snap_radius)
                            var best_score = best_scores.getValue(id0)
                            if (best_score == 0.0 || best_score > score) {
                                best_scores[id0] = score
                                best_configs[id0] = Pair(id1, id2)
                            }
                        }
                    }
                    id0 = id0.next()
                }
            }
            // Now sort the entries, print out the "num_to_print" best ones, and
            // generate a set of candidates for the next round by generating all the
            // 8-way neighbors of the best candidates, and keeping up to"num_to_keep" of
            // them.  The results vary slightly according to how many candidates we
            // keep, but the variations are much smaller than the conservative
            // assumptions made by the S2CellIdSnapFunction implementation.
            var num_to_keep = 100
            var num_to_print = 3
            val sorted = mutableListOf<Pair<Double, S2CellId>>()
            for (entry in best_scores) {
                sorted.add(Pair(entry.value, entry.key))
            }
            sorted.sortWith { e1, e2 -> ComparisonChain.start().compare(e1.first, e2.first).compare(e1.second, e2.second).result() }
            best_cells.clear()
            println("Level %d:".format(level))
            for (entry in sorted) {
                val id = entry.second
                if (--num_to_print >= 0) {
                    val uv = id.getCenterUV()
                    val nbrs = best_configs.getValue(id)
                    println("  %s = %.15f u=%7.4f v=%7.4f %s %s %s".format(
                            label, entry.first, uv[0], uv[1], id.toToken(),
                            nbrs.first.toToken(),
                            nbrs.second.toToken())
                    )
                }
                val nbrs = mutableListOf(id)
                id.appendAllNeighbors(id.level(), nbrs)
                for (nbr in nbrs) {
                    // The S2Cell hierarchy has many regions that are symmetrical.  We can
                    // eliminate most of the "duplicates" by restricting the search to cells
                    // in kS2CellIdFocus.
                    if (kSearchFocusId.contains(nbr) || nbr.contains(kSearchFocusId)) {
                        if (best_cells.add(nbr) && --num_to_keep <= 0) {
                            return sorted[0].first
                        }
                    }
                }
            }
            return sorted[0].first
        }

        fun getS2CellIdMinEdgeSeparation(label: String, objective: S2CellIdMinEdgeSeparationFunction): Double {
            var best_score = 1e10;
            val best_cells = mutableSetOf<S2CellId>()
            best_cells.add(kSearchRootId)
            for (level in 0..S2CellId.kMaxLevel) {
                val score = getS2CellIdMinEdgeSeparation(label, objective, level, best_cells)
                best_score = min(best_score, score)
            }
            return best_score
        }


        fun isValid(ll: IntLatLng, scale: Long): Boolean {
            // A coordinate value of "scale" corresponds to 180 degrees.
            return (abs(ll[0]) <= scale / 2 && abs(ll[1]) <= scale)
        }

        fun hasValidVertices(ll: IntLatLng, scale: Long): Boolean {
            // Like IsValid, but excludes latitudes of 90 and longitudes of 180.
            // A coordinate value of "scale" corresponds to 180 degrees.
            return (abs(ll[0]) < scale / 2 && abs(ll[1]) < scale);
        }

        fun rescale(ll: IntLatLng, scale_factor: Double): IntLatLng {
            return IntLatLng(
                    (scale_factor * ll[0]).roundToLong(),
                    (scale_factor * ll[1]).roundToLong()
            )
        }

        fun toPoint(ll: IntLatLng, scale: Long): S2Point {
            return S2LatLng.fromRadians(
                    ll[0] * (M_PI / scale),
                    ll[1] * (M_PI / scale)
            ).toPoint()
        }

        fun getVertex(ll: IntLatLng, scale: Long, i: Int): S2Point {
            // Return the points in CCW order starting from the lower left.
            val dlat = if (i == 0 || i == 3) -1 else 1
            val dlng = if (i == 0 || i == 1) -1 else 1
            return toPoint(IntLatLng(2 * ll.first + dlat, 2 * ll.second + dlng), 2 * scale)
        }

        fun getMaxVertexDistance(p: S2Point, ll: IntLatLng, scale: Long): S1Angle {
            return maxOf(
                    maxOf(
                            S1Angle(p, getVertex(ll, scale, 0)),
                            S1Angle(p, getVertex(ll, scale, 1))
                    ),
                    maxOf(
                            S1Angle(p, getVertex(ll, scale, 2)),
                            S1Angle(p, getVertex(ll, scale, 3))
                    )
            )
        }

        fun getLatLngMinVertexSeparation(old_scale: Long, scale: Long, best_configs: MutableSet<IntLatLng>): Double {
            // The worst-case separation ratios always occur when the snap_radius is not
            // much larger than the minimum, since this allows the site spacing to be
            // reduced by as large a fraction as possible.
            //
            // For the minimum vertex separation ratio, we choose a site and one of its
            // 8-way neighbors, then look at the ratio of the distance to the center of
            // that neighbor to the distance to the furthest corner of that neighbor
            // (which is the largest possible snap radius for this configuration).
            val min_snap_radius_at_scale = S1Angle.radians(M_SQRT1_2 * M_PI / scale)
            var scores = mutableListOf<Pair<Double, IntLatLng>>()
            val scale_factor = scale.toDouble() / old_scale
            for (parent in best_configs) {
                val new_parent: IntLatLng = rescale(parent, scale_factor)
                for (dlat0 in -7L..7L) {
                    val ll0 = new_parent + IntLatLng(dlat0, 0)
                    if (!isValid(ll0, scale) || ll0[0] < 0) continue
                    val site0 = toPoint(ll0, scale)
                    for (dlat1 in 0L..2L) {
                        for (dlng1 in 0L..5L) {
                            val ll1 = ll0 + IntLatLng(dlat1, dlng1)
                            if (ll1 == ll0 || !hasValidVertices(ll1, scale)) continue
                            val max_snap_radius = getMaxVertexDistance(site0, ll1, scale)
                            if (max_snap_radius < min_snap_radius_at_scale) continue
                            val site1 = toPoint(ll1, scale)
                            val vertex_sep = S1Angle(site0, site1)
                            val r = vertex_sep / max_snap_radius;
                            scores.add(Pair(r.radians, ll0))
                        }
                    }
                }
            }
            // Now sort the entries, print out the "num_to_print" best ones, and keep
            // the best "num_to_keep" of them to seed the next round.
            scores = scores.distinct().toMutableList()
            scores.sortWith { e1, e2 -> ComparisonChain.start()
                    .compare(e1.first, e2.first)
                    .compare(e1.second.first, e2.second.first)
                    .compare(e1.second.second, e2.second.second)
                    .result()
            }
            best_configs.clear()
            var num_to_keep = 100
            var num_to_print = 1
            for (entry in scores) {
                if (--num_to_print >= 0) {
                    println("Scale %14d: min_vertex_sep_ratio = %.15f, %s".format(
                            scale, entry.first,
                            toPoint(entry.second, scale).toString())
                    )
                }
                if (best_configs.add(entry.second) && --num_to_keep <= 0) break;
            }
            return scores[0].first
        }

        // A triple of scaled S2LatLng coordinates.  The coordinates are multiplied by
        // (M_PI / scale) to convert them to radians.
        data class LatLngConfig(var scale: Long, var ll0: IntLatLng, var ll1: IntLatLng, var ll2: IntLatLng) : Comparable<LatLngConfig> {

            override fun compareTo(other: LatLngConfig): Int {
                Assertions.assertEQ(scale, other.scale)
                return ComparisonChain.start()
                        .compare(ll0.first, other.ll0.first)
                        .compare(ll0.second, other.ll0.second)
                        .compare(ll1.first, other.ll1.first)
                        .compare(ll1.second, other.ll1.second)
                        .compare(ll2.first, other.ll2.first)
                        .compare(ll2.second, other.ll2.second)
                        .result()
            }

        }

        interface LatLngMinEdgeSeparationFunction {
            fun apply(scale: Long, edge_sep: S1Angle, max_snap_radius: S1Angle): Double
        }


        fun getLatLngMinEdgeSeparation(label: String, objective: LatLngMinEdgeSeparationFunction, scale: Long, best_configs: MutableList<LatLngConfig>): Double {
            val min_snap_radius_at_scale = S1Angle.radians(M_SQRT1_2 * M_PI / scale)
            var scores = mutableListOf<Pair<Double, LatLngConfig>>()
            for (parent in best_configs) {
                // To reduce duplicates, we require that site0 always has longitude 0.
                Assertions.assertEQ(0L, parent.ll0[1])
                val scale_factor = scale.toDouble() / parent.scale
                parent.ll0 = rescale(parent.ll0, scale_factor);
                parent.ll1 = rescale(parent.ll1, scale_factor);
                parent.ll2 = rescale(parent.ll2, scale_factor);
                for (dlat0 in -1L..1L) {
                    val ll0 = parent.ll0 + IntLatLng(dlat0, 0)
                    // To reduce duplicates, we require that site0.latitude >= 0.
                    if (!isValid(ll0, scale) || ll0[0] < 0) continue
                    val site0 = toPoint(ll0, scale)
                    for (dlat1 in -1L..1L) {
                        for (dlng1 in -2L..2L) {
                            val ll1 = parent.ll1 + IntLatLng(dlat0 + dlat1, dlng1)
                            if (ll1 == ll0 || !hasValidVertices(ll1, scale)) continue
                            // Only consider neighbors within 2 latitude units of site0.
                            if (abs(ll1[0] - ll0[0]) > 2) continue

                            val site1 = toPoint(ll1, scale);
                            val max_v1 = getMaxVertexDistance(site0, ll1, scale)
                            for (dlat2 in -1L..1L) {
                                for (dlng2 in -2L..2L) {
                                    val ll2 = parent.ll2 + IntLatLng(dlat0 + dlat2, dlng2)
                                    if (!hasValidVertices(ll2, scale)) continue
                                    // Only consider neighbors within 2 latitude units of site0.
                                    if (abs(ll2[0] - ll0[0]) > 2) continue
                                    // To reduce duplicates, we require ll1 < ll2 lexicographically
                                    // and site2.longitude >= 0.  (It's *not* okay to
                                    // require site1.longitude >= 0, because then some configurations
                                    // with site1.latitude == site2.latitude would be missed.)
                                    if (ll2[0] < ll1[0] || (ll2[0] == ll1[0] && ll2[1] <= ll1[1]) || ll2[1] < 0) continue

                                    val site2 = toPoint(ll2, scale)
                                    val min_snap_radius = getCircumRadius(site0, site1, site2)
                                    if (min_snap_radius > SnapFunction.kMaxSnapRadius) {
                                        continue
                                    }
                                    // Only the original points *before* snapping that need to be at
                                    // least "snap_radius" away from "site0".  The points after
                                    // snapping ("site1" and "site2") may be closer.
                                    val max_v2 = getMaxVertexDistance(site0, ll2, scale)
                                    val max_snap_radius = minOf(max_v1, max_v2)
                                    if (min_snap_radius > max_snap_radius) continue
                                    if (max_snap_radius < min_snap_radius_at_scale) continue

                                    // This is a valid configuration, so evaluate it.
                                    val edge_sep = S2EdgeDistances.getDistance(site0, site1, site2)
                                    val score = objective.apply(scale, edge_sep, max_snap_radius)
                                    val config = LatLngConfig(scale, ll0, ll1, ll2)
                                    scores.add(Pair(score, config))
                                }
                            }
                        }
                    }
                }
            }
            // Now sort the entries, print out the "num_to_print" best ones, and keep
            // the best "num_to_keep" of them to seed the next round.
            scores = scores.distinct().toMutableList()
            scores.sortWith { e1, e2 -> ComparisonChain.start().compare(e1.first, e2.first).compare(e1.second, e2.second).result() }
            best_configs.clear()
            var num_to_keep = 100
            var num_to_print = 3
            println("Scale %d:".format(scale))
            for (entry in scores) {
                val config = entry.second
                val scale = config.scale
                if (--num_to_print >= 0) {
                    println("  %s = %.15f %s %s %s".format(
                            label, entry.first,
                            toPoint(config.ll0, scale).toString(),
                            toPoint(config.ll1, scale).toString(),
                            toPoint(config.ll2, scale).toString())
                    )
                }
                // Optional: filter the candidates to concentrate on a specific region
                // (e.g., the north pole).
                best_configs.add(config)
                if (--num_to_keep <= 0) break;
            }
            return scores[0].first;
        }

        fun getLatLngMinEdgeSeparation(label: String, objective: LatLngMinEdgeSeparationFunction): Double {
            var best_score = 1e10;
            val best_configs = mutableListOf<LatLngConfig>()
            var scale = 6L;  // Initially points are 30 degrees apart.
            val max_lng = scale
            val max_lat = scale / 2
            for (lat0 in 0L..max_lat) {
                for (lat1 in (lat0 - 2)..min(max_lat, lat0 + 2)) {
                    for (lng1 in 0L..max_lng) {
                        for (lat2 in lat1..min(max_lat, lat0 + 2)) {
                            for (lng2 in 0L..max_lng) {
                            val ll0 = IntLatLng(lat0, 0)
                            val ll1 = IntLatLng(lat1, lng1)
                            val ll2 = IntLatLng(lat2, lng2)
                            if (ll2[0] < ll1[0] || (ll2[0] == ll1[0] && ll2[1] <= ll1[1])) continue
                            best_configs.add(LatLngConfig(scale, ll0, ll1, ll2))
                        }
                        }
                    }
                }
            }
            println("Starting with ${best_configs.size} configurations")
            var target_scale = 180L;
            for (exp in 0..10) {
                while (scale < target_scale) {
                    scale = min((1.8 * scale).toLong(), target_scale)
                    val score = getLatLngMinEdgeSeparation (label, objective, scale, best_configs)
                    if (scale == target_scale) {
                        best_score = min(best_score, score);
                    }
                }
                target_scale *= 10
            }
            return best_score
        }

    }

}
