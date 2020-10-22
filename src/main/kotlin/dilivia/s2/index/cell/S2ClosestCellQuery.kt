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
package dilivia.s2.index.cell

import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2EdgeDistances
import dilivia.s2.S2Point
import dilivia.s2.index.S2MinDistance
import dilivia.s2.index.S2MinDistanceCellTarget
import dilivia.s2.index.S2MinDistanceCellUnionTarget
import dilivia.s2.index.S2MinDistanceEdgeTarget
import dilivia.s2.index.S2MinDistanceFactory
import dilivia.s2.index.S2MinDistancePointTarget
import dilivia.s2.index.S2MinDistanceShapeIndexTarget
import dilivia.s2.index.S2MinDistanceTarget
import dilivia.s2.index.cell.S2ClosestCellQueryBase.Result
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion


// S2ClosestCellQuery is a helper class for finding the closest cell(s) to a
// given point, edge, S2Cell, S2CellUnion, or geometry collection.  A typical
// use case would be to add a collection of S2Cell coverings to an S2CellIndex
// (representing a collection of original geometry), and then use
// S2ClosestCellQuery to find all coverings that are within a given distance
// of some target geometry (which could be represented exactly, or could also
// be a covering).  The distance to the original geometry corresponding to
// each covering could then be measured more precisely if desired.
//
// For example, here is how to find all cells that are closer than
// "distance_limit" to a given target point:
//
//   S2ClosestCellQuery query(&cell_index);
//   query.mutable_options()->set_max_distance(distance_limit);
//   S2ClosestCellQuery::PointTarget target(target_point);
//   for (const auto& result : query.FindClosestCells(&target)) {
//     // result.distance() is the distance to the target.
//     // result.cell_id() is the indexed S2CellId.
//     // result.label() is the integer label associated with the S2CellId.
//     DoSomething(target_point, result);
//   }
//
// You can find either the k closest cells, or all cells within a given
// radius, or both (i.e., the k closest cells up to a given maximum radius).
// By default *all* cells are returned, so you should always specify either
// max_results() or max_distance() or both.  You can also restrict the results
// to cells that intersect a given S2Region; for example:
//
//   S2LatLngRect rect(...);
//   query.mutable_options()->set_region(&rect);  // Does *not* take ownership.
//
// There is a FindClosestCell() convenience method that returns the closest
// cell.  However, if you only need to test whether the distance is above or
// below a given threshold (e.g., 10 km), it is typically much faster to use
// the IsDistanceLess() method instead.  Unlike FindClosestCell(), this method
// stops as soon as it can prove that the minimum distance is either above or
// below the threshold.  Example usage:
//
//   if (query.IsDistanceLess(&target, limit_distance)) ...
//
// To find the closest cells to a query edge rather than a point, use:
//
//   S2ClosestCellQuery::EdgeTarget target(v0, v1);
//   query.FindClosestCells(&target);
//
// Similarly you can find the closest cells to an S2Cell using an
// S2ClosestCellQuery::CellTarget, you can find the closest cells to an
// S2CellUnion using an S2ClosestCellQuery::CellUnionTarget, and you can find
// the closest cells to an arbitrary collection of points, polylines, and
// polygons by using an S2ClosestCellQuery::ShapeIndexTarget.
//
// The implementation is designed to be fast for both simple and complex
// geometric objects.
class S2ClosestCellQuery {
    // See S2ClosestCellQueryBase for full documentation.

    // Convenience constructor that calls Init().  Options may be specified here
    // or changed at any time using the mutable_options() accessor method.
    //
    // REQUIRES: "index" must persist for the lifetime of this object.
    // REQUIRES: ReInit() must be called if "index" is modified.
    constructor(index: S2CellIndex, options: Options = Options()) {
        init(index, options)
    }

    // Default constructor; requires Init() to be called.
    constructor() {

    }

    lateinit var options: Options
    private val base = S2ClosestCellQueryBase(distanceFactory = S2MinDistanceFactory)


    // Initializes the query.  Options may be specified here or changed at any
    // time using the mutable_options() accessor method.
    //
    // REQUIRES: "index" must persist for the lifetime of this object.
    // REQUIRES: ReInit() must be called if "index" is modified.
    fun init(index: S2CellIndex, options: Options = Options()) {
        this.options = options
        base.init(index)
    }

    // Reinitializes the query.  This method must be called if the underlying
    // S2CellIndex is modified (by calling Clear() and Build() again).
    fun reInit() {
        base.reInit()
    }

    // Returns a reference to the underlying S2CellIndex.
    fun index(): S2CellIndex = base.index()

    // Returns the closest cells to the given target that satisfy the current
    // options.  This method may be called multiple times.
    fun findClosestCells(target: S2MinDistanceTarget): List<Result<S2MinDistance>> = base.findClosestCells(target, options)

    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    fun findClosestCells(target: S2MinDistanceTarget, results: MutableList<Result<S2MinDistance>>) {
        base.findClosestCells(target, options, results)
    }

    //////////////////////// Convenience Methods ////////////////////////

    // Returns the closest cell to the target.  If no cell satisfies the search
    // criteria, then the Result object will have distance == Infinity() and
    // is_empty() == true.
    fun findClosestCell(target: S2MinDistanceTarget): Result<S2MinDistance>? {
        val tmpOptions = options.clone()
        tmpOptions.setMaxResults(1)
        return base.findClosestCell(target, tmpOptions)
    }

    // Returns the minimum distance to the target.  If the index or target is
    // empty, returns S1ChordAngle::Infinity().
    //
    // Use IsDistanceLess() if you only want to compare the distance against a
    // threshold value, since it is often much faster.
    fun getDistance(target: S2MinDistanceTarget): S1ChordAngle? = findClosestCell(target)?.distance?.value

    // Returns true if the distance to "target" is less than "limit".
    //
    // This method is usually much faster than GetDistance(), since it is much
    // less work to determine whether the minimum distance is above or below a
    // threshold than it is to calculate the actual minimum distance.
    fun isDistanceLess(target: S2MinDistanceTarget, limit: S1ChordAngle): Boolean {
        val tmpOptions = options.clone()
        tmpOptions.setMaxResults(1)
        tmpOptions.setMaxDistance(limit)
        tmpOptions.maxError = S1ChordAngle.straight
        return !base.findClosestCell(target, tmpOptions).isEmpty()
    }

    // Like IsDistanceLess(), but also returns true if the distance to "target"
    // is exactly equal to "limit".
    fun isDistanceLessOrEqual(target: S2MinDistanceTarget, limit: S1ChordAngle): Boolean {
        val tmpOptions = options.clone()
        tmpOptions.setMaxResults(1)
        tmpOptions.setInclusiveMaxDistance(limit)
        tmpOptions.maxError = S1ChordAngle.straight
        return !base.findClosestCell(target, tmpOptions).isEmpty()
    }

    // Like IsDistanceLessOrEqual(), except that "limit" is increased by the
    // maximum error in the distance calculation.  This ensures that this
    // function returns true whenever the true, exact distance is less than
    // or equal to "limit".
    fun isConservativeDistanceLessOrEqual(target: S2MinDistanceTarget, limit: S1ChordAngle): Boolean {
        val tmpOptions = options.clone()
        tmpOptions.setMaxResults(1)
        tmpOptions.setConservativeMaxDistance(limit)
        tmpOptions.maxError = S1ChordAngle.straight
        return !base.findClosestCell(target, tmpOptions).isEmpty()
    }

    // Options that control the set of cells returned.  Note that by default
    // *all* cells are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    class Options : S2ClosestCellQueryBase.Options<S2MinDistance>(distanceFactory = S2MinDistanceFactory) {
        // See S2ClosestCellQueryBase::Options for the full set of options.

        // Specifies that only cells whose distance to the target is less than
        // "max_distance" should be returned.
        //
        // Note that cells whose distance is exactly equal to "max_distance" are
        // not returned.  Normally this doesn't matter, because distances are not
        // computed exactly in the first place, but if such cells are needed then
        // see set_inclusive_max_distance() below.
        //
        // DEFAULT: Distance::Infinity()
        fun setMaxDistance(maxDistance: S1ChordAngle) {
            this.maxDistance = S2MinDistance(maxDistance)
        }

        // Like set_max_distance(), except that cells whose distance is exactly
        // equal to "max_distance" are also returned.  Equivalent to calling
        // set_max_distance(max_distance.Successor()).
        fun setInclusiveMaxDistance(maxDistance: S1ChordAngle) {
            setMaxDistance(maxDistance.successor())
        }

        // Like set_inclusive_max_distance(), except that "max_distance" is also
        // increased by the maximum error in the distance calculation.  This
        // ensures that all cells whose true distance is less than or equal to
        // "max_distance" will be returned (along with some cells whose true
        // distance is slightly greater).
        //
        // Algorithms that need to do exact distance comparisons can use this
        // option to find a set of candidate cells that can then be filtered
        // further (e.g., using s2pred::CompareDistance).
        fun setConservativeMaxDistance(maxDistance: S1ChordAngle) {
            this.maxDistance = S2MinDistance(maxDistance.plusError(
                    S2EdgeDistances.getUpdateMinDistanceMaxError(maxDistance)).successor())
        }

        // Versions of set_max_distance that take an S1Angle argument.  (Note that
        // these functions require a conversion, and that the S1ChordAngle versions
        // are preferred.)
        fun setMaxDistance(maxDistance: S1Angle) {
            this.maxDistance = S2MinDistance(maxDistance)
        }

        fun setInclusiveMaxDistance(maxDistance: S1Angle) {
            setInclusiveMaxDistance(S1ChordAngle(maxDistance))
        }

        fun setConservativeMaxDistance(max_distance: S1Angle) {
            setConservativeMaxDistance(S1ChordAngle(max_distance))
        }

        // See S2ClosestCellQueryBase::Options for documentation.
        fun setMaxError(maxError: S1Angle) {
            this.maxError = S1ChordAngle(maxError)
        }

        public override fun clone(): Options {
            val clone = Options()
            clone.setMaxResults(getMaxResults())
            clone.maxDistance = maxDistance
            clone.maxError = maxError
            clone.region = region
            clone.useBruteForce = useBruteForce
            return clone
        }

    }

    // The thresholds for using the brute force algorithm are generally tuned to
    // optimize IsDistanceLess (which compares the distance against a threshold)
    // rather than FindClosest (which actually computes the minimum distance).
    // This is because the former operation is (1) more common, (2) inherently
    // faster, and (3) closely related to finding all cells within a given
    // distance, which is also very common.

    // Target subtype that computes the closest distance to a point.
    class PointTarget(point: S2Point) : S2MinDistancePointTarget(point) {

        override fun maxBruteForceIndexSize(): Int {
            // Break-even points:                   Point cloud      Cap coverings
            // BM_FindClosest                                18                 16
            // BM_IsDistanceLess                              8                  9
            return 9
        }

    }

    // Target subtype that computes the closest distance to an edge.
    class EdgeTarget(a: S2Point, b: S2Point) : S2MinDistanceEdgeTarget(a, b) {

        override fun maxBruteForceIndexSize(): Int {
            // Break-even points:                   Point cloud      Cap coverings
            // BM_FindClosestToLongEdge                      14                 16
            // BM_IsDistanceLessToLongEdge                    5                  5
            return 5
        }

    }

    // Target subtype that computes the closest distance to an S2Cell
    // (including the interior of the cell).
    class CellTarget(cell: S2Cell) : S2MinDistanceCellTarget(cell) {

        override fun maxBruteForceIndexSize(): Int {
            // Break-even points:                   Point cloud      Cap coverings
            // BM_FindClosestToSmallCell                     12                 13
            // BM_IsDistanceLessToSmallCell                   6                  6
            //
            // Note that the primary use of CellTarget is to implement CellUnionTarget,
            // and therefore it is very important to optimize for the case where a
            // distance limit has been specified.
            return 6
        }

    }

    // Target subtype that computes the closest distance to an S2CellUnion.
    class CellUnionTarget(cellUnion: S2CellUnion) : S2MinDistanceCellUnionTarget(cellUnion) {

        override fun maxBruteForceIndexSize(): Int {
            // Break-even points:                   Point cloud      Cap coverings
            // BM_FindClosestToSmallCoarseCellUnion          12                 10
            // BM_IsDistanceLessToSmallCoarseCellUnion        7                  6
            return 8
        }

    }

    // Target subtype that computes the closest distance to an S2ShapeIndex
    // (an arbitrary collection of points, polylines, and/or polygons).
    //
    // By default, distances are measured to the boundary and interior of
    // polygons in the S2ShapeIndex rather than to polygon boundaries only.
    // If you wish to change this behavior, you may call
    //
    //   target.set_include_interiors(false);
    //
    // (see S2MinDistanceShapeIndexTarget for details).
    class ShapeIndexTarget(index: S2ShapeIndex) : S2MinDistanceShapeIndexTarget(index) {

        override fun maxBruteForceIndexSize(): Int {
            // Break-even points:                   Point cloud      Cap coverings
            // BM_FindClosestToSmallCoarseShapeIndex         10                  8
            // BM_IsDistanceLessToSmallCoarseShapeIndex       7                  6
            return 7
        }

    }

}
