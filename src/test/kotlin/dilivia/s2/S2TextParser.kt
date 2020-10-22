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

import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.shape.S2LaxPolygonShape
import dilivia.s2.shape.S2LaxPolylineShape
import dilivia.s2.shape.S2PointVectorShape
import mu.KotlinLogging

object S2TextParser {

    private val logger = KotlinLogging.logger {  }

    // Parses a string of one or more latitude-longitude coordinates in degrees,
    // and return the corresponding vector of S2LatLng points.
    // Examples of the input format:
    //     ""                            // no points
    //     "-20:150"                     // one point
    //     "-20:150, -20:151, -19:150"   // three points
    fun parseLatLngs(str: String): List<S2LatLng> {
        val latlngs = mutableListOf<S2LatLng>()
        check(parseLatLngs(str, latlngs)) { ": str == \"$str\"" }
        return latlngs
    }

    // As above, but does not S2_CHECK-fail on invalid input. Returns true if
    // conversion is successful.
    fun parseLatLngs(str: String, latlngs: MutableList<S2LatLng>): Boolean {
        val ps = str.split(',').map { pointStr -> pointStr.trim() }.map { pointStr ->
            val coordsStr = pointStr.split(':')
            if (coordsStr.size != 2) return false
            Pair(coordsStr[0].trim(), coordsStr[1].trim())
        }

        for (p in ps) {
            val lat = p.first.toDoubleOrNull() ?: return false
            val lng = p.second.toDoubleOrNull() ?: return false
            latlngs.add(S2LatLng.fromDegrees(lat, lng))
        }
        return true
    }

    // Parses a string in the same format as ParseLatLngs, and return the
    // corresponding vector of S2Point values.
    fun parsePoints(str: String): List<S2Point> {
        val vertices = mutableListOf<S2Point>()
        check(parsePoints(str, vertices)) { ": str == \"$str\"" }
        return vertices
    }

    // As above, but does not S2_CHECK-fail on invalid input. Returns true if
    // conversion is successful.
    fun parsePoints(str: String, vertices: MutableList<S2Point>): Boolean {
        val latlngs = mutableListOf<S2LatLng>()
        if (!parseLatLngs(str, latlngs)) return false
        for (latlng in latlngs) {
            vertices.add(latlng.toPoint())
        }
        return true
    }

    fun makePoint(str: String): S2Point {
        val vertices = mutableListOf<S2Point>()
        check (parsePoints(str, vertices) && vertices.size == 1) { ": str == \"$str\"" }
        return vertices[0];
    }

    fun makeLoop(str: String, check: Boolean = true): S2Loop {
        return when (str) {
            "empty" -> S2Loop(S2Loop.kEmpty)
            "full" -> S2Loop(S2Loop.kFull)
            else -> {
                val vertices = mutableListOf<S2Point>()
                check (parsePoints(str, vertices)) { ": str == \"$str\"" }
                S2Loop(vertices, check = check);
            }
        }
    }

    // Like MakePolyline, but returns an S2LaxPolylineShape instead.
    fun makeLaxPolyline(str: String): S2LaxPolylineShape {
        val vertices = mutableListOf<S2Point>()
        check (parsePoints(str, vertices))  { ": str == \"$str\"" }
        return S2LaxPolylineShape(vertices)
    }

    // Parses an S2CellId in the format "f/dd..d" where "f" is a digit in the
    // range [0-5] representing the S2CellId face, and "dd..d" is a string of
    // digits in the range [0-3] representing each child's position with respect
    // to its parent.  (Note that the latter string may be empty.)
    //
    // For example "4/" represents S2CellId::FromFace(4), and "3/02" represents
    // S2CellId::FromFace(3).child(0).child(2).
    //
    // This function is a wrapper for S2CellId::FromDebugString().
    fun makeCellId(str: String): S2CellId {
        val cellId = S2CellId.fromDebugString(str)
        check(cellId != S2CellId.none()) { "Invalid cell id: str == \"$str\"" }
        return cellId;
    }

    // Parses a comma-separated list of S2CellIds in the format above, and returns
    // the corresponding S2CellUnion.  (Note that S2CellUnions are automatically
    // normalized by sorting, removing duplicates, and replacing groups of 4 child
    // cells by their parent cell.)
    fun makeCellUnion(str: String): S2CellUnion {
        val cellIds = mutableListOf<S2CellId>()
        str.split(",").map { it.trim() }.forEach { cell_str -> cellIds.add(makeCellId(cell_str)) }
        return S2CellUnion(cellIds)
    }


    // Returns a MutableS2ShapeIndex containing the points, polylines, and loops
    // (in the form of a single polygon) described by the following format:
    //
    //   point1|point2|... # line1|line2|... # polygon1|polygon2|...
    //
    // Examples:
    //   1:2 | 2:3 # #                     // Two points
    //   # 0:0, 1:1, 2:2 | 3:3, 4:4 #      // Two polylines
    //   # # 0:0, 0:3, 3:0; 1:1, 2:1, 1:2  // Two nested loops (one polygon)
    //   5:5 # 6:6, 7:7 # 0:0, 0:1, 1:0    // One of each
    //   # # empty                         // One empty polygon
    //   # # empty | full                  // One empty polygon, one full polygon
    //
    // Loops should be directed so that the region's interior is on the left.
    // Loops can be degenerate (they do not need to meet S2Loop requirements).
    //
    // CAVEAT: Because whitespace is ignored, empty polygons must be specified
    //         as the string "empty" rather than as the empty string ("").
    fun makeIndex(str: String): MutableS2ShapeIndex {
        val index = MutableS2ShapeIndex()
        check(makeIndex(str, index)) { ": str == \"$str\"" }
        return index;
    }

    // As above, but does not S2_CHECK-fail on invalid input. Returns true if
    // conversion is successful.
    fun makeIndex(str: String, index: MutableS2ShapeIndex): Boolean {
        val strs = str.split('#').map { it.trim() }
        check(3 == strs.size) { "Must contain two # characters: $str" }

        val points = mutableListOf<S2Point>()
        for (point_str in strs[0].split('|').map { it.trim() }.filter { it.isNotEmpty() }) {
            val point = makePoint(point_str)
            points.add(point)
        }
        if (points.isNotEmpty()) {
            index.add(S2PointVectorShape(points = points))
        }

        for (line_str in strs[1].split('|').map { it.trim() }.filter { it.isNotEmpty() }) {
            val laxPolyline = makeLaxPolyline(line_str)
            index.add(laxPolyline)
        }

        for (polygon_str in strs[2].split('|').map { it.trim() }.filter { it.isNotEmpty() }) {
            val lax_polygon = makeLaxPolygon(polygon_str)
            index.add(lax_polygon)
        }
        return true
    }

    fun makePolygon(str: String, debug: S2Debug = S2Debug.ALLOW): S2Polygon = internalMakePolygon(str, debug, true)

   fun makeLaxPolygon(str: String): S2LaxPolygonShape {
       val loopStrs = str.trim().let {
           var s = it
           while(s.endsWith(";")) s = s.removeSuffix(";").trim()
           s
       }.split(';')
       val loops = mutableListOf<MutableList<S2Point>>()
       for (loop_str in loopStrs) {
           if (loop_str == "full") {
               loops.add(mutableListOf<S2Point>())
           } else if (loop_str != "empty") {
               val points = mutableListOf<S2Point>()
               check (parsePoints(loop_str, points)) { "Fail to parse loop: $loop_str" }
               loops.add(points)
           }
       }
       return S2LaxPolygonShape(loops)
    }


    fun makeVerbatimPolygon(str: String): S2Polygon {
        logger.debug { "Create verbatim Polygon: $str" }
        return internalMakePolygon(str, S2Debug.ALLOW, false)
    }

    private fun internalMakePolygon(str: String, debug: S2Debug, normalize_loops: Boolean): S2Polygon {
        val loopStrs = str.trim().let {
            var s = it
            while(s.endsWith(";")) s = s.removeSuffix(";").trim()
            s
        }.split(';')
        val loops = mutableListOf<S2Loop>()
        for (loop_str in loopStrs) {
            logger.trace { "internalMakePolygon: create loop $loop_str" }
            var s = loop_str.trim()
            if (s.isBlank()) s = "empty"
            val loop = makeLoop(s, debug == S2Debug.ALLOW)
            // Don't normalize loops that were explicitly specified as "full".
            if (normalize_loops && !loop.isFull()) loop.normalize()
            loops.add(loop)
        }
        return S2Polygon(loops, debug == S2Debug.ALLOW)
    }

    fun toString(polygon: S2Polygon, loop_separator: String = " ; "): String {
        if (polygon.isEmpty()) {
            return "empty"
        } else if (polygon.isFull()) {
            return "full"
        }
        var out = ""
        for (i in 0 until polygon.numLoops()) {
            if (i > 0) out += loop_separator
            val loop = polygon.loop(i)
            out += loop.vertices().joinToString(", ")
        }
        return polygon.toDebugString(loop_separator)
    }
}
