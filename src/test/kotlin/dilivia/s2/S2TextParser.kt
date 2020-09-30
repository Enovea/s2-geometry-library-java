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

import dilivia.s2.region.S2Loop
import dilivia.s2.shape.S2LaxPolylineShape

object S2TextParser {

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

}