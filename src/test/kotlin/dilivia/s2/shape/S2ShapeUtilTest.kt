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
package dilivia.s2.shape

import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2TextParser


class S2ShapeUtilTest : S2GeometryTestCase() {

    fun testContainsBruteForceNoInterior() {
        // Defines a polyline that almost entirely encloses the point 0:0.
        val polyline = S2TextParser.makeLaxPolyline("0:0, 0:1, 1:-1, -1:-1, -1e9:1");
        assertFalse(S2ShapeUtil.containsBruteForce(polyline, makePoint("0:0")))
    }
/*
    fun testContainsBruteForceContainsReferencePoint() {
        // Checks that ContainsBruteForce agrees with GetReferencePoint.
        val polygon = makeLaxPolygon("0:0, 0:1, 1:-1, -1:-1, -1e9:1");
        val ref = polygon.getReferencePoint()
        assertEquals(ref.contained, containsBruteForce(polygon, ref.point))
    }

    fun testContainsBruteForceConsistentWithS2Loop() {
        // Checks that ContainsBruteForce agrees with S2Loop::Contains().
        val loop = S2Loop.makeRegularLoop(makePoint("89:-179"), S1Angle.degrees(10), 100)
        val shape = S2Loop.Shape(0 , loop)
        for (i in 0 until loop.numVertices()) {
            assertEquals(loop.contains(loop.vertex(i)), S2ShapeUtil.containsBruteForce(shape, loop.vertex(i)))
        }
    }

 */
}

