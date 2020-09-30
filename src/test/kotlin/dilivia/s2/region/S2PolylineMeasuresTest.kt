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
package dilivia.s2.region

import com.google.common.geometry.S2.M_PI
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.times
import java.lang.Math.pow
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.sin

class S2PolylineMeasuresTest : S2GeometryTestCase() {

    fun testGetLengthAndCentroidGreatCircles() {
        // Construct random great circles and divide them randomly into segments.
        // Then make sure that the length and centroid are correct.  Note that
        // because of the way the centroid is computed, it does not matter how
        // we split the great circle into segments.

        repeat(100) {
          // Choose a coordinate frame for the great circle.
          val (x, y, z) = S2Random.randomFrame()

            val line = mutableListOf<S2Point>()
            var theta = 0.0
            while (theta < 2 * M_PI) {
                line.add(cos(theta) * x + sin(theta) * y);
                theta += pow(S2Random.randomDouble(), 10.0)
            }
            // Close the circle.
            line.add(line[0]);
            val length = S2Polyline.getLength(line);
            assertTrue(abs(length.radians - 2 * M_PI) <= 2e-14);
            val centroid = S2Polyline.getCentroid(line);
            assertTrue(centroid.norm() <= 2e-14);
        }
    }

}
