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
package dilivia.s2.index.point

import dilivia.s2.math.Matrix3x3
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.region.S2Cap
import kotlin.math.ceil
import kotlin.math.sqrt
import kotlin.math.tan

// Generates vertices on a square grid that includes the entire given cap.
class GridPointIndexFactory : PointIndexFactory {

    override fun addPoints(indexCap: S2Cap, numPoints: Int, index: S2PointIndex<Int>) {
        val sqrt_num_points = ceil(sqrt(numPoints.toDouble())).toInt()
        val frame = S2Random.randomFrameAt(indexCap.center)
        val radius = indexCap.radius().radians
        val spacing = 2 * radius / sqrt_num_points;
        for (i in 0 until sqrt_num_points) {
            for (j in 0 until sqrt_num_points) {
                val point = S2Point(tan((i + 0.5) * spacing - radius), tan((j + 0.5) * spacing - radius), 1.0)
                index.add(S2Point.fromFrame(Matrix3x3.fromCols(frame), point.normalize()), i * sqrt_num_points + j)
            }
        }
    }
}
