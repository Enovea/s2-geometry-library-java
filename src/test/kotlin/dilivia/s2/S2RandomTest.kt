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

import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2LatLngRect
import junit.framework.TestCase

class S2RandomTest: TestCase() {

    fun testSamplePointInCap() {
        var cap = S2Cap.fromCenterAngle(S2Point(1.0, 0.0, 0.0), S1Angle.degrees(4.0))
        repeat(10000) {
            val point = S2Random.samplePoint(cap)
            assertTrue(cap.contains(point))
        }

        repeat(100) {
            val center = S2Random.randomPoint()
            val angle = S2Random.randomDouble(1.0, 10.0)
            cap = S2Cap.fromCenterAngle(center, S1Angle.degrees(angle))
            repeat(10000) {
                val point = S2Random.samplePoint(cap)
                assertTrue(cap.contains(point))
            }
        }
    }

    fun testSampleCapInRect() {
        repeat(100) {
            val p1 = S2LatLng.fromPoint(S2Random.randomPoint())
            val p2 = S2LatLng.fromPoint(S2Random.randomPoint())
            val rect = S2LatLngRect.fromPointPair(p1, p2)
            repeat(10000) {
                val point = S2Random.samplePoint(rect)
                assertTrue(rect.contains(point))
            }
        }
    }
}
