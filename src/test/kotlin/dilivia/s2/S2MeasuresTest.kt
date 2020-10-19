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

import dilivia.s2.S2.M_PI
import dilivia.s2.S2.M_PI_2
import dilivia.s2.S2.M_PI_4
import mu.KotlinLogging
import kotlin.math.abs
import kotlin.math.max

class S2MeasuresTest : S2GeometryTestCase() {

  private val logger = KotlinLogging.logger {  }

  fun testAngleMethods() {
    val pz = S2Point(0, 0, 1)
    val p000 = S2Point(1, 0, 0)
    val p045 = S2Point(1, 1, 0).normalize()
    val p090 = S2Point(0, 1, 0)
    val p180 = S2Point(-1, 0, 0)

    assertEquals(S2Measures.angle(p000, pz, p045), M_PI_4)
    assertEquals(S2Measures.turnAngle(p000, pz, p045), -3 * M_PI_4)

    assertEquals(S2Measures.angle(p045, pz, p180), 3 * M_PI_4)
    assertEquals(S2Measures.turnAngle(p045, pz, p180), -M_PI_4)

    assertEquals(S2Measures.angle(p000, pz, p180), M_PI)
    assertEquals(S2Measures.turnAngle(p000, pz, p180), -0.0)

    assertEquals(S2Measures.angle(pz, p000, p045), M_PI_2)
    assertEquals(S2Measures.turnAngle(pz, p000, p045), M_PI_2)

    assertEquals(S2Measures.angle(pz, p000, pz), 0.0)
    assertEquals(abs(S2Measures.turnAngle(pz, p000, pz)), M_PI)
  }

  fun testAreaMethods() {
    val pz = S2Point(0, 0, 1)
    val p000 = S2Point(1, 0, 0)
    val p045 = S2Point(1, 1, 0).normalize()
    val p090 = S2Point(0, 1, 0)
    val p180 = S2Point(-1, 0, 0)

    assertEquals(S2Measures.area(p000, p090, pz), M_PI_2)
    assertEquals(S2Measures.area(p045, pz, p180), 3 * M_PI_4)

    // Make sure that Area() has good *relative* accuracy even for
    // very small areas.
    val eps = 1e-10
    val pepsx = S2Point(eps, 0.0, 1.0).normalize()
    val pepsy = S2Point(0.0, eps, 1.0).normalize()
    val expected1 = 0.5 * eps * eps
    assertEquals(S2Measures.area(pepsx, pepsy, pz), expected1, 1e-14 * expected1)

    // Make sure that it can handle degenerate triangles.
    val pr = S2Point(0.257, -0.5723, 0.112).normalize()
    val pq = S2Point(-0.747, 0.401, 0.2235).normalize()
    assertEquals(S2Measures.area(pr, pr, pr), 0.0)
    // The following test is not exact due to rounding error.
    assertEquals(S2Measures.area(pr, pq, pr), 0.0, 1e-15)
    assertEquals(S2Measures.area(p000, p045, p090), 0.0)

    var max_girard = 0.0
    repeat(10000) {
      val p0 = S2Random.randomPoint()
      val d1 = S2Random.randomPoint()
      val d2 = S2Random.randomPoint()
      val p1 = (p0 + 1e-15 * d1).normalize()
      val p2 = (p0 + 1e-15 * d2).normalize()
      // The actual displacement can be as much as 1.2e-15 due to roundoff.
      // This yields a maximum triangle area of about 0.7e-30.
      assertTrue(S2Measures.area(p0, p1, p2) <= 0.7e-30)
      max_girard = max(max_girard, S2Measures.girardArea(p0, p1, p2))
    }
    // This check only passes if GirardArea() uses RobustCrossProd().
    logger.info { "Worst case Girard for triangle area 1e-30: $max_girard" }
    assertTrue(max_girard <= 1e-14)

    // Try a very long and skinny triangle.
    val p045eps = S2Point (1.0, 1.0, eps).normalize()
    val expected2 = 5.8578643762690495119753e-11;  // Mathematica.
    assertEquals(S2Measures.area(p000, p045eps, p090), expected2, 1e-9 * expected2)

    // Triangles with near-180 degree edges that sum to a quarter-sphere.
    val eps2 = 1e-14
    val p000eps2 = S2Point (1.0, 0.1*eps2, eps2).normalize()
    val quarter_area1 = S2Measures.area(p000eps2, p000, p045) +
            S2Measures.area(p000eps2, p045, p180) +
            S2Measures.area(p000eps2, p180, pz) +
            S2Measures.area(p000eps2, pz, p000)
    assertEquals(quarter_area1, M_PI, 1e-15)

    // Four other triangles that sum to a quarter-sphere.
    val p045eps2 = S2Point (1.0, 1.0, eps2).normalize()
    val quarter_area2 = S2Measures.area(p045eps2, p000, p045) +
            S2Measures.area(p045eps2, p045, p180) +
            S2Measures.area(p045eps2, p180, pz) +
            S2Measures.area(p045eps2, pz, p000)
    assertEquals(quarter_area2, M_PI)

    // Compute the area of a hemisphere using four triangles with one near-180
    // degree edge and one near-degenerate edge.
    repeat(100) {
      val lng = 2 * M_PI * S2Random.randomDouble()
      val p0 = S2LatLng .fromRadians(1e-20, lng).normalized().toPoint()
      val p1 = S2LatLng .fromRadians(0.0, lng).normalized().toPoint()
      val p2_lng = lng +S2Random.randomDouble()
      val p2 = S2LatLng .fromRadians(0.0, p2_lng).normalized().toPoint()
      val p3 = S2LatLng .fromRadians(0.0, lng + M_PI).normalized().toPoint()
      val p4 = S2LatLng .fromRadians(0.0, lng + 5.0).normalized().toPoint()
      val area =(S2Measures.area(p0, p1, p2) + S2Measures.area(p0, p2, p3) + S2Measures.area(p0, p3, p4) + S2Measures.area(p0, p4, p1))
      assertEquals(area, 2 * M_PI, 2e-15)
    }

    // This tests a case where the triangle has zero area, but S2Measures.area()
    // computes (dmin > 0) due to rounding errors.
    assertEquals(0.0, S2Measures.area(S2LatLng.fromDegrees(-45, -170).toPoint(),
            S2LatLng.fromDegrees(45, -170).toPoint(),
            S2LatLng.fromDegrees(0, -170).toPoint()))
  }

}
