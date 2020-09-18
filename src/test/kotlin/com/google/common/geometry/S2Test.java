/*
 * Copyright 2005 Google Inc.
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
package com.google.common.geometry;

import dilivia.s2.S2GeometryTestCase;
import dilivia.s2.S2Point;
import dilivia.s2.S2Random;

import java.util.logging.Logger;

public strictfp class S2Test extends S2GeometryTestCase {

    private static Logger logger = Logger.getLogger(S2Test.class.getName());

    // Test helper methods for testing the traversal order.
    private static int swapAxes(int ij) {
        return ((ij >> 1) & 1) + ((ij & 1) << 1);
    }

    private static int invertBits(int ij) {
        return ij ^ 3;
    }


    public void testAngleArea() {
        S2Point pz = new S2Point(0, 0, 1);
        S2Point p000 = new S2Point(1, 0, 0);
        S2Point p045 = new S2Point(1, 1, 0);
        S2Point p090 = new S2Point(0, 1, 0);
        S2Point p180 = new S2Point(-1, 0, 0);
        assertDoubleNear(S2.angle(p000, pz, p045), S2.M_PI_4);
        assertDoubleNear(S2.angle(p045, pz, p180), 3 * S2.M_PI_4);
        assertDoubleNear(S2.angle(p000, pz, p180), S2.M_PI);
        assertDoubleNear(S2.angle(pz, p000, pz), 0);
        assertDoubleNear(S2.angle(pz, p000, p045), S2.M_PI_2);

        assertDoubleNear(S2.area(p000, p090, pz), S2.M_PI_2);
        assertDoubleNear(S2.area(p045, pz, p180), 3 * S2.M_PI_4);

        // Make sure that area() has good *relative* accuracy even for
        // very small areas.
        final double eps = 1e-10;
        S2Point pepsx = new S2Point(eps, 0, 1);
        S2Point pepsy = new S2Point(0, eps, 1);
        double expected1 = 0.5 * eps * eps;
        assertDoubleNear(S2.area(pepsx, pepsy, pz), expected1, 1e-14 * expected1);

        // Make sure that it can handle degenerate triangles.
        S2Point pr = new S2Point(0.257, -0.5723, 0.112);
        S2Point pq = new S2Point(-0.747, 0.401, 0.2235);
        assertEquals(S2.area(pr, pr, pr), 0.0);
        // TODO: The following test is not exact in optimized mode because the
        // compiler chooses to mix 64-bit and 80-bit intermediate results.
        assertDoubleNear(S2.area(pr, pq, pr), 0);
        assertEquals(S2.area(p000, p045, p090), 0.0);

        double maxGirard = 0;
        for (int i = 0; i < 10000; ++i) {
            S2Point p0 = S2Random.randomPoint();
            S2Point d1 = S2Random.randomPoint();
            S2Point d2 = S2Random.randomPoint();
            S2Point p1 = S2Point.plus(p0, S2Point.times(d1, 1e-15));
            S2Point p2 = S2Point.plus(p0, S2Point.times(d2, 1e-15));
            // The actual displacement can be as much as 1.2e-15 due to roundoff.
            // This yields a maximum triangle area of about 0.7e-30.
            assertTrue(S2.area(p0, p1, p2) < 0.7e-30);
            maxGirard = Math.max(maxGirard, S2.girardArea(p0, p1, p2));
        }
        logger.info("Worst case Girard for triangle area 1e-30: " + maxGirard);

        // Try a very long and skinny triangle.
        S2Point p045eps = new S2Point(1, 1, eps);
        double expected2 = 5.8578643762690495119753e-11; // Mathematica.
        assertDoubleNear(S2.area(p000, p045eps, p090), expected2, 1e-9 * expected2);

        // Triangles with near-180 degree edges that sum to a quarter-sphere.
        final double eps2 = 1e-10;
        S2Point p000eps2 = new S2Point(1, 0.1 * eps2, eps2);
        double quarterArea1 =
                S2.area(p000eps2, p000, p090) + S2.area(p000eps2, p090, p180) + S2.area(p000eps2, p180, pz)
                        + S2.area(p000eps2, pz, p000);
        assertDoubleNear(quarterArea1, S2.M_PI);

        // Four other triangles that sum to a quarter-sphere.
        S2Point p045eps2 = new S2Point(1, 1, eps2);
        double quarterArea2 =
                S2.area(p045eps2, p000, p090) + S2.area(p045eps2, p090, p180) + S2.area(p045eps2, p180, pz)
                        + S2.area(p045eps2, pz, p000);
        assertDoubleNear(quarterArea2, S2.M_PI);
    }

    public void testCCW() {
        S2Point a = new S2Point(0.72571927877036835, 0.46058825605889098, 0.51106749730504852);
        S2Point b = new S2Point(0.7257192746638208, 0.46058826573818168, 0.51106749441312738);
        S2Point c = new S2Point(0.72571927671709457, 0.46058826089853633, 0.51106749585908795);
        assertTrue(S2.robustCCW(a, b, c) != 0);
    }

}
