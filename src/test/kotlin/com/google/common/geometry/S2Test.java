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

    public void testCCW() {
        S2Point a = new S2Point(0.72571927877036835, 0.46058825605889098, 0.51106749730504852);
        S2Point b = new S2Point(0.7257192746638208, 0.46058826573818168, 0.51106749441312738);
        S2Point c = new S2Point(0.72571927671709457, 0.46058826089853633, 0.51106749585908795);
        assertTrue(S2.robustCCW(a, b, c) != 0);
    }

}
