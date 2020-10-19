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

import kotlin.math.IEEErem

@Strictfp
class MutableS1Angle(override var radians: Double = 0.0) : S1Angle(radians) {

    /**
     * Normalize this angle to the range (-180, 180] degrees.
     */
    fun normalize() {
        radians = radians.IEEErem(2.0 * S2.M_PI)
        if (radians <= -S2.M_PI) radians = S2.M_PI
    }

}
