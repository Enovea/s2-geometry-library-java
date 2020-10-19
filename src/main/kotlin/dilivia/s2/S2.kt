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

import kotlin.math.sqrt

object S2 {

    // Declare some frequently used constants
    const val M_PI = Math.PI
    const val M_1_PI = 1.0 / Math.PI
    const val M_PI_2 = Math.PI / 2.0
    const val M_PI_4 = Math.PI / 4.0
    val M_SQRT1_2 = sqrt(0.5)
    val M_SQRT2 = sqrt(2.0)
    val M_SQRT3 = sqrt(3.0)
    const val M_E = Math.E


    val DBL_EPSILON: Double by lazy {
        var d = 1.0
        while (1.0 + d / 2 != 1.0) {
            d /= 2.0
        }
        d
    }
    
}
