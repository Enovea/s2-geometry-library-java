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

class MutableS2Point(override val coords: MutableList<Double>): S2Point(coords) {

    constructor(x: Double, y: Double, z: Double): this(mutableListOf(x, y, z))

    operator fun set(idx: Int, value: Double) {
        coords[idx] = value
    }

    fun x(x: Double) {
        set(0, x)
    }

    fun y(y: Double) {
        set(1, x)
    }

    fun z(z: Double) {
        set(2, x)
    }

    override fun newInstance(coords: List<Double>): S2Point = S2Point(if (coords is MutableList) coords else coords.toMutableList())


}