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
package dilivia.s2.math

@Strictfp
class DoubleType() : FloatingPointType<Double> {

    override val zero: Double = 0.0

    override val one: Double = 1.0

    override fun fromDouble(v: Double): Double  = v

    override fun sqrt(v: Double): Double = kotlin.math.sqrt(v)

    override fun plus(a: Double, b: Double): Double = a + b

    override fun minus(a: Double, b: Double): Double = a - b

    override fun inv(a: Double): Double = 1.0 / a

    override fun times(a: Double, b: Double): Double = a * b

    override fun div(a: Double, b: Double): Double = a / b

    override fun abs(a: Double): Double = kotlin.math.abs(a)

    override fun sum(values: List<Double>): Double = values.sum()

    override fun equals(v1: Double, v2: Double): Boolean {
        return v1 == v2
    }

    override fun atan2(y: Double, x: Double): Double = kotlin.math.atan2(y, x)
}