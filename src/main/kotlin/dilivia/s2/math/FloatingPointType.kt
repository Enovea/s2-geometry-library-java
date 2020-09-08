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

interface FloatingPointType<T: Number> {

    val zero: T

    val one: T

    fun fromDouble(v: Double): T

    fun sqrt(v: T): T

    fun plus(a: T, b: T): T

    fun minus(a: T, b: T): T

    fun inv(a: T): T

    fun times(a: T, b: T): T

    fun div(a: T, b: T): T

    fun abs(a: T): T

    fun sum(values: List<T>): T

    fun equals(v1: T, v2: T): Boolean

    fun atan2(y: T, x: T): T

}