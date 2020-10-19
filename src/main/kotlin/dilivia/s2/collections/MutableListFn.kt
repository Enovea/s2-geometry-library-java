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
package dilivia.s2.collections

import java.util.ArrayList
import kotlin.math.min


fun <T> MutableList<T>.remove(fromIndex: Int, toIndex: Int) {
    require(fromIndex >= 0 && fromIndex <= this.lastIndex)
    require(toIndex >= fromIndex && toIndex <= this.size)
    val nbElementToRemove = toIndex - fromIndex
    repeat(nbElementToRemove) { this.removeAt(fromIndex) }
}

fun <T> MutableList<T>.assign(size: Int, value: T) {
    if (this is ArrayList) {
        this.ensureCapacity(size)
    }
    (0 until min(this.size, size)).forEach { i -> this[i] = value }
    while(this.size < size) this.add(value)
    while (this.size > size) this.removeLast()
}


fun <T : Comparable<T>> MutableList<T>.sortAndRemoveDuplicates() {
    this.sort()
    var idx = this.lastIndex - 1
    while(idx > 1) {
        if (this[idx] == this[idx - 1]) this.removeAt(idx)
        --idx
    }
}

fun <T> MutableList<T>.reverse(startIdx: Int, endIdx: Int) {
    this.subList(startIdx, endIdx).reverse()
}
