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

// SequenceLexicon is a class for compactly representing sequences of values
// (e.g., tuples).  It automatically eliminates duplicates, and maps the
// remaining sequences to sequentially increasing integer ids.  See also
// ValueLexicon and IdSetLexicon.
//
// Each distinct sequence is mapped to a 32-bit integer.  The space used for
// each sequence is approximately 11 bytes plus the memory needed to represent
// the sequence elements.  For example, a sequence of three "double"s would
// need about 11 + 3*8 = 35 bytes.  Note also that sequences are referred to
// using 32-bit ids rather than 64-bit pointers.
//
// This class has the same thread-safety properties as "string": const methods
// are thread safe, and non-const methods are not thread safe.
//
// Example usage:
//
//   SequenceLexicon<string> lexicon;
//   vector<string> pets {"cat", "dog", "parrot"};
//   uint32 pets_id = lexicon.Add(pets);
//   S2_CHECK_EQ(pets_id, lexicon.Add(pets));
//   string values;
//   for (const auto& pet : lexicon.sequence(pets_id)) {
//     values += pet;
//   }
//   S2_CHECK_EQ("catdogparrot", values);
//

class SequenceLexicon<T> {

    // A class representing a sequence of values.
    inner class Sequence(val beginIdx: Int, val endIdx: Int) {
        fun begin(): T = values[beginIdx]
        fun end(): T = values[endIdx]
        fun size() = endIdx - beginIdx
        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (javaClass != other?.javaClass) return false

            other as SequenceLexicon<*>.Sequence

            return values() == other.values()
        }

        override fun hashCode(): Int = values().hashCode()

        fun values(): List<T> {
            return values.slice(beginIdx until endIdx)
        }

    };

    private inner class SequenceKey(val id: Int) {

        fun sequence(): Sequence = sequence(id)

        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (javaClass != other?.javaClass) return false

            other as SequenceLexicon<*>.SequenceKey

            return sequence() == other.sequence()
        }

        override fun hashCode(): Int {
            return sequence().hashCode()
        }


    }

    private val values: MutableList<T> = mutableListOf()
    private val begins: MutableList<Int> = mutableListOf(0)
    private val idSet: MutableMap<SequenceKey, Int> = mutableMapOf()

    // Clears all data from the lexicon.
    fun clear() {
        values.clear()
        begins.clear()
        idSet.clear()
        begins.add(0)
    }

    // Add the given sequence of values to the lexicon if it is not already
    // present, and return its integer id.  Ids are assigned sequentially
    // starting from zero.  This is a convenience method equivalent to
    // Add(std::begin(container), std::end(container)).
    fun add(container: Iterable<T>): Int {
        values.addAll(container)
        begins.add(values.size)
        val sequenceKey = SequenceKey(begins.size - 2)
        var id = idSet[sequenceKey]
        if (id == null) {
            id = begins.size - 2
            idSet[sequenceKey] = id
            return id
        } else {
            begins.removeLast()
            while (values.size > begins.last()) values.removeLast()
            return id
        }
    }

    // Return the number of value sequences in the lexicon.
    fun size(): Int = begins.size - 1

    // Return the value sequence with the given id.  This method can be used
    // with range-based for loops as follows:
    //   for (const auto& value : lexicon.sequence(id)) { ... }
    fun sequence(id: Int): Sequence = Sequence(begins[id], begins[id + 1])

}

