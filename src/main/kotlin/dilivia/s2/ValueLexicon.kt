package dilivia.s2


// ValueLexicon is a class that maps distinct values to sequentially numbered
// integer identifiers.  It automatically eliminates duplicates and uses a
// compact representation.  See also SequenceLexicon.
//
// Each distinct value is mapped to a 32-bit integer.  The space used for each
// value is approximately 7 bytes plus the space needed for the value itself.
// For example, int64 values would need approximately 15 bytes each.  Note
// also that values are referred to using 32-bit ids rather than 64-bit
// pointers.
//
// This class has the same thread-safety properties as "string": const methods
// are thread safe, and non-const methods are not thread safe.
//
// Example usage:
//
//   ValueLexicon<string> lexicon;
//   uint32 cat_id = lexicon.Add("cat");
//   EXPECT_EQ(cat_id, lexicon.Add("cat"));
//   EXPECT_EQ("cat", lexicon.value(cat_id));
//
class ValueLexicon<T> {

    private inner class ValueKey(val id: Int) {

        fun value(): T = values[id]

        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (javaClass != other?.javaClass) return false

            other as ValueLexicon<*>.ValueKey

            val value = value()
            val otherValue = other.value() as T

            if (value != otherValue) return false

            return true
        }

        override fun hashCode(): Int {
            return value().hashCode()
        }


    }

    private val values: MutableList<T> = mutableListOf()
    private val idSet: MutableMap<ValueKey, Int> = mutableMapOf()

    // Clears all data from the lexicon.
    @Synchronized
    fun clear() {
        values.clear()
        idSet.clear()
    }

    // Add the given value to the lexicon if it is not already present, and
    // return its integer id.  Ids are assigned sequentially starting from zero.
    @Synchronized
    fun add(value: T): Int {
        if (values.isNotEmpty() && value == values.last()) {
            return values.lastIndex
        }
        values.add(value)
        val valueKey = ValueKey(values.lastIndex)
        val id = idSet[valueKey]
        return if (id == null) {
            idSet[valueKey] = values.lastIndex
            values.lastIndex
        }
        else {
            values.removeLast()
            id
        }
    }

    // Return the number of values in the lexicon.
    @Synchronized
    fun size(): Int = values.size

    // Return the value with the given id.
    @Synchronized
    fun value(id: Int) = values[id]

}
