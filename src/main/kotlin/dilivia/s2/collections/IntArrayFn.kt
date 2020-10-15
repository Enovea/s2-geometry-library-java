package dilivia.s2.collections


/**
 * Gets the index of the first element in index range [first, last) that is not less (i.e. greater or equal to) value,
 * org last if no such element is found.
 *
 * @param first The first index of the research range.
 * @param last The last index (exclusive) of the research range.
 * @param value The value to compare the elements to.
 * @return The index of the first not less element.
 */
fun IntArray.lowerBound(first: Int = 0, last: Int = this.size, value: Int): Int {
    var f = kotlin.math.min(first, this.size)
    val l = kotlin.math.min(last, this.size)
    var idx: Int
    var step: Int
    var count = l - f

    while (count > 0) {
        step = count / 2
        idx = f + step
        if (this[idx] < value) {
            f = ++idx
            count -= step + 1
        }
        else
            count = step
    }
    return f
}


/**
 * Gets the index of the first element in range [first, last) that is greater than value, or last if no such element
 * is found.
 *
 * @param first The first index of the research range.
 * @param last The last index (exclusive) of the research range.
 * @param value The value to compare the elements to.
 * @return The index of the first greater element.
 */
fun IntArray.upperBound(first: Int, last: Int, value: Int): Int {
    var f = kotlin.math.min(first, this.size)
    val l = kotlin.math.min(last, this.size)
    var idx: Int
    var step: Int
    var count = l - f

    while (count > 0) {
        idx = f
        step = count / 2
        idx += step
        if (value >= this[idx]) {
            f = ++idx
            count -= step + 1
        }
        else
            count = step
    }
    return f
}
