package dilivia.s2.collections

import kotlin.math.min


fun <T> List<T>.isSorted(): Boolean where T:Comparable<T>{
    if (this.size <= 1) return true
    val iter = this.iterator()
    var current: T
    var previous = iter.next()
    while (iter.hasNext()) {
        current = iter.next();
        if (previous > current) return false
        previous = current;
    }
    return true
}


/**
 *
 */
fun <V, T: Comparable<V>> List<T>.lowerBound(first: Int = 0, last: Int = this.size, value: V): Int {
    var f = min(first, this.size)
    val l = min(last, this.size)
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
 *
 */
fun <V, T: Comparable<V>> List<T>.upperBound(first: Int = 0, last: Int = this.size, value: V): Int {
    var f = min(first, this.size)
    val l = min(last, this.size)
    var idx: Int
    var step: Int
    var count = l - f

    while (count > 0) {
        idx = f
        step = count / 2
        idx += step
        if (this[idx] <= value) {
            f = ++idx
            count -= step + 1
        }
        else
            count = step
    }
    return f
}
