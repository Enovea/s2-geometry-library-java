package dilivia.s2.collections

fun <A, B> Pair<A, B>.reverse(): Pair<B, A> = Pair(this.second, this.first)
