package dilivia.s2.index

// A parameter that controls the reporting of edge intersections.
//
//  - CrossingType::INTERIOR reports intersections that occur at a point
//    interior to both edges (i.e., not at a vertex).
//
//  - CrossingType::ALL reports all intersections, even those where two edges
//    intersect only because they share a common vertex.
enum class CrossingType { INTERIOR, ALL };