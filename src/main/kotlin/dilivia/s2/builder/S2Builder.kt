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
package dilivia.s2.builder

import com.google.common.collect.ComparisonChain
import dilivia.s2.Assertions.assertEQ
import dilivia.s2.Assertions.assertGE
import dilivia.s2.Assertions.assertLE
import dilivia.s2.Assertions.assertNE
import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2
import dilivia.s2.S2.DBL_EPSILON
import dilivia.s2.S2CellId
import dilivia.s2.S2EdgeCrossings
import dilivia.s2.S2EdgeDistances
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.graph.VertexId
import dilivia.s2.builder.layers.Layer
import dilivia.s2.builder.snap.IdentitySnapFunction
import dilivia.s2.builder.snap.SnapFunction
import dilivia.s2.collections.assign
import dilivia.s2.collections.assignWith
import dilivia.s2.collections.sortAndRemoveDuplicates
import dilivia.s2.index.S2ClosestEdgeQuery
import dilivia.s2.index.S2MinDistance
import dilivia.s2.index.point.S2ClosestPointQuery
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryEdgeTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryPointTarget
import dilivia.s2.index.point.S2ClosestPointQueryBase
import dilivia.s2.index.point.S2PointIndex
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.region.S2Polyline
import dilivia.s2.shape.S2Shape
import dilivia.s2.sin
import mu.KotlinLogging
import kotlin.math.acos
import kotlin.math.max
import kotlin.math.sqrt


//////////////////////  Input Types  /////////////////////////
// All types associated with the S2Builder inputs are prefixed with "Input".

// Identifies an input vertex.
typealias InputVertexId = Int

// Defines an input edge.
typealias InputEdge = Pair<InputVertexId, InputVertexId>

// Identifies an input edge.
typealias InputEdgeId = Int

// Identifies the set of input edge ids that were snapped to a given edge.
typealias InputEdgeIdSetId = Int

// Sort key for prioritizing input vertices.  (Note that keys are *not*
// compared using std::less; see SortInputVertices for details.)
typealias InputVertexKey = Pair<S2CellId, InputVertexId>

//////////////////////  Output Types  /////////////////////////
// These types define the output vertices and edges.

// Identifies a snapped vertex ("snap site").  If there is only one layer,
// than SiteId is the same as Graph::VertexId, but if there are many layers
// then each Graph may contain only a subset of the sites.  Also see
// GraphOptions::allow_vertex_filtering().
typealias SiteId = Int

// Defines an output edge.
typealias Edge = Pair<SiteId, SiteId>

fun min(edge1: Edge, edge2: Edge) = when {
    edge1.first < edge2.first -> edge1
    edge1.first > edge2.first -> edge2
    edge1.second > edge2.second -> edge2
    else -> edge1
}

// Identifies an output edge.
typealias EdgeId = Int

// Identifies an output edge in a particular layer.
typealias LayerEdgeId = Pair<Int, EdgeId>

// Every edge can have a set of non-negative integer labels attached to it.
// When used with an appropriate layer type, you can then retrieve the
// labels associated with each output edge.  This can be useful when merging
// or combining data from several sources.  (Note that in many cases it is
// easier to use separate output layers rather than labels.)
//
// Labels are 32-bit non-negative integers.  To support other label types,
// you can use ValueLexicon to store the set of unique labels seen so far:
//
//   ValueLexicon<MyLabel> my_label_lexicon;
//   builder.set_label(my_label_lexicon.Add(label));
//
// The current set of labels is represented as a stack.  This makes it easy
// to add and remove labels hierarchically (e.g., polygon 5, loop 2).  Use
// set_label() and clear_labels() if you need at most one label per edge.
//
typealias Label = Int

// Each input edge has "label set id" (an int32) representing the set of
// labels attached to that edge.  This vector is populated only if at least
// one label is used.
typealias LabelSetId = Int;

// For output layers that represent polygons, there is an ambiguity inherent
// in spherical geometry that does not exist in planar geometry.  Namely, if
// a polygon has no edges, does it represent the empty polygon (containing
// no points) or the full polygon (containing all points)?  This ambiguity
// also occurs for polygons that consist only of degeneracies, e.g. a
// degenerate loop with only two edges could be either a degenerate shell in
// the empty polygon or a degenerate hole in the full polygon.
//
// To resolve this ambiguity, an IsFullPolygonPredicate may be specified for
// each output layer (see AddIsFullPolygonPredicate below).  If the output
// after snapping consists only of degenerate edges and/or sibling pairs
// (including the case where there are no edges at all), then the layer
// implementation calls the given predicate to determine whether the polygon
// is empty or full except for those degeneracies.  The predicate is given
// an S2Builder::Graph containing the output edges, but note that in general
// the predicate must also have knowledge of the input geometry in order to
// determine the correct result.
//
// This predicate is only needed by layers that are assembled into polygons.
// It is not used by other layer types.
interface IsFullPolygonPredicate {
    fun test(g: Graph, error: S2Error): Boolean
}

// Indicates whether the input edges are undirected.  Typically this is
// specified for each output layer (e.g., s2builderutil::S2PolygonLayer).
//
// Directed edges are preferred, since otherwise the output is ambiguous.
// For example, output polygons may be the *inverse* of the intended result
// (e.g., a polygon intended to represent the world's oceans may instead
// represent the world's land masses).  Directed edges are also somewhat
// more efficient.
//
// However even with undirected edges, most S2Builder layer types try to
// preserve the input edge direction whenever possible.  Generally, edges
// are reversed only when it would yield a simpler output.  For example,
// S2PolygonLayer assumes that polygons created from undirected edges should
// cover at most half of the sphere.  Similarly, S2PolylineVectorLayer
// assembles edges into as few polylines as possible, even if this means
// reversing some of the "undirected" input edges.
//
// For shapes with interiors, directed edges should be oriented so that the
// interior is to the left of all edges.  This means that for a polygon with
// holes, the outer loops ("shells") should be directed counter-clockwise
// while the inner loops ("holes") should be directed clockwise.  Note that
// S2Builder::AddPolygon() follows this convention automatically.
enum class EdgeType { DIRECTED, UNDIRECTED }

// An S2Shape used to represent the entire collection of S2Builder input edges.
// Vertices are specified as indices into a vertex vector to save space.
// Requires that "edges" is constant for the lifetime of this object.
class VertexIdEdgeVectorShape(private val edges: List<Pair<Int, Int>>, private val vertices: List<S2Point>) : S2Shape() {

    fun vertex0(e: Int): S2Point = vertex(edges[e].first)
    fun vertex1(e: Int): S2Point = vertex(edges[e].second)

    // S2Shape interface:

    override val dimension: Int = 1
    override val numEdges: Int = edges.size
    override fun edge(edgeId: Int): dilivia.s2.shape.Edge = dilivia.s2.shape.Edge(vertices[edges[edgeId].first], vertices[edges[edgeId].second])

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

    override val numChains: Int = edges.size
    override fun chain(chain_id: Int): Chain = Chain(chain_id, 1)
    override fun chainEdge(chainId: Int, offset: Int): dilivia.s2.shape.Edge = edge(chainId)

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(edgeId, 0)

    private fun vertex(i: Int): S2Point = vertices[i]

}


//
// This class is a replacement for S2PolygonBuilder.  Once all clients have
// been updated to use this class, S2PolygonBuilder will be removed.


// S2Builder is a tool for assembling polygonal geometry from edges.  Here are
// some of the things it is designed for:
//
// 1. Building polygons, polylines, and polygon meshes from unsorted
//    collections of edges.
//
// 2. Snapping geometry to discrete representations (such as S2CellId centers
//    or E7 lat/lng coordinates) while preserving the input topology and with
//    guaranteed error bounds.
//
// 3. Simplifying geometry (e.g. for indexing, display, or storage).
//
// 4. Importing geometry from other formats, including repairing geometry
//    that has errors.
//
// 5. As a tool for implementing more complex operations such as polygon
//    intersections and unions.
//
// The implementation is based on the framework of "snap rounding".  Unlike
// most snap rounding implementations, S2Builder defines edges as geodesics on
// the sphere (straight lines) and uses the topology of the sphere (i.e.,
// there are no "seams" at the poles or 180th meridian).  The algorithm is
// designed to be 100% robust for arbitrary input geometry.  It offers the
// following properties:
//
//   - Guaranteed bounds on how far input vertices and edges can move during
//     the snapping process (i.e., at most the given "snap_radius").
//
//   - Guaranteed minimum separation between edges and vertices other than
//     their endpoints (similar to the goals of Iterated Snap Rounding).  In
//     other words, edges that do not intersect in the output are guaranteed
//     to have a minimum separation between them.
//
//   - Idempotency (similar to the goals of Stable Snap Rounding), i.e. if the
//     input already meets the output criteria then it will not be modified.
//
//   - Preservation of the input topology (up to the creation of
//     degeneracies).  This means that there exists a continuous deformation
//     from the input to the output such that no vertex crosses an edge.  In
//     other words, self-intersections won't be created, loops won't change
//     orientation, etc.
//
//   - The ability to snap to arbitrary discrete point sets (such as S2CellId
//     centers, E7 lat/lng points on the sphere, or simply a subset of the
//     input vertices), rather than being limited to an integer grid.
//
// Here are some of its other features:
//
//  - It can handle both directed and undirected edges.  Undirected edges can
//    be useful for importing data from other formats, e.g. where loops have
//    unspecified orientations.
//
//  - It can eliminate self-intersections by finding all edge pairs that cross
//    and adding a new vertex at each intersection point.
//
//  - It can simplify polygons to within a specified tolerance.  For example,
//    if two vertices are close enough they will be merged, and if an edge
//    passes nearby a vertex then it will be rerouted through that vertex.
//    Optionally, it can also detect nearly straight chains of short edges and
//    replace them with a single long edge, while maintaining the same
//    accuracy, separation, and topology guarantees ("simplify_edge_chains").
//
//  - It supports many different output types through the concept of "layers"
//    (polylines, polygons, polygon meshes, etc).  You can build multiple
//    layers at once in order to ensure that snapping does not create
//    intersections between different objects (for example, you can simplify a
//    set of contour lines without the risk of having them cross each other).
//
//  - It supports edge labels, which allow you to attach arbitrary information
//    to edges and have it preserved during the snapping process.  (This can
//    also be achieved using layers, at a coarser level of granularity.)
//
// Caveats:
//
//  - Because S2Builder only works with edges, it cannot distinguish between
//    the empty and full polygons.  If your application can generate both the
//    empty and full polygons, you must implement logic outside of this class.
//
// Example showing how to snap a polygon to E7 coordinates:
//
//  using s2builderutil::IntLatLngSnapFunction;
//  S2Builder builder(S2Builder::Options(IntLatLngSnapFunction(7)));
//  S2Polygon output;
//  builder.StartLayer(absl::make_unique<s2builderutil::S2PolygonLayer>(&output));
//  builder.AddPolygon(input);
//  S2Error error;
//  if (!builder.Build(&error)) {
//    S2_LOG(ERROR) << error;
//    ...
//  }
class S2Builder(val options: Options = Options()) {

    //////////// Parameters

    // The maximum distance (inclusive) that a vertex can move when snapped,
    // equal to S1ChordAngle(options_.snap_function().snap_radius()).
    private var siteSnapRadiusCa: S1ChordAngle

    // The maximum distance (inclusive) that an edge can move when snapping to a
    // snap site.  It can be slightly larger than the site snap radius when
    // edges are being split at crossings.
    private var edgeSnapRadiusCa: S1ChordAngle

    private var maxEdgeDeviation: S1Angle
    private var edgeSiteQueryRadiusCa: S1ChordAngle
    private var minEdgeLengthToSplitCa: S1ChordAngle

    private var minSiteSeparation: S1Angle
    private var minSiteSeparationCa: S1ChordAngle
    private var minEdgeSiteSeparationCa: S1ChordAngle
    private var minEdgeSiteSeparationCaLimit: S1ChordAngle

    private var maxAdjacentSiteSeparationCa: S1ChordAngle

    // The squared sine of the edge snap radius.  This is equivalent to the snap
    // radius (squared) for distances measured through the interior of the
    // sphere to the plane containing an edge.  This value is used only when
    // interpolating new points along edges (see GetSeparationSite).
    private var edgeSnapRadiusSin2: Double

    // A copy of the argument to Build().
    private var error: S2Error = S2Error(code = S2Error.OK, text = "")

    // True if snapping was requested.  This is true if either snap_radius() is
    // positive, or split_crossing_edges() is true (which implicitly requests
    // snapping to ensure that both crossing edges are snapped to the
    // intersection point).
    private var snappingRequested: Boolean

    // Initially false, and set to true when it is discovered that at least one
    // input vertex or edge does not meet the output guarantees (e.g., that
    // vertices are separated by at least snap_function.min_vertex_separation).
    private var snappingNeeded: Boolean

    //////////// Input Data /////////////

    // A flag indicating whether label_set_ has been modified since the last
    // time label_set_id_ was computed.
    private var labelSetModified: Boolean

    private var inputVertices: MutableList<S2Point> = mutableListOf()
    private var inputEdges: MutableList<InputEdge> = mutableListOf()

    private var layers: MutableList<Layer> = mutableListOf()
    private var layerOptions: MutableList<GraphOptions> = mutableListOf()
    private var layerBegins: MutableList<InputEdgeId> = mutableListOf()
    private var layerIsFullPolygonPredicates: MutableList<IsFullPolygonPredicate> = mutableListOf()

    private var labelSetIds: MutableList<LabelSetId> = mutableListOf()
    private var label_set_lexicon: IdSetLexicon = IdSetLexicon()

    // The current set of labels (represented as a stack).
    private var labelSet: MutableList<Label> = mutableListOf()

    // The LabelSetId corresponding to the current label set, computed on demand
    // (by adding it to label_set_lexicon()).
    private var labelSetId: LabelSetId

    ////////////// Data for Snapping and Simplifying //////////////

    // The number of sites specified using ForceVertex().  These sites are
    // always at the beginning of the sites_ vector.
    private var num_forced_sites: SiteId = -1

    // The set of snapped vertex locations ("sites").
    private var sites: MutableList<S2Point> = mutableListOf()

    // A map from each input edge to the set of sites "nearby" that edge,
    // defined as the set of sites that are candidates for snapping and/or
    // avoidance.  Note that compact_array will inline up to two sites, which
    // usually takes care of the vast majority of edges.  Sites are kept sorted
    // by increasing distance from the origin of the input edge.
    //
    // Once snapping is finished, this field is discarded unless edge chain
    // simplification was requested, in which case instead the sites are
    // filtered by removing the ones that each edge was snapped to, leaving only
    // the "sites to avoid" (needed for simplification).
    private var edgeSites: MutableList<ArrayList<SiteId>> = mutableListOf()

    // Initializes an S2Builder with the given options.
    init {
        val snapFunction = options.snapFunction
        val snapRadius = snapFunction.snapRadius
        assertLE(snapRadius, SnapFunction.kMaxSnapRadius)

        // Convert the snap radius to an S1ChordAngle.  This is the "true snap
        // radius" used when evaluating exact predicates (s2predicates.h).
        siteSnapRadiusCa = S1ChordAngle(snapRadius)

        // When split_crossing_edges() is true, we need to use a larger snap radius
        // for edges than for vertices to ensure that both edges are snapped to the
        // edge intersection location.  This is because the computed intersection
        // point is not exact; it may be up to kIntersectionError away from its true
        // position.  The computed intersection point might then be snapped to some
        // other vertex up to snap_radius away.  So to ensure that both edges are
        // snapped to a common vertex, we need to increase the snap radius for edges
        // to at least the sum of these two values (calculated conservatively).
        var edgeSnapRadius = snapRadius
        if (!options.simplifyEdgeChains) {
            edgeSnapRadiusCa = siteSnapRadiusCa
        } else {
            edgeSnapRadius += S2EdgeCrossings.kIntersectionError
            edgeSnapRadiusCa = roundUp(edgeSnapRadius)
        }
        snappingRequested = (edgeSnapRadius > S1Angle.zero)

        // Compute the maximum distance that a vertex can be separated from an
        // edge while still affecting how that edge is snapped.
        maxEdgeDeviation = snapFunction.maxEdgeDeviation()
        edgeSiteQueryRadiusCa = S1ChordAngle(maxEdgeDeviation + snapFunction.minEdgeVertexSeparation())

        // Compute the maximum edge length such that even if both endpoints move by
        // the maximum distance allowed (i.e., snap_radius), the center of the edge
        // will still move by less than max_edge_deviation().  This saves us a lot
        // of work since then we don't need to check the actual deviation.
        minEdgeLengthToSplitCa = S1ChordAngle.radians(2 * acos(sin(snapRadius) / sin(maxEdgeDeviation)))

        // If the condition below is violated, then AddExtraSites() needs to be
        // modified to check that snapped edges pass on the same side of each "site
        // to avoid" as the input edge.  Currently it doesn't need to do this
        // because the condition below guarantees that if the snapped edge passes on
        // the wrong side of the site then it is also too close, which will cause a
        // separation site to be added.
        //
        // Currently max_edge_deviation() is at most 1.1 * snap_radius(), whereas
        // min_edge_vertex_separation() is at least 0.219 * snap_radius() (based on
        // S2CellIdSnapFunction, which is currently the worst case).
        assertLE(snapFunction.maxEdgeDeviation(), snapFunction.snapRadius + snapFunction.minEdgeVertexSeparation())

        // To implement idempotency, we check whether the input geometry could
        // possibly be the output of a previous S2Builder invocation.  This involves
        // testing whether any site/site or edge/site pairs are too close together.
        // This is done using exact predicates, which require converting the minimum
        // separation values to an S1ChordAngle.
        minSiteSeparation = snapFunction.minVertexSeparation()
        minSiteSeparationCa = S1ChordAngle(minSiteSeparation)
        minEdgeSiteSeparationCa = S1ChordAngle(snapFunction.minEdgeVertexSeparation())

        // This is an upper bound on the distance computed by S2ClosestPointQuery
        // where the true distance might be less than min_edge_site_separation_ca_.
        minEdgeSiteSeparationCaLimit = addPointToEdgeError(minEdgeSiteSeparationCa)

        // Compute the maximum possible distance between two sites whose Voronoi
        // regions touch.  (The maximum radius of each Voronoi region is
        // edge_snap_radius_.)  Then increase this bound to account for errors.
        maxAdjacentSiteSeparationCa = addPointToPointError(roundUp(edgeSnapRadius * 2.0))

        // Finally, we also precompute sin^2(edge_snap_radius), which is simply the
        // squared distance between a vertex and an edge measured perpendicular to
        // the plane containing the edge, and increase this value by the maximum
        // error in the calculation to compare this distance against the bound.
        val d = sin(edgeSnapRadius)
        edgeSnapRadiusSin2 = d * d;
        edgeSnapRadiusSin2 += ((9.5 * d + 2.5 + 2 * S2.M_SQRT3) * d + 9 * DBL_EPSILON) * DBL_EPSILON

        // Initialize the current label set.
        labelSetId = IdSetLexicon.emptySetId()
        labelSetModified = false

        // If snapping was requested, we try to determine whether the input geometry
        // already meets the output requirements.  This is necessary for
        // idempotency, and can also save work.  If we discover any reason that the
        // input geometry needs to be modified, snapping_needed_ is set to true.
        snappingNeeded = false;
    }

    data class Options(

            // Sets the desired snap function.  The snap function is copied
            // internally, so you can safely pass a temporary object.
            //
            // Note that if your input data includes vertices that were created using
            // S2::GetIntersection(), then you should use a "snap_radius" of
            // at least S2::kIntersectionSnapRadius, e.g. by calling
            //
            //  options.set_snap_function(s2builderutil::IdentitySnapFunction(
            //      S2::kIntersectionSnapRadius));
            //
            // DEFAULT: s2builderutil::IdentitySnapFunction(S1Angle::Zero())
            // [This does no snapping and preserves all input vertices exactly.]
            val snapFunction: SnapFunction = IdentitySnapFunction(S1Angle.zero),

            // If true, then detect all pairs of crossing edges and eliminate them by
            // adding a new vertex at their intersection point.
            //
            // When this option is true, the effective snap_radius() for edges is
            // increased by S2::kIntersectionError to take into account the
            // additional error when computing intersection points.  In other words,
            // edges may move by up to snap_radius() + S2::kIntersectionError.
            //
            // Undirected edges should always be used when the output is a polygon,
            // since splitting a directed loop at a self-intersection converts it into
            // two loops that don't define a consistent interior according to the
            // "interior is on the left" rule.  (On the other hand, it is fine to use
            // directed edges when defining a polygon *mesh* because in that case the
            // input consists of sibling edge pairs.)
            //
            // Self-intersections can also arise when importing data from a 2D
            // projection.  You can minimize this problem by subdividing the input
            // edges so that the S2 edges (which are geodesics) stay close to the
            // original projected edges (which are curves on the sphere).  This can
            // be done using s2builderutil::EdgeSplitter(), for example.
            //
            // DEFAULT: false
            val splitCrossingEdges: Boolean = false,

            // If true, then simplify the output geometry by replacing nearly straight
            // chains of short edges with a single long edge.
            //
            // The combined effect of snapping and simplifying will not change the
            // input by more than the guaranteed tolerances (see the list documented
            // with the SnapFunction class).  For example, simplified edges are
            // guaranteed to pass within snap_radius() of the *original* positions of
            // all vertices that were removed from that edge.  This is a much tighter
            // guarantee than can be achieved by snapping and simplifying separately.
            //
            // However, note that this option does not guarantee idempotency.  In
            // other words, simplifying geometry that has already been simplified once
            // may simplify it further.  (This is unavoidable, since tolerances are
            // measured with respect to the original geometry, which is no longer
            // available when the geometry is simplified a second time.)
            //
            // When the output consists of multiple layers, simplification is
            // guaranteed to be consistent: for example, edge chains are simplified in
            // the same way across layers, and simplification preserves topological
            // relationships between layers (e.g., no crossing edges will be created).
            // Note that edge chains in different layers do not need to be identical
            // (or even have the same number of vertices, etc) in order to be
            // simplified together.  All that is required is that they are close
            // enough together so that the same simplified edge can meet all of their
            // individual snapping guarantees.
            //
            // Note that edge chains are approximated as parametric curves rather than
            // point sets.  This means that if an edge chain backtracks on itself (for
            // example, ABCDEFEDCDEFGH) then such backtracking will be preserved to
            // within snap_radius() (for example, if the preceding point were all in a
            // straight line then the edge chain would be simplified to ACFCFH, noting
            // that C and F have degree > 2 and therefore can't be simplified away).
            //
            // Simplified edges are assigned all labels associated with the edges of
            // the simplified chain.
            //
            // For this option to have any effect, a SnapFunction with a non-zero
            // snap_radius() must be specified.  Also note that vertices specified
            // using ForceVertex are never simplified away.
            //
            // DEFAULT: false
            val simplifyEdgeChains: Boolean = false,

            // If true, then snapping occurs only when the input geometry does not
            // already meet the S2Builder output guarantees (see the SnapFunction
            // class description for details).  This means that if all input vertices
            // are at snapped locations, all vertex pairs are separated by at least
            // min_vertex_separation(), and all edge-vertex pairs are separated by at
            // least min_edge_vertex_separation(), then no snapping is done.
            //
            // If false, then all vertex pairs and edge-vertex pairs closer than
            // "snap_radius" will be considered for snapping.  This can be useful, for
            // example, if you know that your geometry contains errors and you want to
            // make sure that features closer together than "snap_radius" are merged.
            //
            // This option is automatically turned off by simplify_edge_chains(),
            // since simplifying edge chains is never guaranteed to be idempotent.
            //
            // DEFAULT: true
            val idempotent: Boolean = true,

            val verbose: Boolean = false
    )


    // Starts a new output layer.  This method must be called before adding any
    // edges to the S2Builder.  You may call this method multiple times to build
    // multiple geometric objects that are snapped to the same set of sites.
    //
    // For example, if you have a set of contour lines, then you could put each
    // contour line in a separate layer.  This keeps the contour lines separate
    // from each other, while also ensuring that no crossing edges are created
    // when they are snapped and/or simplified.  (This is not true if the
    // contour lines are snapped or simplified independently.)
    //
    // Similarly, if you have a set of polygons that share common boundaries
    // (e.g., countries), you can snap and/or simplify them at the same time by
    // putting them in different layers, while ensuring that their boundaries
    // remain consistent (i.e., no crossing edges or T-vertices are introduced).
    //
    // Ownership of the layer is transferred to the S2Builder.  Example usage:
    //
    // S2Polyline line1, line2;
    // builder.StartLayer(make_unique<s2builderutil::S2PolylineLayer>(&line1)));
    // ... Add edges using builder.AddEdge(), etc ...
    // builder.StartLayer(make_unique<s2builderutil::S2PolylineLayer>(&line2)));
    // ... Add edges using builder.AddEdge(), etc ...
    // S2Error error;
    // S2_CHECK(builder.Build(&error)) << error;  // Builds "line1" & "line2"
    fun startLayer(layer: Layer) {
        logger.debug { "Start layer ${layer.javaClass.simpleName}: ${layer.graphOptions()}" }
        layerOptions.add(layer.graphOptions())
        layerBegins.add(inputEdges.size)
        layerIsFullPolygonPredicates.add(IsFullPolygon(false))
        layers.add(layer)
    }

    // Adds a degenerate edge (representing a point) to the current layer.
    fun addPoint(v: S2Point) {
        addEdge(v, v)
    }

    // Adds the given edge to the current layer.
    fun addEdge(v0: S2Point, v1: S2Point) {
        check(layers.isNotEmpty()) { "Call StartLayer before adding any edges" }
        logger.trace { "Add edge to ${layers.last().javaClass.simpleName}: (${v0.toDegreesString()}, ${v1.toDegreesString()})" }

        if (v0 == v1 && (layerOptions.last().degenerate_edges == DegenerateEdges.DISCARD)) {
            return
        }
        val j0 = addVertex(v0)
        val j1 = addVertex(v1)
        inputEdges.add(InputEdge(j0, j1))

        // If there are any labels, then attach them to this input edge.
        if (labelSetModified) {
            if (labelSetIds.isEmpty()) {
                // Populate the missing entries with empty label sets.
                labelSetIds.assign(inputEdges.size - 1, labelSetId)
            }
            labelSetId = label_set_lexicon.add(labelSet)
            labelSetIds.add(labelSetId)
            labelSetModified = false
        } else if (labelSetIds.isNotEmpty()) {
            labelSetIds.add(labelSetId)
        }
    }

    // Adds the edges in the given polyline.  (Note that if the polyline
    // consists of 0 or 1 vertices, this method does nothing.)
    fun addPolyline(polyline: S2Polyline) {
        val n = polyline.numVertices()
        for (i in 1 until n) {
            addEdge(polyline.vertex(i - 1), polyline.vertex(i))
        }
    }

    // Adds the edges in the given loop.  If the sign() of the loop is negative
    // (i.e. this loop represents a hole within a polygon), the edge directions
    // are automatically reversed to ensure that the polygon interior is always
    // to the left of every edge.
    fun addLoop(loop: S2Loop) {
        logger.debug { "Add Loop to ${layers.last().javaClass.simpleName}:\n${loop.toDebugString()}" }

        // Ignore loops that do not have a boundary.
        if (loop.isEmptyOrFull()) return

        // For loops that represent holes, we add the edge from vertex n-1 to vertex
        // n-2 first.  This is because these edges will be assembled into a
        // clockwise loop, which will eventually be normalized in S2Polygon by
        // calling S2Loop::Invert().  S2Loop::Invert() reverses the order of the
        // vertices, so to end up with the original vertex order (0, 1, ..., n-1) we
        // need to build a clockwise loop with vertex order (n-1, n-2, ..., 0).
        // This is done by adding the edge (n-1, n-2) first, and then ensuring that
        // Build() assembles loops starting from edges in the order they were added.
        val n = loop.numVertices()
        for (i in 0 until n) {
            addEdge(loop.orientedVertex(i), loop.orientedVertex(i + 1))
        }
    }

    // Adds the loops in the given polygon.  Loops representing holes have their
    // edge directions automatically reversed as described for AddLoop().  Note
    // that this method does not distinguish between the empty and full polygons,
    // i.e. adding a full polygon has the same effect as adding an empty one.
    fun addPolygon(polygon: S2Polygon) {
        logger.debug { "Add polygon to ${layers.last().javaClass.simpleName}:\n${polygon.toDebugString(";\n")}" }
        for (i in 0 until polygon.numLoops()) {
            addLoop(polygon.loop(i))
        }

        logger.trace {
            """
                |Polygon added:
                |----------------
                | - input vertices: ${inputVertices.map { p -> p.toDegreesString() }}
                | - input edges: $inputEdges
                | - label set ids: $labelSetIds
                | - layer begins: $layerBegins
                |----------------
            """.trimMargin()
        }
    }

    // Adds the edges of the given shape to the current layer.
    fun addShape(shape: S2Shape) {
        var e = 0
        val n = shape.numEdges
        while (e < n) {
            val edge = shape.edge(e)
            addEdge(edge.v0, edge.v1)
            ++e
        }
    }

    // For layers that are assembled into polygons, this method specifies a
    // predicate that is called when the output consists entirely of degenerate
    // edges and/or sibling pairs.  The predicate is given an S2Builder::Graph
    // containing the output edges (if any) and is responsible for deciding
    // whether this graph represents the empty polygon (possibly with degenerate
    // shells) or the full polygon (possibly with degenerate holes).  Note that
    // this cannot be determined from the output edges alone; it also requires
    // knowledge of the input geometry.  (Also see IsFullPolygonPredicate above.)
    //
    // This method should be called at most once per layer; additional calls
    // simply overwrite the previous value for the current layer.
    //
    // The default predicate simply returns false (i.e., degenerate polygons are
    // assumed to be empty).  Arguably it would better to return an error in
    // this case, but the fact is that relatively few clients need to be able to
    // construct full polygons, and it is unreasonable to expect all such
    // clients to supply an appropriate predicate.
    //
    // The reason for having a predicate rather than a boolean value is that the
    // predicate is responsible for determining whether the output polygon is
    // empty or full.  In general the input geometry is not degenerate, but
    // rather collapses into a degenerate configuration due to snapping and/or
    // simplification.
    //
    // TODO(ericv): Provide standard predicates to handle common cases,
    // e.g. valid input geometry that becomes degenerate due to snapping.
    fun addIsFullPolygonPredicate(predicate: IsFullPolygonPredicate) {
        layerIsFullPolygonPredicates[layerIsFullPolygonPredicates.lastIndex] = predicate
    }

    // A predicate that returns an error indicating that no polygon predicate
    // has been specified.}
    class IsFullPolygonUnspecified : IsFullPolygonPredicate {
        override fun test(g: Graph, error: S2Error): Boolean {
            error.code = S2Error.BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED
            return false
        }
    }

    // Returns a predicate that returns a constant value (true or false);
    class IsFullPolygon(val isFull: Boolean) : IsFullPolygonPredicate {
        override fun test(g: Graph, error: S2Error): Boolean {
            return isFull
        }
    }

    // Forces a vertex to be located at the given position.  This can be used to
    // prevent certain input vertices from moving.  However if you are trying to
    // preserve part of the input boundary, be aware that this option does not
    // prevent edges from being split by new vertices.
    //
    // Forced vertices are never snapped; if this is desired then you need to
    // call options().snap_function().SnapPoint() explicitly.  Forced vertices
    // are also never simplified away (if simplify_edge_chains() is used).
    //
    // Caveat: Since this method can place vertices arbitrarily close together,
    // S2Builder makes no minimum separation guaranteees with forced vertices.
    fun forceVertex(vertex: S2Point) {
        sites.add(vertex)
    }

    // Clear the stack of labels.
    fun clearLabels() {
        labelSet.clear();
        labelSetModified = true
    }

    // Add a label to the stack.
    // REQUIRES: label >= 0.
    fun pushLabel(label: Label) {
        assertGE(label, 0)
        labelSet.add(label)
        labelSetModified = true
    }

    // Remove a label from the stack.
    fun popLabel() {
        labelSet.removeLast()
        labelSetModified = true
    }

    // Convenience function that clears the stack and adds a single label.
    // REQUIRES: label >= 0.
    fun setLabel(label: Label) {
        assertGE(label, 0)
        labelSet.clear()
        labelSet.add(label)
        labelSetModified = true
    }

    // Performs the requested edge splitting, snapping, simplification, etc, and
    // then assembles the resulting edges into the requested output layers.
    //
    // Returns true if all edges were assembled; otherwise sets "error"
    // appropriately.  Depending on the error, some or all output layers may
    // have been created.  Automatically resets the S2Builder state so that it
    // can be reused.
    //
    // REQUIRES: error != nullptr.
    fun build(): S2Error {
        // S2_CHECK rather than S2_DCHECK because this is friendlier than crashing on the
        // "error->ok()" call below.  It would be easy to allow (error == nullptr)
        // by declaring a local "tmp_error", but it seems better to make clients
        // think about error handling.
        error.code = S2Error.OK
        error.text = ""

        // Mark the end of the last layer.
        layerBegins.add(inputEdges.size)

        // See the algorithm overview at the top of this file.
        if (snappingRequested && !options.idempotent) {
            snappingNeeded = true
        }

        logger.trace {
            """
                |Build
                |---------------------------
                | - options: $options
                | - input vertices: ${inputVertices.map { p -> p.toDegreesString() }}
                | - input edges: $inputEdges
                | - label set ids: $labelSetIds
                | - layer begins: $layerBegins
                | - snapping requested: $snappingRequested
                | - snapping needed: $snappingNeeded
                |---------------------------
            """.trimMargin() }

        chooseSites()
        buildLayers()
        reset();
        return error.copy()
    }

    fun build(error: S2Error): Boolean {
        error.init(build())
        return error.isOk()
    }

    // Clears all input data and resets the builder state.  Any options
    // specified are preserved.
    fun reset() {
        inputVertices.clear()
        inputEdges.clear()
        layers.clear()
        layerOptions.clear()
        layerBegins.clear()
        layerIsFullPolygonPredicates.clear()
        labelSetIds.clear()
        label_set_lexicon.clear()
        labelSet.clear()
        labelSetModified = false
        sites.clear()
        edgeSites.clear()
        snappingNeeded = false
    }

    private fun addVertex(v: S2Point): InputVertexId {
        // Remove duplicate vertices that follow the pattern AB, BC, CD.  If we want
        // to do anything more sophisticated, either use a ValueLexicon, or sort the
        // vertices once they have all been added, remove duplicates, and update the
        // edges.
        if (inputVertices.isEmpty() || v != inputVertices.last()) {
            inputVertices.add(v);
        }
        return inputVertices.size - 1
    }

    private fun chooseSites(): Unit {
        if (inputVertices.isEmpty()) return

        val inputEdgeIndex = MutableS2ShapeIndex()
        inputEdgeIndex.add(VertexIdEdgeVectorShape(inputEdges, inputVertices))
        if (options.splitCrossingEdges) {
            addEdgeCrossings(inputEdgeIndex)
        }
        if (snappingRequested) {
            val siteIndex = S2PointIndex<SiteId>()
            addForcedSites(siteIndex)
            chooseInitialSites(siteIndex)
            collectSiteEdges(siteIndex)
        }
        if (snappingNeeded) {
            addExtraSites(inputEdgeIndex)
        } else {
            copyInputEdges()
        }
    }

    private fun copyInputEdges(): Unit {
        // Sort the input vertices, discard duplicates, and update the input edges
        // to refer to the pruned vertex list.  (We sort in the same order used by
        // ChooseInitialSites() to avoid inconsistencies in tests.)
        val sorted = sortInputVertices()
        val vmap = mutableListOf<InputVertexId>()
        sites.clear()
        if (sites is ArrayList) {
            (sites as ArrayList<S2Point>).ensureCapacity(inputVertices.size)
        }
        var i = 0
        while (i < sorted.size) {
            val site = inputVertices[sorted[i].second]
            vmap[sorted[i].second] = sites.size
            while (++i < sorted.size && inputVertices[sorted[i].second] == site) {
                vmap[sorted[i].second] = sites.size
            }
            sites.add(site)
        }
        inputVertices = sites

        for (e in inputEdges.indices) {
            val edge = inputEdges[e]
            inputEdges[e] = Pair(vmap[edge.first], vmap[edge.second])
        }
    }

    private fun sortInputVertices(): List<InputVertexKey> {
        // Sort all the input vertices in the order that we wish to consider them as
        // candidate Voronoi sites.  Any sort order will produce correct output, so
        // we have complete flexibility in choosing the sort key.  We could even
        // leave them unsorted, although this would have the disadvantage that
        // changing the order of the input edges could cause S2Builder to snap to a
        // different set of Voronoi sites.
        //
        // We have chosen to sort them primarily by S2CellId since this improves the
        // performance of many S2Builder phases (due to better spatial locality).
        // It also allows the possibility of replacing the current S2PointIndex
        // approach with a more efficient recursive divide-and-conquer algorithm.
        //
        // However, sorting by leaf S2CellId alone has two small disadvantages in
        // the case where the candidate sites are densely spaced relative to the
        // snap radius (e.g., when using the IdentitySnapFunction, or when snapping
        // to E6/E7 near the poles, or snapping to S2CellId/E6/E7 using a snap
        // radius larger than the minimum value required):
        //
        //  - First, it tends to bias the Voronoi site locations towards points that
        //    are earlier on the S2CellId Hilbert curve.  For example, suppose that
        //    there are two parallel rows of input vertices on opposite sides of the
        //    edge between two large S2Cells, and the rows are separated by less
        //    than the snap radius.  Then only vertices from the cell with the
        //    smaller S2CellId are selected, because they are considered first and
        //    prevent us from selecting the sites from the other cell (because they
        //    are closer than "snap_radius" to an existing site).
        //
        //  - Second, it tends to choose more Voronoi sites than necessary, because
        //    at each step we choose the first site along the Hilbert curve that is
        //    at least "snap_radius" away from all previously selected sites.  This
        //    tends to yield sites whose "coverage discs" overlap quite a bit,
        //    whereas it would be better to cover all the input vertices with a
        //    smaller set of coverage discs that don't overlap as much.  (This is
        //    the "geometric set cover problem", which is NP-hard.)
        //
        // It is not worth going to much trouble to fix these problems, because they
        // really aren't that important (and don't affect the guarantees made by the
        // algorithm), but here are a couple of heuristics that might help:
        //
        // 1. Sort the input vertices by S2CellId at a coarse level (down to cells
        // that are O(snap_radius) in size), and then sort by a fingerprint of the
        // S2Point coordinates (i.e., quasi-randomly).  This would retain most of
        // the advantages of S2CellId sorting, but makes it more likely that we will
        // select sites that are further apart.
        //
        // 2. Rather than choosing the first uncovered input vertex and snapping it
        // to obtain the next Voronoi site, instead look ahead through following
        // candidates in S2CellId order and choose the furthest candidate whose
        // snapped location covers all previous uncovered input vertices.
        //
        // TODO(ericv): Experiment with these approaches.

        val keys = ArrayList<InputVertexKey>(inputVertices.size)
        for (i in 0 until inputVertices.size) {
            keys.add(InputVertexKey(S2CellId.fromPoint(inputVertices[i]), i))
        }

        keys.sortWith { a, b ->
            ComparisonChain.start()
                    .compare(a.first, b.first)
                    .compare(a.second, b.second)
                    .result()
        }

        logger.trace { "Sorted input vertex key: $keys" }
        return keys
    }

    private fun addEdgeCrossings(input_edge_index: MutableS2ShapeIndex): Unit = TODO()

    private fun addForcedSites(siteIndex: S2PointIndex<SiteId>) {
        logger.trace { "Add forced sites to point index: $sites" }
        // Sort the forced sites and remove duplicates.
        sites.sortAndRemoveDuplicates()
        // Add the forced sites to the index.
        for (id in 0 until sites.size) {
            siteIndex.add(sites[id], id)
        }
        num_forced_sites = sites.size
    }

    private fun isForced(v: SiteId): Boolean = v < num_forced_sites

    private fun chooseInitialSites(siteIndex: S2PointIndex<SiteId>) {
        // Find all points whose distance is <= min_site_separation_ca_.
        val options = S2ClosestPointQuery.Options()
        options.setConservativeMaxDistance(minSiteSeparationCa)
        val siteQuery = S2ClosestPointQuery(siteIndex, options)
        val results = mutableListOf<S2ClosestPointQueryBase.Result<S2MinDistance, SiteId>>()

        // Apply the snap_function() to each input vertex, then check whether any
        // existing site is closer than min_vertex_separation().  If not, then add a
        // new site.
        //
        // NOTE(ericv): There are actually two reasonable algorithms, which we call
        // "snap first" (the one above) and "snap last".  The latter checks for each
        // input vertex whether any existing site is closer than snap_radius(), and
        // only then applies the snap_function() and adds a new site.  "Snap last"
        // can yield slightly fewer sites in some cases, but it is also more
        // expensive and can produce surprising results.  For example, if you snap
        // the polyline "0:0, 0:0.7" using IntLatLngSnapFunction(0), the result is
        // "0:0, 0:0" rather than the expected "0:0, 0:1", because the snap radius
        // is approximately sqrt(2) degrees and therefore it is legal to snap both
        // input points to "0:0".  "Snap first" produces "0:0, 0:1" as expected.
        for (key in sortInputVertices()) {
            val vertex = inputVertices[key.second]
            val site = snapSite(vertex)
            // If any vertex moves when snapped, the output cannot be idempotent.
            snappingNeeded = snappingNeeded || site != vertex

            // FindClosestPoints() measures distances conservatively, so we need to
            // recheck the distances using exact predicates.
            //
            // NOTE(ericv): When the snap radius is large compared to the average
            // vertex spacing, we could possibly avoid the call the FindClosestPoints
            // by checking whether sites_.back() is close enough.
            val target = S2ClosestPointQueryPointTarget(site)
            siteQuery.findClosestPoints(target, results)
            var addSite = true
            for (result in results) {
                if (S2Predicates.compareDistance(site, result.point(), minEdgeSiteSeparationCa) <= 0) {
                    addSite = false
                    // This pair of sites is too close.  If the sites are distinct, then
                    // the output cannot be idempotent.
                    snappingNeeded = snappingNeeded || site != result.point()
                }
            }
            if (addSite) {
                logger.trace { "Add site: ${site.toDegreesString()}" }
                siteIndex.add(site, sites.size)
                sites.add(site)
                siteQuery.reInit()
            }
        }
    }

    private fun snapSite(point: S2Point): S2Point {
        if (!snappingRequested) return point
        val site = options.snapFunction.snapPoint(point)
        val distMoved = S1ChordAngle.between(site, point)
        if (distMoved > siteSnapRadiusCa) {
            error.init(S2Error.BUILDER_SNAP_RADIUS_TOO_SMALL,
                    String.format("Snap function moved vertex (%.15g, %.15g, %.15g) by %.15g, which is more than the specified snap radius of %.15g",
                            point.x(), point.y(), point.z(), distMoved.toAngle().radians, siteSnapRadiusCa.toAngle().radians))
        }
        logger.trace { "Snap site ${point.toDegreesString()} => ${site.toDegreesString()}" }
        return site;
    }

    // For each edge, find all sites within min_edge_site_query_radius_ca_ and
    // store them in edge_sites_.  Also, to implement idempotency this method also
    // checks whether the input vertices and edges may already satisfy the output
    // criteria.  If any problems are found then snapping_needed_ is set to true.
    private fun collectSiteEdges(siteIndex: S2PointIndex<SiteId>) {
        // Find all points whose distance is <= edge_site_query_radius_ca_.
        val options = S2ClosestPointQuery.Options()
        options.setConservativeMaxDistance(edgeSiteQueryRadiusCa)
        logger.trace { "options = $options" }
        val siteQuery = S2ClosestPointQuery(siteIndex, options)
        val results = ArrayList<S2ClosestPointQueryBase.Result<S2MinDistance, SiteId>>()
        edgeSites.assignWith(inputEdges.size) { ArrayList() }
        for (e in 0 until inputEdges.size) {
            val edge = inputEdges[e]
            val v0 = inputVertices[edge.first]
            val v1 = inputVertices[edge.second]

            logger.trace { "options = $options" }
            val target = S2ClosestPointQueryEdgeTarget(v0, v1)
            siteQuery.findClosestPoints(target, results)
            logger.trace { "Closest sites from edge $e = (${v0.toDegreesString()}, ${v1.toDegreesString()}): $results" }
            val sites = edgeSites[e]
            sites.ensureCapacity(results.size)
            for (result in results) {
                sites.add(result.data())
                if (!snappingNeeded &&
                        result.distance.value < minEdgeSiteSeparationCaLimit &&
                        result.point() != v0 && result.point() != v1 &&
                        S2Predicates.compareEdgeDistance(result.point(), v0, v1, minEdgeSiteSeparationCa) < 0) {
                    snappingNeeded = true
                }
            }
            sortSitesByDistance(v0, sites)
            println(edgeSites)
        }

        logger.trace { "Collected site edges: $edgeSites" }
    }

    private fun sortSitesByDistance(x: S2Point, sites: MutableList<SiteId>) {
        // Sort sites in increasing order of distance to X.
        sites.sortWith { i, j -> S2Predicates.compareDistances(x, this.sites[i], this.sites[j]) }
    }

    // There are two situatons where we need to add extra Voronoi sites in order
    // to ensure that the snapped edges meet the output requirements:
    //
    //  (1) If a snapped edge deviates from its input edge by more than
    //      max_edge_deviation(), we add a new site on the input edge near the
    //      middle of the snapped edge.  This causes the snapped edge to split
    //      into two pieces, so that it follows the input edge more closely.
    //
    //  (2) If a snapped edge is closer than min_edge_vertex_separation() to any
    //      nearby site (the "site to avoid"), then we add a new site (the
    //      "separation site") on the input edge near the site to avoid.  This
    //      causes the snapped edge to follow the input edge more closely and is
    //      guaranteed to increase the separation to the required distance.
    //
    // We check these conditions by snapping all the input edges to a chain of
    // Voronoi sites and then testing each edge in the chain.  If a site needs to
    // be added, we mark all nearby edges for re-snapping.
    private fun addExtraSites(inputEdgeIndex: MutableS2ShapeIndex) {
        // When options_.split_crossing_edges() is true, this function may be called
        // even when site_snap_radius_ca_ == 0 (because edge_snap_radius_ca_ > 0).
        // However neither of the conditions above apply in that case.
        if (siteSnapRadiusCa == S1ChordAngle.zero) return

        val chain = mutableListOf<SiteId>()  // Temporary
        val snapQueue = mutableListOf<InputEdgeId>()
        for (max_e in 0 until inputEdges.size) {
            snapQueue.add(max_e)
            logger.trace { "AddExtraSites: process edge $max_e => snapQueue = $snapQueue" }
            while (snapQueue.isNotEmpty()) {
                val e = snapQueue.removeLast()
                snapEdge(e, chain)
                // We could save the snapped chain here in a snapped_chains_ vector, to
                // avoid resnapping it in AddSnappedEdges() below, however currently
                // SnapEdge only accounts for less than 5% of the runtime.
                maybeAddExtraSites(e, max_e, chain, inputEdgeIndex, snapQueue)
            }
        }
    }

    private fun maybeAddExtraSites(edgeId: InputEdgeId, maxEdgeId: InputEdgeId, chain: List<SiteId>,
                                   inputEdgeIndex: MutableS2ShapeIndex, snapQueue: MutableList<InputEdgeId>) {
        // The snapped chain is always a *subsequence* of the nearby sites
        // (edge_sites_), so we walk through the two arrays in parallel looking for
        // sites that weren't snapped.  We also keep track of the current snapped
        // edge, since it is the only edge that can be too close.
        var i = 0
        for (id in edgeSites[edgeId]) {
            if (id == chain[i]) {
                if (++i == chain.size) break
                // Check whether this snapped edge deviates too far from its original
                // position.  If so, we split the edge by adding an extra site.
                val v0 = sites[chain[i - 1]];
                val v1 = sites[chain[i]];
                if (S1ChordAngle.between(v0, v1) < minEdgeLengthToSplitCa) continue

                val edge = inputEdges[edgeId]
                val a0 = inputVertices[edge.first]
                val a1 = inputVertices[edge.second]
                if (!S2EdgeDistances.isEdgeBNearEdgeA(a0, a1, v0, v1, maxEdgeDeviation)) {
                    // Add a new site on the input edge, positioned so that it splits the
                    // snapped edge into two approximately equal pieces.  Then we find all
                    // the edges near the new site (including this one) and add them to
                    // the snap queue.
                    //
                    // Note that with large snap radii, it is possible that the snapped
                    // edge wraps around the sphere the "wrong way".  To handle this we
                    // find the preferred split location by projecting both endpoints onto
                    // the input edge and taking their midpoint.
                    val mid = (S2EdgeDistances.project(v0, a0, a1) + S2EdgeDistances.project(v1, a0, a1)).normalize()
                    val newSite = getSeparationSite(mid, v0, v1, edgeId)
                    addExtraSite(newSite, maxEdgeId, inputEdgeIndex, snapQueue)
                    return;
                }
            } else if (i > 0 && id >= num_forced_sites) {
                // Check whether this "site to avoid" is closer to the snapped edge than
                // min_edge_vertex_separation().  Note that this is the only edge of the
                // chain that can be too close because its vertices must span the point
                // where "site_to_avoid" projects onto the input edge XY (this claim
                // relies on the fact that all sites are separated by at least the snap
                // radius).  We don't try to avoid sites added using ForceVertex()
                // because we don't guarantee any minimum separation from such sites.
                val siteToAvoid = sites[id];
                val v0 = sites[chain[i - 1]];
                val v1 = sites[chain[i]];
                if (S2Predicates.compareEdgeDistance(siteToAvoid, v0, v1, minEdgeSiteSeparationCa) < 0) {
                    // A snapped edge can only approach a site too closely when there are
                    // no sites near the input edge near that point.  We fix that by
                    // adding a new site along the input edge (a "separation site"), then
                    // we find all the edges near the new site (including this one) and
                    // add them to the snap queue.
                    val newSite = getSeparationSite(siteToAvoid, v0, v1, edgeId)
                    assertNE(siteToAvoid, newSite)
                    addExtraSite(newSite, maxEdgeId, inputEdgeIndex, snapQueue);
                    return;
                }
            }
        }
    }

    // Adds a new site, then updates "edge_sites"_ for all edges near the new site
    // and adds them to "snap_queue" for resnapping (unless their edge id exceeds
    // "max_edge_id", since those edges have not been snapped the first time yet).
    private fun addExtraSite(newSite: S2Point, maxEdgeId: InputEdgeId, inputEdgeIndex: MutableS2ShapeIndex,
                             snapQueue: MutableList<InputEdgeId>) {
        val newSiteId = sites.size
        sites.add(newSite)
        // Find all edges whose distance is <= edge_site_query_radius_ca_.
        val options = S2ClosestEdgeQuery.Options()
        options.setConservativeMaxDistance(edgeSiteQueryRadiusCa)
        options.includeInteriors = false
        val query = S2ClosestEdgeQuery(inputEdgeIndex, options)
        val target = S2ClosestEdgeQuery.PointTarget(newSite)
        for (result in query.findClosestEdges(target)) {
            val e = result.edgeId
            val siteIds = edgeSites[e]
            siteIds.add(newSiteId)
            sortSitesByDistance(inputVertices[inputEdges[e].first], siteIds)
            if (e <= maxEdgeId) snapQueue.add(e)
        }
    }

    private fun getSeparationSite(siteToAvoid: S2Point, v0: S2Point, v1: S2Point, inputEdgeId: InputEdgeId): S2Point {
        // Define the "coverage disc" of a site S to be the disc centered at S with
        // radius "snap_radius".  Similarly, define the "coverage interval" of S for
        // an edge XY to be the intersection of XY with the coverage disc of S.  The
        // SnapFunction implementations guarantee that the only way that a snapped
        // edge can be closer than min_edge_vertex_separation() to a non-snapped
        // site (i.e., site_to_avoid) if is there is a gap in the coverage of XY
        // near this site.  We can fix this problem simply by adding a new site to
        // fill this gap, located as closely as possible to the site to avoid.
        //
        // To calculate the coverage gap, we look at the two snapped sites on
        // either side of site_to_avoid, and find the endpoints of their coverage
        // intervals.  The we place a new site in the gap, located as closely as
        // possible to the site to avoid.  Note that the new site may move when it
        // is snapped by the snap_function, but it is guaranteed not to move by
        // more than snap_radius and therefore its coverage interval will still
        // intersect the gap.
        val edge = inputEdges[inputEdgeId]
        val x = inputVertices[edge.first]
        val y = inputVertices[edge.second]
        val xyDir = y - x
        val n = S2Point.robustCrossProd(x, y)
        var newSite = S2EdgeDistances.project(siteToAvoid, x, y, n)
        val gapMin = getCoverageEndpoint(v0, x, y, n)
        val gapMax = getCoverageEndpoint(v1, y, x, -n)
        if ((newSite - gapMin).dotProd(xyDir) < 0) {
            newSite = gapMin
        } else if ((gapMax - newSite).dotProd(xyDir) < 0) {
            newSite = gapMax
        }
        newSite = snapSite(newSite)
        assertNE(v0, newSite)
        assertNE(v1, newSite)
        return newSite
    }

    // Given a site P and an edge XY with normal N, intersect XY with the disc of
    // radius snap_radius() around P, and return the intersection point that is
    // further along the edge XY toward Y.
    private fun getCoverageEndpoint(p: S2Point, x: S2Point, y: S2Point, n: S2Point): S2Point {
        // Consider the plane perpendicular to P that cuts off a spherical cap of
        // radius snap_radius().  This plane intersects the plane through the edge
        // XY (perpendicular to N) along a line, and that line intersects the unit
        // sphere at two points Q and R, and we want to return the point R that is
        // further along the edge XY toward Y.
        //
        // Let M be the midpoint of QR.  This is the point along QR that is closest
        // to P.  We can now express R as the sum of two perpendicular vectors OM
        // and MR in the plane XY.  Vector MR is in the direction N x P, while
        // vector OM is in the direction (N x P) x N, where N = X x Y.
        //
        // The length of OM can be found using the Pythagorean theorem on triangle
        // OPM, and the length of MR can be found using the Pythagorean theorem on
        // triangle OMR.
        //
        // In the calculations below, we save some work by scaling all the vectors
        // by n.CrossProd(p).Norm2(), and normalizing at the end.
        val n2 = n.norm2()
        val nDp = n.dotProd(p)
        val nXp = n.crossProd(p)
        val nXpXn = p * n2 - n * nDp
        val om = nXpXn * sqrt(1 - edgeSnapRadiusSin2)
        val mr2 = edgeSnapRadiusSin2 * n2 - nDp * nDp

        // MR is constructed so that it points toward Y (rather than X).
        val mr = nXp * sqrt(max(0.0, mr2))
        return (om + mr).normalize()
    }

    private fun snapEdge(e: InputEdgeId, chain: MutableList<SiteId>) {
        chain.clear()
        val edge = inputEdges[e]
        if (!snappingNeeded) {
            chain.add(edge.first)
            chain.add(edge.second)
            logger.trace { "Snap edge $e: snapping not need, chain = $chain" }
            return
        }

        val x = inputVertices[edge.first]
        val y = inputVertices[edge.second]

        // Optimization: if there is only one nearby site, return.
        // Optimization: if there are exactly two nearby sites, and one is close
        // enough to each vertex, then return.

        // Now iterate through the sites.  We keep track of the sequence of sites
        // that are visited.
        val candidates = edgeSites[e]
        logger.trace { "Snap edge: $e (${x.toDegreesString()}, ${y.toDegreesString()}), candidate = $candidates" }
        for (site_id in candidates) {
            val c = sites[site_id]
            // Skip any sites that are too far away.  (There will be some of these,
            // because we also keep track of "sites to avoid".)  Note that some sites
            // may be close enough to the line containing the edge, but not to the
            // edge itself, so we can just use the dot product with the edge normal.
            if (S2Predicates.compareEdgeDistance(c, x, y, edgeSnapRadiusCa) > 0) {
                logger.trace { "Candidate $site_id ${c.toDegreesString()} is too far from edge $e" }
                continue
            }
            // Check whether the new site C excludes the previous site B.  If so,
            // repeat with the previous site, and so on.
            var addSiteC = true
            while (chain.isNotEmpty()) {
                val b = sites[chain.last()]

                // First, check whether B and C are so far apart that their clipped
                // Voronoi regions can't intersect.
                val bc = S1ChordAngle.between(b, c)
                if (bc >= maxAdjacentSiteSeparationCa) {
                    break
                }

                // Otherwise, we want to check whether site C prevents the Voronoi
                // region of B from intersecting XY, or vice versa.  This can be
                // determined by computing the "coverage interval" (the segment of XY
                // intersected by the coverage disc of radius snap_radius) for each
                // site.  If the coverage interval of one site contains the coverage
                // interval of the other, then the contained site can be excluded.
                val result = S2Predicates.getVoronoiSiteExclusion(b, c, x, y, edgeSnapRadiusCa)
                if (result == S2Predicates.Excluded.FIRST) {
                    chain.removeLast()
                    continue // Site B excluded by C
                }
                if (result == S2Predicates.Excluded.SECOND) {
                    addSiteC = false  // Site C is excluded by B.
                    break
                }
                assertEQ(S2Predicates.Excluded.NEITHER, result)

                // Otherwise check whether the previous site A is close enough to B and
                // C that it might further clip the Voronoi region of B.
                if (chain.size < 2) break
                val a = sites[chain[chain.lastIndex - 1]]
                val ac = S1ChordAngle.between(a, c)
                if (ac >= maxAdjacentSiteSeparationCa) break

                // If triangles ABC and XYB have the same orientation, the circumcenter
                // Z of ABC is guaranteed to be on the same side of XY as B.
                val xyb = S2Predicates.sign(x, y, b)
                if (S2Predicates.sign(a, b, c) == xyb) {
                    break  // The circumcenter is on the same side as B but further away.
                }
                // Other possible optimizations:
                //  - if AB > max_adjacent_site_separation_ca_ then keep B.
                //  - if d(B, XY) < 0.5 * min(AB, BC) then keep B.

                // If the circumcenter of ABC is on the same side of XY as B, then B is
                // excluded by A and C combined.  Otherwise B is needed and we can exit.
                if (S2Predicates.edgeCircumcenterSign(x, y, a, b, c) != xyb) break

            }
            if (addSiteC) {
                chain.add(site_id)
            }
        }
        assert(chain.isNotEmpty())
        /*
        if (google::DEBUG_MODE) {
            for (SiteId site_id : candidates) {
                if (s2pred::CompareDistances(y, sites_[chain->back()],
                sites_[site_id]) > 0) {
                S2_LOG(ERROR) << "Snapping invariant broken!";
            }
            }
        }
         */

        if (options.verbose) {
            print("(${edge.first}, ${edge.second}): ")
            for (id in chain) print("$id ")
            println()
        }
    }

    /** Build the layers. */
    private fun buildLayers() {
        logger.debug {
            """
                |Start build layers:
                |--------------------------
                |Layers:
                |${layers.joinToString(",\n", prefix = " - ")}
            """.trimMargin()
        }

        // Each output edge has an "input edge id set id" (an Int) representing the set of input edge ids that were
        // snapped to this edge.  The actual InputEdgeIds can be retrieved using "inputEdgeIdSetLexicon".
        val layerEdges = ArrayList<ArrayList<Edge>>()
        val layerInputEdgeIds = ArrayList<ArrayList<InputEdgeIdSetId>>()
        val inputEdgeIdSetLexicon = IdSetLexicon()
        buildLayerEdges(layerEdges, layerInputEdgeIds, inputEdgeIdSetLexicon)

        // At this point we have no further need for the input geometry or nearby site data, so we clear those fields to
        // save space.
        inputVertices.clear()
        inputEdges.clear()
        edgeSites.clear()

        // If there are a large number of layers, then we build a minimal subset of vertices for each layer. This ensures
        // that layer types that iterate over vertices will run in time proportional to the size of that layer rather
        // than the size of all layers combined.
        val layerVertices = ArrayList<List<S2Point>>()
        val kMinLayersForVertexFiltering = 10
        if (layers.size >= kMinLayersForVertexFiltering) {
            // Disable vertex filtering if it is disallowed by any layer.  (This could
            // be optimized, but in current applications either all layers allow
            // filtering or none of them do.)
            var allowVertexFiltering = false
            for (options in layerOptions) {
                allowVertexFiltering = allowVertexFiltering and options.allowVertexFiltering
            }
            if (allowVertexFiltering) {
                val filterTmp = mutableListOf<VertexId>()  // Temporary used by FilterVertices.
                layerVertices.ensureCapacity(layers.size)
                for (i in 0 until layers.size) {
                    layerVertices[i] = Graph.filterVertices(sites, layerEdges[i], filterTmp)
                }
                sites.clear()  // Release memory
            }
        }
        for (i in 0 until layers.size) {
            val vertices = if (layerVertices.isEmpty()) sites else layerVertices[i]
            val graph = Graph(layerOptions[i], vertices, layerEdges[i], layerInputEdgeIds[i],
                    inputEdgeIdSetLexicon, labelSetIds, label_set_lexicon, layerIsFullPolygonPredicates[i]
            )
            layers[i].build(graph, error)
            // Don't free the layer data until all layers have been built, in order to
            // support building multiple layers at once (e.g. ClosedSetNormalizer).
        }
    }

    /**
     * Snaps and possibly simplifies the edges for each layer, populating the given output arguments. The resulting
     * edges can be used to construct an Graph directly (no further processing is necessary).
     *
     */
    private fun buildLayerEdges(
            layerEdges: ArrayList<ArrayList<Edge>>,
            layerInputEdgeIds: ArrayList<ArrayList<InputEdgeIdSetId>>,
            inputEdgeIdSetLexicon: IdSetLexicon) {
        // Edge chains are simplified only when a non-zero snap radius is specified.
        // If so, we build a map from each site to the set of input vertices that
        // snapped to that site.
        val siteVertices = ArrayList<MutableList<InputVertexId>>()
        val simplify = snappingNeeded && options.splitCrossingEdges
        if (simplify) siteVertices.ensureCapacity(sites.size)

        layerEdges.assign(layers.size, ArrayList())
        layerInputEdgeIds.assign(layers.size, ArrayList())
        for (i in 0 until layers.size) {
            addSnappedEdges(layerBegins[i], layerBegins[i + 1], layerOptions[i], layerEdges[i], layerInputEdgeIds[i], inputEdgeIdSetLexicon, siteVertices)
        }
        if (simplify) {
            simplifyEdgeChains(siteVertices, layerEdges, layerInputEdgeIds, inputEdgeIdSetLexicon)
        }
        // We simplify edge chains before processing the per-layer GraphOptions
        // because simplification can create duplicate edges and/or sibling edge
        // pairs which may need to be removed.
        for (i in 0 until layers.size) {
            // The errors generated by ProcessEdges are really warnings, so we simply
            // record them and continue.
            Graph.processEdges(layerOptions[i], layerEdges[i], layerInputEdgeIds[i], inputEdgeIdSetLexicon, error)
        }
    }
    // Snaps all the input edges for a given layer, populating the given output
    // arguments.  If (*site_vertices) is non-empty then it is updated so that
    // (*site_vertices)[site] contains a list of all input vertices that were
    // snapped to that site.
    private fun addSnappedEdges(begin: InputEdgeId, end: InputEdgeId, options: GraphOptions, edges: MutableList<Edge>,
                                inputEdgeIds: MutableList<InputEdgeIdSetId>, inputEdgeIdSetLexicon: IdSetLexicon,
                                siteVertices: ArrayList<MutableList<InputVertexId>>) {
        val discardDegenerateEdges = options.degenerate_edges == DegenerateEdges.DISCARD
        val chain = mutableListOf<SiteId>()
        for (e in begin until end) {
            val id = inputEdgeIdSetLexicon.addSingleton(e)
            snapEdge(e, chain)
            maybeAddInputVertex(inputEdges[e].first, chain[0], siteVertices)
            if (chain.size == 1) {
                if (discardDegenerateEdges) continue
                addSnappedEdge(chain[0], chain[0], id, options.edge_type, edges, inputEdgeIds)
            } else {
                maybeAddInputVertex(inputEdges[e].second, chain.last(), siteVertices)
                for (i in 1 until chain.size) {
                    addSnappedEdge(chain[i - 1], chain[i], id, options.edge_type, edges, inputEdgeIds)
                }
            }
        }
        if (this.options.verbose) dumpEdges(edges, sites)
    }

    // If "site_vertices" is non-empty, ensures that (*site_vertices)[id] contains
    // "v".  Duplicate entries are allowed.
    private fun maybeAddInputVertex(v: InputVertexId, id: SiteId, site_vertices: MutableList<MutableList<InputVertexId>>) {
        if (site_vertices.isEmpty()) return

        // Optimization: check if we just added this vertex.  This is worthwhile
        // because the input edges usually form a continuous chain, i.e. the
        // destination of one edge is the same as the source of the next edge.
        val vertices = site_vertices[id]
        if (vertices.isEmpty() || vertices.last() != v) {
            vertices.add(v)
        }
    }

    // Adds the given edge to "edges" and "input_edge_ids".  If undirected edges
    // are being used, also adds an edge in the opposite direction.
    private fun addSnappedEdge(src: SiteId, dst: SiteId, id: InputEdgeIdSetId, edge_type: EdgeType,
                               edges: MutableList<Edge>, input_edge_ids: MutableList<InputEdgeIdSetId>) {
        edges.add(Edge(src, dst))
        input_edge_ids.add(id)
        if (edge_type == EdgeType.UNDIRECTED) {
            edges.add(Edge(dst, src))
            // Automatically created edges do not have input edge ids or labels.  This
            // can be used to distinguish the original direction of the undirected edge.
            input_edge_ids.add(IdSetLexicon.emptySetId())
        }
    }

    private fun simplifyEdgeChains(site_vertices: MutableList<MutableList<InputVertexId>>, layer_edges: ArrayList<ArrayList<Edge>>,
                                   layer_input_edge_ids: ArrayList<ArrayList<InputEdgeIdSetId>>,
                                   input_edge_id_set_lexicon: IdSetLexicon): Unit = TODO()

    private fun mergeLayerEdges(layer_edges: MutableList<MutableList<Edge>>, layer_input_edge_ids: MutableList<MutableList<InputEdgeIdSetId>>,
                                edges: MutableList<Edge>, input_edge_ids: MutableList<InputEdgeIdSetId>,
                                edge_layers: MutableList<Int>): Unit = TODO()


    companion object {

        private val logger = KotlinLogging.logger(S2Builder::class.java.name)

        fun stableLessThan(a: Edge, b: Edge, ai: LayerEdgeId, bi: LayerEdgeId): Boolean = TODO()

        private fun roundUp(a: S1Angle): S1ChordAngle {
            val ca = S1ChordAngle(a)
            return ca.plusError(ca.getS1AngleConstructorMaxError())
        }

        private fun addPointToPointError(ca: S1ChordAngle): S1ChordAngle {
            return ca.plusError(ca.getS2PointConstructorMaxError())
        }

        private fun addPointToEdgeError(ca: S1ChordAngle): S1ChordAngle {
            return ca.plusError(S2EdgeDistances.getUpdateMinDistanceMaxError(ca))
        }

        private fun dumpEdges(edges: MutableList<Pair<SiteId, SiteId>>, vertices: MutableList<S2Point>) {
            edges.map { edge -> "S2Polyline: ${vertices[edge.first]}:${vertices[edge.second]} ($edge)" }
                    .forEach { println(it) }
        }
    }

}


