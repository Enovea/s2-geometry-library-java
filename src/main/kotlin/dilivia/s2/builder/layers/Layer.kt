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
package dilivia.s2.builder.layers

import dilivia.s2.S2Error
import dilivia.s2.builder.Graph
import dilivia.s2.builder.GraphOptions


// This class is not needed by ordinary S2Builder clients.  It is only
// necessary if you wish to implement a new S2Builder::Layer subtype.
abstract class Layer {

  // Defines options for building the edge graph that is passed to Build().
  abstract fun graphOptions(): GraphOptions

  // Assembles a graph of snapped edges into the geometry type implemented by
  // this layer.  If an error is encountered, sets "error" appropriately.
  //
  // Note that when there are multiple layers, the Graph objects passed to all
  // layers are guaranteed to be valid until the last Build() method returns.
  // This makes it easier to write algorithms that gather the output graphs
  // from several layers and process them all at once (such as
  // s2builderutil::ClosedSetNormalizer).
  abstract fun build(g: Graph, error: S2Error)

}

