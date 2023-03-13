/*
This is a web application.  The backend server is written in Go and uses the
html/package to create the html used by the web browser, which points to localhost:8080/fordfulkersonmaxflow.
Ford-Fulkerson finds the max flow between a source vertex and the sink vertex in a st-flow network.
Plot the st-flow network showing the vertices and edges connecting the source and target.
The user enters the number of rows and columns in the flow network.  The number of flow edges is
given by [(rows-1)*3 + 1]*(cols-1) + 2*rows.  The user can change the capacity of each edge.  The capacity
is the maximum flow that the edge can support.  The flow for each edge as well as the total flow is shown.
*/

package main

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"net/http"
	"strconv"
	"strings"
	"text/template"
)

const (
	addr                         = "127.0.0.1:8080"                         // http server listen address
	fileFordFulkersonCapacities  = "templates/fordfulkersoncapacities.html" // html for Ford-Fulkerson capacities
	fileFordFulkersonFlows       = "templates/fordfulkersonflows.html"      // html for Ford-Fulkerson flows
	fileGraphOptions             = "templates/graphoptions.html"            // html for Graph Options
	patternFordFulkersonCapacity = "/fordfulkersoncapacity"                 // http handler for Ford-Fulkerson capacity connections
	patternFordFulkersonFlow     = "/fordfulkersonflow"                     // http handler for Ford-Fulkerson flow connections
	patternGraphOptions          = "/graphoptions"                          // http handler for Graph Options
	rows                         = 300                                      // #rows in grid
	columns                      = rows                                     // #columns in grid
	xlabels                      = 11                                       // # labels on x axis
	ylabels                      = 11                                       // # labels on y axis
	fileEdges                    = "flowedges.csv"                          // flow edge v, w, capacity, flow
	networkColsMax               = 10                                       // flow network limits on vertices and edges
	networkRowsMax               = 10
	networkColsMin               = 2
	networkRowsMin               = 2
)

// FlowEdges are the links between vertices v and w in the s-t flow network
type FlowEdge struct {
	v        int // one vertex
	w        int // the other vertex
	capacity int // maxmimum flow
	flow     int // flow
}

// Type to contain all the HTML template actions
type PlotT struct {
	Grid             []string // plotting grid
	Status           string   // status of the plot
	Xlabel           []string // x-axis labels
	Ylabel           []string // y-axis labels
	Flow             string   // flow in the s-t network
	NetRows          string   // number of rows in the s-t network, 2-10
	NetCols          string   // number of columns in the s-t network, 2-10
	FlowEdges        string   // number of flow edges in the s-t network
	NetFlows         []int    // flows for edges not from/to the source/sink
	SourceFlows      []int    // flows for the source edges
	SinkFlows        []int    // flows for the sink edges
	NetCapacities    [][]int  // capacities for edges not from/to the source/sink
	SourceCapacities []int    // capacities for the source edges
	SinkCapacities   []int    // capacities for the sink edges
}

// FordFulkerson type for Max Flow methods
type FordFulkerson struct {
	edgeTo      []*FlowEdge         // last edge on shortest s->v path
	marked      []bool              // is s->v path in residual graph
	adj         []map[int]*FlowEdge // adjacency list is a list of maps of int to FlowEdge pointers
	flow        int                 // flow in the s-t network
	plot        *PlotT              // data to be distributed in the HTML template
	networkRows int                 // number of rows in the flow network
	networkCols int                 // number of columns in the flow network
}

// global variables for parse and execution of the html template
var (
	tmplFormCapacities *template.Template
	tmplFormFlows      *template.Template
)

// init parses the html template fileS
func init() {
	tmplFormCapacities = template.Must(template.ParseFiles(fileFordFulkersonCapacities))
	tmplFormFlows = template.Must(template.ParseFiles(fileFordFulkersonFlows))
}

// residualCapacityTo returns the residual capacity to vertex
func (fe *FlowEdge) residualCapacityTo(vertex int) int {
	if vertex == fe.v {
		return fe.flow
	} else if vertex == fe.w {
		return fe.capacity - fe.flow
	} else {
		fmt.Printf("residualCapacityTo to vertex %d invalid\n", vertex)
		return math.MaxInt
	}
}

// addResidualFlowTo adds flow to vertex
func (fe *FlowEdge) addResidualFlowTo(vertex int, delta int) {
	if vertex == fe.v {
		fe.flow -= delta
	} else if vertex == fe.w {
		fe.flow += delta
	} else {
		fmt.Printf("addResidualFlow %d to vertix %d invalid\n", delta, vertex)
	}
}

// other returns the vertex connected to vert
func (fe *FlowEdge) other(vert int) int {
	if vert == fe.v {
		return fe.w
	} else {
		return fe.v
	}
}

// findMaxFlow finds the max flow from source to sink vertices in the network
func (ff *FordFulkerson) findMaxFlow(r *http.Request) error {
	// need the rows and columns in the flow network
	networkRows := r.PostFormValue("netrows")
	networkCols := r.PostFormValue("netcols")
	var err error
	if len(networkRows) == 0 || len(networkCols) == 0 {
		return fmt.Errorf("network rows or network columns not set")
	}
	ff.networkRows, err = strconv.Atoi(networkRows)
	if err != nil {
		fmt.Printf("network rows error: %v\n", err)
		return err
	}
	ff.networkCols, err = strconv.Atoi(networkCols)
	if err != nil {
		fmt.Printf("network columns Atoi error: %v\n", err)
		return err
	}

	if ff.networkCols < networkColsMin || ff.networkRows < networkRowsMin ||
		ff.networkCols > networkColsMax || ff.networkRows > networkRowsMax {
		return fmt.Errorf("number of network rows and/or columns are invalid")
	}

	// number of vertices is source + sink + networkRows*networkCols
	numVertices := ff.networkRows*ff.networkCols + 2

	// Create the adjacency list, source vertex is the first one, sink is the last one,
	// row-wise numbering of vertices in the flow network between source and sink
	source := 0
	sink := numVertices - 1

	// make adjacency list of FlowEdges, each vertex gets a map of edges
	ff.adj = make([]map[int]*FlowEdge, numVertices)

	ff.adj[source] = make(map[int]*FlowEdge)
	ff.adj[sink] = make(map[int]*FlowEdge)

	// Connections for source and sink vertex
	for i := source; i < ff.networkRows; i++ {
		otherSource := i*ff.networkCols + 1
		// get the capacity in the webpage, if any
		val := r.PostFormValue(fmt.Sprintf("%d-%d", i, otherSource))
		// default capacity if none assigned by user
		cap := 1
		if len(val) > 0 {
			cap, err = strconv.Atoi(val)
			if err != nil {
				fmt.Errorf("capacity string conversion to int error: %v\n", err)
				cap = 1
			}
		}
		ff.adj[source][otherSource] = &FlowEdge{v: source, w: otherSource, capacity: cap}

		otherSink := otherSource + ff.networkCols - 1
		val = r.PostFormValue(fmt.Sprintf("%d-%d", i, otherSink))
		// default capacity if none assigned by user
		cap = 1
		if len(val) > 0 {
			cap, err = strconv.Atoi(val)
			if err != nil {
				fmt.Errorf("capacity string conversion to int error: %v\n", err)
				cap = 1
			}
		}
		ff.adj[sink][otherSink] = &FlowEdge{v: otherSink, w: sink, capacity: cap}
	}

	// Make the network connections for forward flow direction (to the right)
	for i := source + 1; i < sink; i++ {
		ff.adj[i] = make(map[int]*FlowEdge)
		// first column
		if i%ff.networkCols == 1 {
			val := r.PostFormValue(fmt.Sprintf("%d-%d", i, i+1))
			cap := 1
			if len(val) > 0 {
				cap, err = strconv.Atoi(val)
				if err != nil {
					fmt.Errorf("capacity string conversion to int error: %v\n", err)
					cap = 1
				}
			}
			// right connection
			ff.adj[i][i+1] = &FlowEdge{v: i, w: i + 1, capacity: cap}
			// down right connection
			if i < numVertices-1-ff.networkCols {
				val := r.PostFormValue(fmt.Sprintf("%d-%d", i, i+ff.networkCols+1))
				cap := 1
				if len(val) > 0 {
					cap, err = strconv.Atoi(val)
					if err != nil {
						fmt.Errorf("capacity string conversion to int error: %v\n", err)
						cap = 1
					}
				}
				ff.adj[i][i+ff.networkCols+1] = &FlowEdge{v: i, w: i + ff.networkCols + 1, capacity: cap}
			}
			// up right connection
			if i > ff.networkCols {
				val := r.PostFormValue(fmt.Sprintf("%d-%d", i, i-ff.networkCols+1))
				cap := 1
				if len(val) > 0 {
					cap, err = strconv.Atoi(val)
					if err != nil {
						fmt.Errorf("capacity string conversion to int error: %v\n", err)
						cap = 1
					}
				}
				ff.adj[i][i-ff.networkCols+1] = &FlowEdge{v: i, w: i - ff.networkCols + 1, capacity: cap}
			}
		} else if i%ff.networkCols == 0 { // last column
			// right to sink reference
			ff.adj[i][sink] = ff.adj[sink][i]
		} else { // all other columns
			val := r.PostFormValue(fmt.Sprintf("%d-%d", i, i+1))
			cap := 1
			if len(val) > 0 {
				cap, err = strconv.Atoi(val)
				if err != nil {
					fmt.Errorf("capacity string conversion to int error: %v\n", err)
					cap = 1
				}
			}
			// right connection
			ff.adj[i][i+1] = &FlowEdge{v: i, w: i + 1, capacity: cap}
			// down right connection
			if i < numVertices-1-ff.networkCols {
				val := r.PostFormValue(fmt.Sprintf("%d-%d", i, i+ff.networkCols+1))
				cap := 1
				if len(val) > 0 {
					cap, err = strconv.Atoi(val)
					if err != nil {
						fmt.Errorf("capacity string conversion to int error: %v\n", err)
						cap = 1
					}
				}
				ff.adj[i][i+ff.networkCols+1] = &FlowEdge{v: i, w: i + ff.networkCols + 1, capacity: cap}
			}
			// up right connection
			if i > ff.networkCols {
				val := r.PostFormValue(fmt.Sprintf("%d-%d", i, i-ff.networkCols+1))
				cap := 1
				if len(val) > 0 {
					cap, err = strconv.Atoi(val)
					if err != nil {
						fmt.Errorf("capacity string conversion to int error: %v\n", err)
						cap = 1
					}
				}
				ff.adj[i][i-ff.networkCols+1] = &FlowEdge{v: i, w: i - ff.networkCols + 1, capacity: cap}
			}
		}
	}

	// Make the network connections for backward flow direction, use references
	for i := source + 1; i < sink; i++ {
		ff.adj[i] = make(map[int]*FlowEdge)
		// first column
		if i%ff.networkCols == 1 {
			// left to source reference
			ff.adj[i][source] = ff.adj[source][i]
		} else if i%ff.networkCols == 0 { // last column
			// left connection reference
			ff.adj[i][i-1] = ff.adj[i-1][i]
			// down left connection reference
			if i < numVertices-1-ff.networkCols {
				ff.adj[i][i+ff.networkCols-1] = ff.adj[i+ff.networkCols-1][i]
			}
			// up left connection reference
			if i > ff.networkCols {
				ff.adj[i][i-ff.networkCols-1] = ff.adj[i-ff.networkCols-1][i]
			}
		} else { // all other columns
			// left connection reference
			ff.adj[i][i-1] = ff.adj[i-1][i]
			// down left connection reference
			if i < numVertices-1-ff.networkCols {
				ff.adj[i][i+ff.networkCols-1] = ff.adj[i+ff.networkCols-1][i]
			}
			// up left connection reference
			if i > ff.networkCols {
				ff.adj[i][i-ff.networkCols-1] = ff.adj[i-ff.networkCols-1][i]
			}
		}
	}

	// Find an augmenting path in the residual network via the breadth-first search
	hasAugmentPath := func(s int, t int) bool {
		// initialize marked and edgeTo
		for i := range ff.marked {
			ff.marked[i] = false
		}
		for i := range ff.edgeTo {
			ff.edgeTo[i] = nil
		}

		// BFS queue
		queue := make([]int, 0)

		// Mark the source and put it in the queue
		ff.marked[s] = true
		queue = append(queue, s)
		for len(queue) > 0 {
			v := queue[0]
			queue = queue[1:]
			for _, e := range ff.adj[v] {
				w := e.other(v)
				if e.residualCapacityTo(w) > 0 && !ff.marked[w] {
					ff.edgeTo[w] = e
					ff.marked[w] = true
					queue = append(queue, w)
				}
			}
		}
		return ff.marked[t]
	}

	ff.edgeTo = make([]*FlowEdge, numVertices)
	ff.marked = make([]bool, numVertices)

	// Find max flow in the flow network from source to sink
	for hasAugmentPath(source, sink) {

		// While there exists an augmenting path, use it.

		// Compute the bottleneck capacity
		bottle := math.MaxInt
		for v := sink; v != source; v = ff.edgeTo[v].other(v) {
			result := ff.edgeTo[v].residualCapacityTo(v)
			if result < bottle {
				bottle = result
			}
		}

		// Augment flow
		for v := sink; v != source; v = ff.edgeTo[v].other(v) {
			ff.edgeTo[v].addResidualFlowTo(v, bottle)
		}

		ff.flow += bottle

	}

	return nil
}

// plotFlowNetwork draws the st-flow network in the grid and
func (ff *FordFulkerson) plotFlowNetwork() error {
	// check if the target was found in findSP
	if len(bfsp.distTo) == 0 || bfsp.distTo[bfsp.target] == math.MaxFloat64 {
		return fmt.Errorf("distance to vertex %d not found", bfsp.target)
	}

	var (
		distance float64 = 0.0
		edges    []*Edge = make([]*Edge, 0)
	)

	// Calculate scale factors for x and y
	xscale := (columns - 1) / (bfsp.xmax - bfsp.xmin)
	yscale := (rows - 1) / (bfsp.ymax - bfsp.ymin)

	beginEP := complex(bfsp.xmin, bfsp.ymin) // beginning of the Euclidean graph
	endEP := complex(bfsp.xmax, bfsp.ymax)   // end of the Euclidean graph
	lenEP := cmplx.Abs(endEP - beginEP)      // length of the Euclidean graph

	e := bfsp.edgeTo[bfsp.target]
	if e.w != bfsp.target {
		e.v, e.w = e.w, e.v
	}
	// start at the target and loop until source vertex is plotted to the grid
	nedges := 0
	for {
		v := e.v
		w := e.w

		edges = append(edges, e)

		start := bfsp.location[v]
		end := bfsp.location[w]
		x1 := real(start)
		y1 := imag(start)
		x2 := real(end)
		y2 := imag(end)
		lenEdge := cmplx.Abs(end - start)
		distance += bfsp.graph[v][w]
		ncells := int(columns * lenEdge / lenEP) // number of points to plot in the edge

		deltaX := x2 - x1
		stepX := deltaX / float64(ncells)

		deltaY := y2 - y1
		stepY := deltaY / float64(ncells)

		// loop to draw the SP edge; CSS colors the edge Orange
		x := x1
		y := y1
		for i := 0; i < ncells; i++ {
			row := int((bfsp.ymax-y)*yscale + .5)
			col := int((x-bfsp.xmin)*xscale + .5)
			bfsp.plot.Grid[row*columns+col] = "edgeSP"
			x += stepX
			y += stepY
		}

		// Mark the edge start vertex v.  CSS colors the vertex Black.
		row := int((bfsp.ymax-y1)*yscale + .5)
		col := int((x1-bfsp.xmin)*xscale + .5)
		bfsp.plot.Grid[row*columns+col] = "vertex"

		// Mark the edge end vertex w.  CSS colors the vertex Black.
		row = int((bfsp.ymax-y2)*yscale + .5)
		col = int((x2-bfsp.xmin)*xscale + .5)
		bfsp.plot.Grid[row*columns+col] = "vertex"

		vertices := len(bfsp.location)
		nedges++

		// Exit the loop if source is reached, we have the SP.
		// Or an infinite loop, stop after the number of vertices
		// in the graph is acquired.
		if e.v == bfsp.source || nedges == vertices-1 {
			break
		}

		// move forward to the next edge
		e = bfsp.edgeTo[v]
		if e.w != v {
			e.v, e.w = e.w, e.v
		}
	}

	// Draw the negative-weight edge
	v := bfsp.negEdgeFrom
	w := bfsp.negEdgeTo
	start := bfsp.location[v]
	end := bfsp.location[w]
	x1 := real(start)
	y1 := imag(start)
	x2 := real(end)
	y2 := imag(end)
	lenEdge := cmplx.Abs(end - start)
	ncells := int(columns * lenEdge / lenEP) // number of points to plot in the edge

	deltaX := x2 - x1
	stepX := deltaX / float64(ncells)

	deltaY := y2 - y1
	stepY := deltaY / float64(ncells)

	// loop to draw the edge; CSS colors the cycle edge Yellow
	x := x1
	y := y1
	for i := 0; i < ncells; i++ {
		row := int((bfsp.ymax-y)*yscale + .5)
		col := int((x-bfsp.xmin)*xscale + .5)
		bfsp.plot.Grid[row*columns+col] = "edgeCycle"
		x += stepX
		y += stepY
	}

	// Mark the end vertices of the shortest path
	e = bfsp.edgeTo[bfsp.target]
	x = real(bfsp.location[e.w])
	y = imag(bfsp.location[e.w])
	// Mark the SP end vertex.  CSS colors the vertex Red.
	row := int((bfsp.ymax-y)*yscale + .5)
	col := int((x-bfsp.xmin)*xscale + .5)
	bfsp.plot.Grid[row*columns+col] = "vertexSP2"
	bfsp.plot.Grid[(row+1)*columns+col] = "vertexSP2"
	bfsp.plot.Grid[(row-1)*columns+col] = "vertexSP2"
	bfsp.plot.Grid[row*columns+col+1] = "vertexSP2"
	bfsp.plot.Grid[row*columns+col-1] = "vertexSP2"

	bfsp.plot.TargetLocation = fmt.Sprintf("(%.2f, %.2f)", x, y)
	bfsp.plot.Target = strconv.Itoa(e.w)

	// Mark the SP start vertex.  CSS colors the vertex Blue.
	x = real(bfsp.location[bfsp.source])
	y = imag(bfsp.location[bfsp.source])
	row = int((bfsp.ymax-y)*yscale + .5)
	col = int((x-bfsp.xmin)*xscale + .5)
	bfsp.plot.Grid[row*columns+col] = "vertexSP1"
	bfsp.plot.Grid[(row+1)*columns+col] = "vertexSP1"
	bfsp.plot.Grid[(row-1)*columns+col] = "vertexSP1"
	bfsp.plot.Grid[row*columns+col+1] = "vertexSP1"
	bfsp.plot.Grid[row*columns+col-1] = "vertexSP1"

	bfsp.plot.SourceLocation = fmt.Sprintf("(%.2f, %.2f)", x, y)
	bfsp.plot.Source = strconv.Itoa(bfsp.source)

	// Distance of the SP
	bfsp.plot.DistanceSP = fmt.Sprintf("%.2f", distance)
	if distance < 0.0 {
		bfsp.plot.NegDistance = "negativedistance"
	}

	// Negative-edge vertices and weight
	bfsp.plot.NegativeEdgeFrom = strconv.Itoa(bfsp.negEdgeFrom)
	bfsp.plot.NegativeEdgeTo = strconv.Itoa(bfsp.negEdgeTo)
	bfsp.plot.NegativeWeight = fmt.Sprintf("%.2f", bfsp.negWeight)

	// list of vertices on the Shortest Path
	n := len(edges)
	verts := make([]string, n+1)
	for i := 0; i < n; i++ {
		verts[n-i] = strconv.Itoa(edges[i].w)
	}
	verts[0] = strconv.Itoa(edges[n-1].v)

	bfsp.plot.PathSP = strings.Join(verts, ", ")

	return nil

}

// HTTP handler for /graphoptions connections
func handleGraphOptions(w http.ResponseWriter, r *http.Request) {
	http.ServeFile(w, r, "templates/graphoptions.html")
}

// showCapacities displays the FlowEdge capacities
func (ff *FordFulkerson) showCapacities() error {
	return nil
}

// HTTP handler for /fordfulkersoncapacity connections
func handleFordFulkersonCapacity(w http.ResponseWriter, r *http.Request) {

	// Create the Bellman Ford SP instance
	fordFulkerson := &FordFulkerson{}

	// Accumulate error
	status := make([]string, 0)

	// Find the Max Flow with shortest-augmenting path
	err := fordFulkerson.findMaxFlow(r)
	if err != nil {
		fmt.Printf("findSP error: %v\n", err)
		status = append(status, err.Error())
	}

	// Draw Flow Network into 300 x 300 cell 2px grid
	err = fordFulkerson.plotFlowNetwork()
	if err != nil {
		fmt.Printf("plotSP error: %v\n", err)
		status = append(status, err.Error())
	}

	// Show the Flow Capacities of each FlowEdge
	err = fordFulkerson.showCapacities()
	if err != nil {
		fmt.Printf("showCapacities error: %v\n", err)
		status = append(status, err.Error())
	}

	// Status
	if len(status) > 0 {
		fordFulkerson.plot.Status = strings.Join(status, ", ")
	} else {
		fordFulkerson.plot.Status = "Enter Source and Target Vertices (0-V-1) for another SP"
	}

	// Write to HTTP using template and grid
	if err := tmplFormCapacities.Execute(w, fordFulkerson.plot); err != nil {
		log.Fatalf("Write to HTTP output using template with grid error: %v\n", err)
	}
}

// HTTP handler for /fordfulkersonflow connections
func handleFordFulkersonFlow(w http.ResponseWriter, r *http.Request) {

}

// main sets up the http handlers, listens, and serves http clients
func main() {
	// Set up http servers with handler for Graph Options and Ford-Fulkerson Max Flow
	http.HandleFunc(patternFordFulkersonCapacity, handleFordFulkersonCapacity)
	http.HandleFunc(patternFordFulkersonFlow, handleFordFulkersonFlow)
	http.HandleFunc(patternGraphOptions, handleGraphOptions)
	fmt.Printf("Ford-Fulkerson Max Flow Server listening on %v.\n", addr)
	http.ListenAndServe(addr, nil)
}
