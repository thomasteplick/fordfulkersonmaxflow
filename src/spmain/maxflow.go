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
	"bufio"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"net/http"
	"os"
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
	fileFlowGraph                = "flowgraph.csv" // rows, columns, and edges
)

// FlowEdges are the links between vertices v and w in the s-t flow network
type FlowEdge struct {
	v        int // one vertex
	w        int // the other vertex
	capacity int // maxmimum flow
	flow     int // flow
}

// Each cell has a name and a value to describe the edge flow properties
type Cell struct {
	Name    string // vertices of the flow edge, v-w
	Value   string // integer edge flow or capacity
	BGColor string // background color, red=no flow, green=flow
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
	NetFlows         [][]Cell // flows for edges not from/to the source/sink
	SourceFlows      []Cell   // flows for the source edges
	SinkFlows        []Cell   // flows for the sink edges
	NetCapacities    [][]Cell // names and capacities for edges not from/to the source/sink
	SourceCapacities []Cell   // names and capacities for the source edges
	SinkCapacities   []Cell   // names and capacities for the sink edges
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

// addResidualFlowTo adds flow delta to vertex
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

// getMaxFlowResults reads the MaxFlow data from disk
func (ff *FordFulkerson) getMaxFlowResults() error {
	// Open the .csv to get the flow graph data for this row/column configuration
	f, err := os.Open(fileFlowGraph)
	if err != nil {
		fmt.Printf("Open file %s error: %v\n", fileFlowGraph, err)
	}
	defer f.Close()
	input := bufio.NewScanner(f)
	input.Scan()
	line := input.Text()
	// Each line has comma-separated values
	values := strings.Split(line, ",")
	var rows, cols int
	if rows, err = strconv.Atoi(values[0]); err != nil {
		fmt.Printf("String %s conversion to int error: %v\n", values[0], err)
		return err
	}

	if cols, err = strconv.Atoi(values[1]); err != nil {
		fmt.Printf("String %s conversion to int error: %v\n", values[1], err)
		return err
	}

	ff.networkRows = rows
	ff.networkCols = cols

	// number of vertices is source + sink + networkRows*networkCols
	numVertices := ff.networkRows*ff.networkCols + 2

	// make adjacency list of FlowEdges, each vertex gets a map of edges
	ff.adj = make([]map[int]*FlowEdge, numVertices)
	for i := range ff.adj {
		ff.adj[i] = make(map[int]*FlowEdge)
	}

	sinkFlow := 0

	for input.Scan() {
		line := input.Text()
		// Each line has comma-separated values
		values := strings.Split(line, ",")
		var v, w, cap, flow int
		if v, err = strconv.Atoi(values[0]); err != nil {
			fmt.Printf("String %s conversion to int error: %v\n", values[0], err)
			continue
		}
		if w, err = strconv.Atoi(values[1]); err != nil {
			fmt.Printf("String %s conversion to int error: %v\n", values[1], err)
			continue
		}
		if cap, err = strconv.Atoi(values[2]); err != nil {
			fmt.Printf("String %s conversion to int error: %v\n", values[2], err)
			continue
		}
		if flow, err = strconv.Atoi(values[3]); err != nil {
			fmt.Printf("String %s conversion to int error: %v\n", values[3], err)
			continue
		}
		ff.adj[v][w] = &FlowEdge{v: v, w: w, capacity: cap, flow: flow}
		ff.adj[w][v] = ff.adj[v][w]
		// get flow delivered from source which equals flow delivered to sink
		// this is the max flow in the network
		if v == 0 {
			ff.flow += flow
		} else if w == numVertices-1 {
			sinkFlow += flow
		}
	}
	if ff.flow != sinkFlow {
		fmt.Printf("source flow doesn't equal sink flow in %s\n", fileFlowGraph)
	}

	return nil
}

// findMaxFlow finds the max flow from source to sink vertices in the network
func (ff *FordFulkerson) findMaxFlow(r *http.Request) error {
	// need the rows and columns in the flow network
	networkRows := r.PostFormValue("netrows")
	networkCols := r.PostFormValue("netcolumns")

	// If we got here using links then a form was not submitted and
	// these controls will not be present.  Get the saved max flow data.
	if len(networkRows) == 0 || len(networkCols) == 0 {
		fmt.Printf("Read the previously calculated MaxFlow data from disk.\n")
		// Get the previously calculated MaxFlow data from disk
		err := ff.getMaxFlowResults()
		if err != nil {
			fmt.Printf("getMaxFlowResults error: %v\n", err)
			return err
		}
		return nil
	}

	var (
		err error
	)

	// Create the flow graph for this row/column configuration

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
		val := r.PostFormValue(fmt.Sprintf("%d-%d", source, otherSource))
		// default capacity if none assigned by user to maximize flow
		cap := 3
		if i == source || i == ff.networkRows-1 {
			cap = 2
		}
		if len(val) > 0 {
			cap, err = strconv.Atoi(val)
			if err != nil {
				fmt.Printf("capacity string conversion to int error: %v\n", err)
				cap = 1
			}
		}
		ff.adj[source][otherSource] = &FlowEdge{v: source, w: otherSource, capacity: cap}

		otherSink := otherSource + ff.networkCols - 1
		val = r.PostFormValue(fmt.Sprintf("%d-%d", sink, otherSink))
		// default capacity if none assigned by user to maximize the flow
		cap = 3
		if i == source || i == ff.networkRows-1 {
			cap = 2
		}
		if len(val) > 0 {
			cap, err = strconv.Atoi(val)
			if err != nil {
				fmt.Printf("capacity string conversion to int error: %v\n", err)
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
					fmt.Printf("capacity string conversion to int error: %v\n", err)
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
						fmt.Printf("capacity string conversion to int error: %v\n", err)
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
						fmt.Printf("capacity string conversion to int error: %v\n", err)
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
					fmt.Printf("capacity string conversion to int error: %v\n", err)
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
						fmt.Printf("capacity string conversion to int error: %v\n", err)
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
						fmt.Printf("capacity string conversion to int error: %v\n", err)
						cap = 1
					}
				}
				ff.adj[i][i-ff.networkCols+1] = &FlowEdge{v: i, w: i - ff.networkCols + 1, capacity: cap}
			}
		}
	}

	// Make the network connections for backward flow direction, use references
	for i := source + 1; i < sink; i++ {
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

	ff.edgeTo = make([]*FlowEdge, numVertices)
	ff.marked = make([]bool, numVertices)

	// Find an augmenting path in the residual network using a breadth-first search
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

	// Save the rows, columns, and edges to a csv file
	f, err := os.Create(fileFlowGraph)
	if err != nil {
		fmt.Printf("Create file %s error: %v\n", fileFlowGraph, err)
		return err
	}
	defer f.Close()
	// Save the rows and columns
	fmt.Fprintf(f, "%d,%d\n", ff.networkRows, ff.networkCols)
	// Save the flow edges, loop over the list of maps and only save
	// edges coming from this vertex (forward flows)
	for vert, mp := range ff.adj {
		// loop over the edges for this vertex
		for _, e := range mp {
			if e.v == vert {
				fmt.Fprintf(f, "%d,%d,%d,%d\n", e.v, e.w, e.capacity, e.flow)
			}
		}
	}

	return nil
}

// location returns the complex coordinate of the vertex
func (ff *FordFulkerson) location(vert int) complex128 {
	// source
	if vert == 0 {
		return complex(0.0, float64(ff.networkRows+1)/2.0)
	} else if vert == ff.networkCols*ff.networkRows+1 { // sink
		return complex(float64(ff.networkCols+1), float64(ff.networkRows+1)/2.0)
	} else {
		return complex(float64((vert-1)%ff.networkCols+1), float64(ff.networkRows-((vert-1)/ff.networkCols)))
	}
}

// plotFlowNetwork draws the st-flow network in the grid
func (ff *FordFulkerson) plotFlowNetwork() error {

	ff.plot = &PlotT{}
	ff.plot.Grid = make([]string, rows*columns)
	ff.plot.Xlabel = make([]string, xlabels)
	ff.plot.Ylabel = make([]string, ylabels)

	// Calculate scale factors for x and y, x axis consists of source, sink, network columns
	xscale := float64(columns-1) / float64(ff.networkCols+1)
	yscale := float64(rows-1) / float64(ff.networkRows)

	// Insert the flow graph vertices and edges in the grid
	// loop over the flow graph vertices

	// color the vertices black
	// color the edges connecting the vertices gray
	// create the line y = mx + b for each edge
	// translate complex coordinates to row/col on the grid
	// translate row/col to slice data object []string Grid
	// CSS selectors for background-color are "vertex" and "edge"

	beginEP := complex(0, 0)
	endEP := complex(float64(ff.networkCols+1), float64(ff.networkRows))
	lenEP := cmplx.Abs(endEP - beginEP)

	for i, mp := range ff.adj {

		// Insert the edge between the vertices v, w.  Do this before marking the vertices.
		// CSS colors the edge according to its flow content:  gray is no flow, green has flow

		for _, e := range mp {
			// Plot forward flow edges only
			if i != e.v {
				continue
			}

			// Get complex coordinates of the vertex endpoints v and w
			beginEdge := ff.location(i)
			endEdge := ff.location(e.w)
			lenEdge := cmplx.Abs(endEdge - beginEdge)
			ncells := int(columns * lenEdge / lenEP) // number of points to plot in the edge

			beginX := real(beginEdge)
			endX := real(endEdge)
			deltaX := endX - beginX
			stepX := deltaX / float64(ncells)

			beginY := imag(beginEdge)
			endY := imag(endEdge)
			deltaY := endY - beginY
			stepY := deltaY / float64(ncells)

			// loop to draw the edge
			x := beginX
			y := beginY
			color := "edge"
			if e.flow > 0 {
				color = "edgeflow"
			}
			for i := 0; i < ncells; i++ {
				row := int((float64(ff.networkRows)-y)*yscale + .5)
				col := int((x-0.0)*xscale + .5)
				ff.plot.Grid[row*columns+col] = color
				x += stepX
				y += stepY
			}

			// Mark the edge start vertex v.  CSS colors the vertex black.
			row := int((float64(ff.networkRows)-beginY)*yscale + .5)
			col := int((beginX-0.0)*xscale + .5)
			ff.plot.Grid[row*columns+col] = "vertex"

			// Mark the edge end vertex w.  CSS colors the vertex black.
			row = int((float64(ff.networkRows)-endY)*yscale + .5)
			col = int((endX-0.0)*xscale + .5)
			ff.plot.Grid[row*columns+col] = "vertex"
		}
	}

	// Construct x-axis labels
	incr := float64(ff.networkCols+2) / (xlabels - 1)
	x := 0.0
	// First label is empty for alignment purposes
	for i := range ff.plot.Xlabel {
		ff.plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Construct the y-axis labels
	incr = float64(ff.networkRows) / (ylabels - 1)
	y := 0.0
	for i := range ff.plot.Ylabel {
		ff.plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// show number of net rows and columns
	ff.plot.NetRows = strconv.Itoa(ff.networkRows)
	ff.plot.NetCols = strconv.Itoa(ff.networkCols)

	// show number of flow edges
	ff.plot.FlowEdges = strconv.Itoa(((ff.networkRows-1)*3+1)*(ff.networkCols-1) + 2*ff.networkRows)

	return nil
}

// HTTP handler for /graphoptions connections
func handleGraphOptions(w http.ResponseWriter, r *http.Request) {
	http.ServeFile(w, r, "templates/graphoptions.html")
}

// showCapacities displays the FlowEdge capacities in HTML tables
func (ff *FordFulkerson) showCapacities() error {
	// Construct three separate HTML tables for source, sink, and the network
	// Insert the flow capacities for the respective edges in PlotT object and
	// apply to the HTML template

	ff.plot.SourceCapacities = make([]Cell, ff.networkRows)
	// source table with vertex = 0
	for i := 0; i < ff.networkRows; i++ {
		ff.plot.SourceCapacities[i] = Cell{Value: strconv.Itoa(ff.adj[0][1+i*ff.networkCols].capacity),
			Name: fmt.Sprintf("%d-%d", 0, 1+i*ff.networkCols)}
	}

	ff.plot.SinkCapacities = make([]Cell, ff.networkRows)
	//  sink table with vertex = f.networkRows*ff.networkCols+1
	for i := 0; i < ff.networkRows; i++ {
		ff.plot.SinkCapacities[i] = Cell{Value: strconv.Itoa(ff.adj[ff.networkRows*ff.networkCols+1][i*ff.networkCols+ff.networkCols].capacity),
			Name: fmt.Sprintf("%d-%d", ff.networkRows*ff.networkCols+1, i*ff.networkCols+ff.networkCols)}
	}

	// network table forward flow edges only (flow leaving the vertex)
	ff.plot.NetCapacities = make([][]Cell, (ff.networkRows-1)*3+1)
	for row := 0; row < (ff.networkRows-1)*3+1; row++ {
		// one less edge than the number of vertices in each row
		ff.plot.NetCapacities[row] = make([]Cell, ff.networkCols-1)
	}
	vStart := 1
	for col := 0; col < ff.networkCols-1; col++ {
		v := vStart
		for row := 0; row < (ff.networkRows-1)*3+1; {
			// flow edges for the first network row
			if v == vStart {
				ff.plot.NetCapacities[row][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1].capacity), Name: fmt.Sprintf("%d-%d", v, v+1)}
				ff.plot.NetCapacities[row+1][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1+ff.networkCols].capacity),
					Name: fmt.Sprintf("%d-%d", v, v+1+ff.networkCols)}
				v += ff.networkCols
				row += 2
			} else if v == vStart+(ff.networkRows-1)*ff.networkCols { // flow edges for the last network row
				ff.plot.NetCapacities[row][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1].capacity), Name: fmt.Sprintf("%d-%d", v, v+1)}
				ff.plot.NetCapacities[row+1][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1-ff.networkCols].capacity),
					Name: fmt.Sprintf("%d-%d", v, v+1-ff.networkCols)}
				v += ff.networkCols
				row += 2
			} else { // all the other network rows
				ff.plot.NetCapacities[row][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1].capacity), Name: fmt.Sprintf("%d-%d", v, v+1)}
				ff.plot.NetCapacities[row+1][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1+ff.networkCols].capacity),
					Name: fmt.Sprintf("%d-%d", v, v+1+ff.networkCols)}
				ff.plot.NetCapacities[row+2][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1-ff.networkCols].capacity),
					Name: fmt.Sprintf("%d-%d", v, v+1-ff.networkCols)}
				v += ff.networkCols
				row += 3
			}
		}
		vStart++
	}
	return nil
}

// HTTP handler for /fordfulkersoncapacity connections
func handleFordFulkersonCapacity(w http.ResponseWriter, r *http.Request) {

	// Create the Ford Fulkerson maxflow instance
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
		fordFulkerson.plot.Status = "Click 'Flows' link to see edge flows"
	}

	// Write to HTTP using template and grid
	if err := tmplFormCapacities.Execute(w, fordFulkerson.plot); err != nil {
		log.Fatalf("Write to HTTP output using template with grid error: %v\n", err)
	}
}

// showFlows displays the FlowEdge flows in HTML tables
func (ff *FordFulkerson) showFlows() error {
	// Construct three separate HTML tables for source, sink, and the network
	// Insert the flows for the respective edges in PlotT object and
	// apply to the HTML template

	ff.plot.SourceFlows = make([]Cell, ff.networkRows)
	// source table with vertex = 0,
	for i := 0; i < ff.networkRows; i++ {
		v := 0
		w := 1 + i*ff.networkCols
		flow := ff.adj[v][w].flow
		color := "edgenoflow"
		if flow > 0 {
			color = "edgeflow"
		}
		ff.plot.SourceFlows[i] = Cell{Value: strconv.Itoa(flow),
			Name: fmt.Sprintf("%d-%d", v, w), BGColor: color}
	}

	ff.plot.SinkFlows = make([]Cell, ff.networkRows)
	//  sink table with vertex = f.networkRows*ff.networkCols+1
	for i := 0; i < ff.networkRows; i++ {
		v := ff.networkRows*ff.networkCols + 1
		w := i*ff.networkCols + ff.networkCols
		flow := ff.adj[v][w].flow
		color := "edgenoflow"
		if flow > 0 {
			color = "edgeflow"
		}
		ff.plot.SinkFlows[i] = Cell{Value: strconv.Itoa(ff.adj[v][w].flow),
			Name: fmt.Sprintf("%d-%d", v, w), BGColor: color}
	}

	// network table forward flow edges only (flow leaving the vertex)
	ff.plot.NetFlows = make([][]Cell, (ff.networkRows-1)*3+1)
	for row := 0; row < (ff.networkRows-1)*3+1; row++ {
		// one less edge than the number of vertices in each row
		ff.plot.NetFlows[row] = make([]Cell, ff.networkCols-1)
	}
	vStart := 1
	for col := 0; col < ff.networkCols-1; col++ {
		v := vStart
		for row := 0; row < (ff.networkRows-1)*3+1; {
			// flow edges for all network rows
			w := v + 1
			flow := ff.adj[v][w].flow
			color := "edgenoflow"
			if flow > 0 {
				color = "edgeflow"
			}
			ff.plot.NetFlows[row][col] = Cell{Value: strconv.Itoa(flow),
				Name: fmt.Sprintf("%d-%d", v, w), BGColor: color}
			// flow edges for the first network row
			if v == vStart {
				w := v + 1 + ff.networkCols
				flow := ff.adj[v][w].flow
				color := "edgenoflow"
				if flow > 0 {
					color = "edgeflow"
				}
				ff.plot.NetFlows[row+1][col] = Cell{Value: strconv.Itoa(ff.adj[v][w].flow),
					Name: fmt.Sprintf("%d-%d", v, w), BGColor: color}
				v += ff.networkCols
				row += 2
			} else if v == vStart+(ff.networkRows-1)*ff.networkCols { // flow edges for the last network row
				w := v + 1 - ff.networkCols
				flow := ff.adj[v][w].flow
				color := "edgenoflow"
				if flow > 0 {
					color = "edgeflow"
				}
				ff.plot.NetFlows[row+1][col] = Cell{Value: strconv.Itoa(ff.adj[v][w].flow),
					Name: fmt.Sprintf("%d-%d", v, w), BGColor: color}
				v += ff.networkCols
				row += 2
			} else { // all the other network rows in-between first and last rows
				w := v + 1 + ff.networkCols
				flow := ff.adj[v][w].flow
				color := "edgenoflow"
				if flow > 0 {
					color = "edgeflow"
				}
				ff.plot.NetFlows[row+1][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1+ff.networkCols].flow),
					Name: fmt.Sprintf("%d-%d", v, v+1+ff.networkCols), BGColor: color}
				w = v + 1 - ff.networkCols
				flow = ff.adj[v][w].flow
				color = "edgenoflow"
				if flow > 0 {
					color = "edgeflow"
				}
				ff.plot.NetFlows[row+2][col] = Cell{Value: strconv.Itoa(ff.adj[v][v+1-ff.networkCols].flow),
					Name: fmt.Sprintf("%d-%d", v, v+1-ff.networkCols), BGColor: color}
				v += ff.networkCols
				row += 3
			}
		}
		vStart++
	}

	// Display the network flow
	ff.plot.Flow = strconv.Itoa(ff.flow)
	return nil
}

// HTTP handler for /fordfulkersonflow connections
func handleFordFulkersonFlow(w http.ResponseWriter, r *http.Request) {

	// Create the Ford Fulkerson maxflow instance
	fordFulkerson := &FordFulkerson{}

	// Accumulate error
	status := make([]string, 0)

	// Get the previously calculated MaxFlow results from disk
	err := fordFulkerson.getMaxFlowResults()
	if err != nil {
		fmt.Printf("getMaxFlowResults error: %v\n", err)
		status = append(status, err.Error())
	}

	// Draw Flow Network into 300 x 300 cell 2px grid
	err = fordFulkerson.plotFlowNetwork()
	if err != nil {
		fmt.Printf("plotSP error: %v\n", err)
		status = append(status, err.Error())
	}

	// Show the Flows of each FlowEdge
	err = fordFulkerson.showFlows()
	if err != nil {
		fmt.Printf("showFlows error: %v\n", err)
		status = append(status, err.Error())
	}

	// Status
	if len(status) > 0 {
		fordFulkerson.plot.Status = strings.Join(status, ", ")
	} else {
		fordFulkerson.plot.Status = "Click 'Capacities' link to change capacities"
	}

	// Write to HTTP using template and grid
	if err := tmplFormFlows.Execute(w, fordFulkerson.plot); err != nil {
		log.Fatalf("Write to HTTP output using template with grid error: %v\n", err)
	}
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
