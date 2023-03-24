# fordfulkersonmaxflow
Find the maximum flow in a st-flow network using the Ford-Fulkerson shortest-augmenting path maxflow algorithm.
This program is a web application that is written in go and uses the http/template package to create the HTML.  
The user connects to the Ford-Fulkerson MaxFlow server in a web browser at http://127.0.0.1:8080/fordfulkersongraphoptions.
The user enters the number of rows and columns in the st-flow network.  Two to ten rows or columns can be constructed in
addition to the source and sink vertices.  An n by m graph has n rows consisting of m vertices in each row for a total of
n x m vertices.  The network topology has edges connecting the row vertices to each other as well as the diagonally opposite
vertices.  The capacity of each edge is entered by the user as integers.  The source vertex is the origin of the flow, which
supplies the n x m network, which terminates at the sink vertex.  Clicking on the submit button will send the capacity 
selections to the server, which then calculates the maximum flow with the shortest-augmenting path.  A breadth-first search
is employed to find the shortest augmenting paths.  The resulting flow is displayed by clicking on the Flows link.  A plot of
the flow network is diplayed on the left, with green edges carrying flow and red edges carrying no flow.  The tables on the
right side show the amount of flow in each v-w edge, where v and w are the edge vertices.  By selecting the Capacities link,
the Capacity tables are displayed allowing you to change the edge capacities and calculate a new max flow.

![image](https://user-images.githubusercontent.com/117768679/227620601-0aa29927-ba10-401d-9370-0822f2f096ac.png)
![image](https://user-images.githubusercontent.com/117768679/227620363-8333452b-6642-4d57-b648-4bf64135622c.png)
![image](https://user-images.githubusercontent.com/117768679/227643076-7033224b-0159-4a15-9e2e-72bb42f29479.png)
