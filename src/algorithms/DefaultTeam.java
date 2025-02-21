package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.List;
import java.util.Set;
import java.util.Map;
import java.util.HashSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Collections;

public class DefaultTeam {
	final static int INF = 1000000000;
	final static int B = 1664;
	  
	// ******************** STEINER sans budget ********************
	
	/*
	 * 
	 */
	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		int n = points.size();
        int m = hitPoints.size();
        int [][] predecessors = new int[n][n];

        double[][] shortestPaths = computeAllPairShortestPaths(points, edgeThreshold, predecessors);

        double[][] K = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                K[i][j] = shortestPaths[points.indexOf(hitPoints.get(i))][points.indexOf(hitPoints.get(j))];
            }
        }

        List<int[]> mstEdges = kruskalMST(K);

        Set<Point> steinerPoints = new HashSet<>(hitPoints);
        Map<Point, List<Point>> adjacencyList = new HashMap<>();

        for (Point p : points) {
            adjacencyList.put(p, new ArrayList<>());
        }

        for (int[] edge : mstEdges) {
            Point u = hitPoints.get(edge[0]);
            Point v = hitPoints.get(edge[1]);
            List<Point> chemin = reconstructShortestPath(points, predecessors, u, v);

            for (int i = 0; i < chemin.size() - 1; i++) {
                Point a = chemin.get(i);
                Point b = chemin.get(i + 1);
                adjacencyList.get(a).add(b);
                adjacencyList.get(b).add(a);
                steinerPoints.add(a);
                steinerPoints.add(b);
            }
        }

        return buildSteinerTree(hitPoints.get(0), adjacencyList, new HashSet<>());
    }
	
	/*
	 * 
	 */
	private double[][] computeAllPairShortestPaths(ArrayList<Point> points, int edgeThreshold, int[][] predecessors) {
		int n = points.size();
	    double[][] dist = new double[n][n];

	    for (int i = 0; i < n; i++) {
	        int[] predecessor = new int[n];
	        dist[i] = dijkstra(points, i, edgeThreshold, predecessor);
	        predecessors[i] = predecessor;
	    }
	    return dist;
	}

	/*
	 * 
	 */
    private double[] dijkstra(ArrayList<Point> points, int src, int edgeThreshold, int[] predecessor) {
    	int n = points.size();
        double[] dist = new double[n];
        Arrays.fill(dist, INF);
        Arrays.fill(predecessor, -1); // -1 = "no predecessor"
        dist[src] = 0;

        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingDouble(a -> a[1]));
        pq.offer(new int[]{src, 0});

        while (!pq.isEmpty()) {
            int[] node = pq.poll();
            int u = node[0];
            double d = node[1];

            if (d > dist[u]) continue;

            for (int v = 0; v < n; v++) {
                if (u != v && points.get(u).distance(points.get(v)) <= edgeThreshold) {
                    double weight = points.get(u).distance(points.get(v));
                    if (dist[u] + weight < dist[v]) {
                        dist[v] = dist[u] + weight;
                        predecessor[v] = u; 
                        pq.offer(new int[]{v, (int) dist[v]});
                    }
                }
            }
        }
        return dist;
    }

    /*
     * 
     */
    private List<int[]> kruskalMST(double[][] graph) {
        int n = graph.length;
        List<int[]> edges = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                edges.add(new int[]{i, j, (int) graph[i][j]});
            }
        }
        edges.sort(Comparator.comparingInt(e -> e[2]));

        int[] parent = new int[n];
        for (int i = 0; i < n; i++) parent[i] = i;

        List<int[]> mst = new ArrayList<>();
        for (int[] edge : edges) {
            int u = edge[0], v = edge[1];
            if (find(parent, u) != find(parent, v)) {
                mst.add(new int[]{u, v});
                union(parent, u, v);
            }
        }
        return mst;
    }

    /*
     * 
     */
    private int find(int[] parent, int i) {
        if (parent[i] != i) parent[i] = find(parent, parent[i]);
        return parent[i];
    }

    /*
     * 
     */
    private void union(int[] parent, int x, int y) {
        parent[find(parent, x)] = find(parent, y);
    }

    /*
     * 
     */
    private List<Point> reconstructShortestPath(ArrayList<Point> points, int[][] predecessors, Point u, Point v) {
    	List<Point> path = new ArrayList<>();
        int start = points.indexOf(u);
        int end = points.indexOf(v); 

        if (predecessors[start][end] == -1) return path;

        while (end != start) {
            path.add(points.get(end));
            end = predecessors[start][end];
        }
        path.add(points.get(start));

        Collections.reverse(path);
        return path;
    }

    /*
     * 
     */
    private Tree2D buildSteinerTree(Point root, Map<Point, List<Point>> adjacencyList, Set<Point> visited) {
    	if (visited.contains(root)) return null;
        visited.add(root);

        List<Tree2D> subtrees = new ArrayList<>();
        for (Point neighbor : adjacencyList.get(root)) {
            if (!visited.contains(neighbor)) {
                Tree2D subtree = buildSteinerTree(neighbor, adjacencyList, visited);
                if (subtree != null) {
                    subtrees.add(subtree);
                }
            }
        }
        return new Tree2D(root, new ArrayList<>(subtrees));
    }
    
    
    // ******************** STEINER avec budget ********************
    
    /*
     * 
     */
    public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
        System.out.println("Debut budget");
        System.out.println("Nombre total de points : " + points.size());
        System.out.println("Nombre de hitPoints : " + hitPoints.size());

    	int n = points.size();
        Map<Point, Integer> indexMap = new HashMap<>();
        for (int i = 0; i < n; i++) indexMap.put(points.get(i), i);

        double[][] distances = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) distances[i][j] = 0;
                else {
                    double d = points.get(i).distance(points.get(j));
                    distances[i][j] = (d <= edgeThreshold) ? d : INF;
                }
            }
        }
        
        int m = hitPoints.size();
        double[][] shortestPaths = new double[m][m];
        int[][] predecessors = new int[m][m];
        for (int i = 0; i < m; i++) {
        	System.out.println("i : " + i);
        	Point hit = hitPoints.get(i);
            int src = points.indexOf(hit);
            System.out.println("Src : " + src);
            shortestPaths[i] = dijkstraWithEdges(points, src, edgeThreshold, predecessors[i]);
        }
        System.out.println("Sortie de boucle pour dijkstra");

        return constructBudgetSteinerTree(hitPoints, shortestPaths, predecessors, edgeThreshold);
    }
    
    /*
     * 
     */
    private double[] dijkstraWithEdges(ArrayList<Point> points, int src, int edgeThreshold, int[] predecessor) {
    	System.out.println("_______________");
    	System.out.println("Debut dijkstra");
    	int n = points.size();
        double[] dist = new double[n];
        Arrays.fill(dist, Double.POSITIVE_INFINITY);
        dist[src] = 0;
        predecessor[src] = -1; // -1 = "no predecessor"

        PriorityQueue<Edge> pq = new PriorityQueue<>(Comparator.comparingDouble(e -> e.weight));
        pq.add(new Edge(src, src, 0));

        while (!pq.isEmpty()) {
            Edge e = pq.poll();
            int u = e.v;
            if (e.weight > dist[u]) continue;

            for (int v = 0; v < n; v++) {
                if (u != v && points.get(u).distance(points.get(v)) <= edgeThreshold) {
                    double newDist = dist[u] + points.get(u).distance(points.get(v));
                    if (newDist < dist[v]) {
                        dist[v] = newDist;
                        predecessor[v] = u;
                        pq.add(new Edge(u, v, newDist));
                    }
                }
            }
        }
        System.out.println("Fin dijkstra");
        return dist;
    }
    
    /*
     * 
     */
    private Tree2D constructBudgetSteinerTree(ArrayList<Point> hitPoints, double[][] shortestPaths, int[][] predecessors, int budget) {
    	System.out.println("Debut build tree");
    	Set<Integer> connected = new HashSet<>();
        PriorityQueue<Edge> pq = new PriorityQueue<>(Comparator.comparingDouble(e -> e.weight));

        int first = 0;
        connected.add(first);

        for (int v = 1; v < hitPoints.size(); v++) {
            pq.add(new Edge(first, v, shortestPaths[first][v]));
        }

        Tree2D steinerTree = new Tree2D(hitPoints.get(first), new ArrayList<>());
        double totalCost = 0;

        while (!pq.isEmpty() && connected.size() < hitPoints.size()) {
            Edge e = pq.poll();
            if (totalCost + e.weight > budget) continue;

            connected.add(e.v);
            totalCost += e.weight;

            System.out.println("Début reconstruct path");
            List<Point> path = reconstructShortestPath(hitPoints, predecessors, hitPoints.get(e.u), hitPoints.get(e.v));
            System.out.println("Fin reconstruct path");
            Tree2D subtree = new Tree2D(path.get(0), new ArrayList<>());
            Tree2D current = subtree;
            for (int i = 1; i < path.size(); i++) {
                Tree2D next = new Tree2D(path.get(i), new ArrayList<>());
                current.getSubTrees().add(next);
                current = next;
            }
            steinerTree.getSubTrees().add(subtree);

            for (int v = 0; v < hitPoints.size(); v++) {
                if (!connected.contains(v)) {
                    pq.add(new Edge(e.v, v, shortestPaths[e.v][v]));
                }
            }
        }

        System.out.println("Fin build tree");
        return steinerTree;
    }

}
