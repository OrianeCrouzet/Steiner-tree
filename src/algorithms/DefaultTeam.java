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

public class DefaultTeam {
	final static int INF = 1000000000;
	  
	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
        int n = points.size();
        int m = hitPoints.size();

        // Remplace Floyd-Warshall par Dijkstra pour plus d'efficacité
        double[][] shortestPaths = computeAllPairShortestPaths(points, edgeThreshold);

        // Construire le graphe réduit pour hitPoints
        double[][] K = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                K[i][j] = shortestPaths[points.indexOf(hitPoints.get(i))][points.indexOf(hitPoints.get(j))];
            }
        }

        // Construire un MST sur ce graphe
        List<int[]> mstEdges = kruskalMST(K);

        // Expansion du MST en un arbre plus complet
        Set<Point> steinerPoints = new HashSet<>(hitPoints);
        Map<Point, List<Point>> adjacencyList = new HashMap<>();

        for (Point p : points) {
            adjacencyList.put(p, new ArrayList<>());
        }

        // Remplacement des arêtes MST par les plus courts chemins
        for (int[] edge : mstEdges) {
            Point u = hitPoints.get(edge[0]);
            Point v = hitPoints.get(edge[1]);
            List<Point> chemin = reconstructShortestPath(points, shortestPaths, u, v);

            // Ajouter le chemin à l'arbre
            for (int i = 0; i < chemin.size() - 1; i++) {
                Point a = chemin.get(i);
                Point b = chemin.get(i + 1);
                adjacencyList.get(a).add(b);
                adjacencyList.get(b).add(a);
                steinerPoints.add(a);
                steinerPoints.add(b);
            }
        }

        // Construire et retourner un Tree2D valide
        return buildTree(hitPoints.get(0), adjacencyList, new HashSet<>());
    }

    private double[][] computeAllPairShortestPaths(ArrayList<Point> points, int edgeThreshold) {
        int n = points.size();
        double[][] dist = new double[n][n];

        for (int i = 0; i < n; i++) {
            dist[i] = dijkstra(points, i, edgeThreshold);
        }
        return dist;
    }

    private double[] dijkstra(ArrayList<Point> points, int src, int edgeThreshold) {
        int n = points.size();
        double[] dist = new double[n];
        Arrays.fill(dist, INF);
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
                        pq.offer(new int[]{v, (int) dist[v]});
                    }
                }
            }
        }
        return dist;
    }

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

    private int find(int[] parent, int i) {
        if (parent[i] != i) parent[i] = find(parent, parent[i]);
        return parent[i];
    }

    private void union(int[] parent, int x, int y) {
        parent[find(parent, x)] = find(parent, y);
    }

    private List<Point> reconstructShortestPath(ArrayList<Point> points, double[][] shortestPaths, Point u, Point v) {
        List<Point> chemin = new ArrayList<>();
        chemin.add(u);
        chemin.add(v);
        return chemin;
    }

    private Tree2D buildTree(Point root, Map<Point, List<Point>> adjacencyList, Set<Point> visited) {
        if (visited.contains(root)) return null;
        visited.add(root);

        List<Tree2D> subtrees = new ArrayList<>();
        for (Point neighbor : adjacencyList.get(root)) {
            if (!visited.contains(neighbor)) {
                Tree2D subtree = buildTree(neighbor, adjacencyList, visited);
                if (subtree != null) {
                    subtrees.add(subtree);
                }
            }
        }
        return new Tree2D(root, new ArrayList<>(subtrees));
    }
    
    public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
	    return null;
	}
}
