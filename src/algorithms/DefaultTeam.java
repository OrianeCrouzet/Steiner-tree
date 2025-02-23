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

/**
 * Cette classe implémente des algorithmes pour calculer un arbre de Steiner avec et sans contrainte de budget.
 * L'heuristique se base sur celle présentée en TME 5, exercice 2.
 */
public class DefaultTeam {
	
	// ******************** Constantes ********************
	
    final static int INF = 1000000000; // Valeur représentant l'infini pour les distances
    final static int B = 1664; // Budget maximal pour l'arbre de Steiner avec budget

    
    // ******************** STEINER sans budget ********************

    /**
     * Calcule un arbre de Steiner sans contrainte de budget.
     *
     * @param points : La liste de tous les points disponibles.
     * @param edgeThreshold : La distance maximale pour qu'une arête soit valide entre deux points.
     * @param hitPoints : La liste des points à connecter (points de Steiner).
     * @return Un arbre de Steiner sous forme de Tree2D.
     */
    public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
        int n = points.size();
        int m = hitPoints.size();
        int[][] predecessors = new int[n][n];

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

    /**
     * Construit un arbre de Steiner à partir d'une liste d'adjacence.
     *
     * @param root : Le point racine de l'arbre.
     * @param adjacencyList : La liste d'adjacence représentant les connexions entre les points.
     * @param visited : Un ensemble de points déjà visités.
     * @return Un arbre de Steiner sous forme de Tree2D.
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

    /**
     * Calcule un arbre de Steiner avec une contrainte de budget.
     *
     * @param points : La liste de tous les points disponibles.
     * @param edgeThreshold : La distance maximale pour qu'une arête soit valide entre deux points.
     * @param hitPoints : La liste des points à connecter (points de Steiner).
     * @return Un arbre de Steiner sous forme de Tree2D.
     */
    public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
        int n = points.size();
        int m = hitPoints.size();
        int[][] predecessors = new int[n][n];

        double[][] shortestPaths = computeAllPairShortestPaths(points, edgeThreshold, predecessors);

        return constructBudgetSteinerTree(points, hitPoints, shortestPaths, predecessors, B);
    }

    /**
     * Construit un arbre de Steiner avec une contrainte de budget.
     *
     * @param points : La liste de tous les points disponibles.
     * @param hitPoints : La liste des points à connecter (points de Steiner).
     * @param shortestPaths : La matrice des plus courts chemins entre les points.
     * @param predecessors : La matrice des prédécesseurs.
     * @param budget : Le budget maximal pour la construction de l'arbre.
     * @return Un arbre de Steiner sous forme de Tree2D.
     */
    private Tree2D constructBudgetSteinerTree(ArrayList<Point> points, ArrayList<Point> hitPoints, double[][] shortestPaths, int[][] predecessors, int budget) {
        int m = hitPoints.size();
        Set<Integer> connected = new HashSet<>();
        Map<Point, List<Point>> adjacencyList = new HashMap<>();

        for (Point p : points) {
            adjacencyList.put(p, new ArrayList<>());
        }

        int maisonMereIndex = 0;
        connected.add(maisonMereIndex);

        PriorityQueue<Edge> pq = new PriorityQueue<>((edge1, edge2) -> {
            double ratio1 = edge1.weight / countNewPoints(edge1, connected);
            double ratio2 = edge2.weight / countNewPoints(edge2, connected);
            return Double.compare(ratio1, ratio2);
        });

        for (int i = 1; i < m; i++) {
            double weight = shortestPaths[points.indexOf(hitPoints.get(maisonMereIndex))][points.indexOf(hitPoints.get(i))];
            pq.add(new Edge(maisonMereIndex, i, weight));
        }

        double totalCost = 0;

        while (!pq.isEmpty() && totalCost < budget) {
            Edge e = pq.poll();

            if (totalCost + e.weight > budget) {
                continue;
            }

            if (connected.contains(e.v)) {
                continue;
            }

            connected.add(e.v);
            totalCost += e.weight;

            Point u = hitPoints.get(e.u);
            Point v = hitPoints.get(e.v);
            List<Point> path = reconstructShortestPath(points, predecessors, u, v);

            for (int i = 0; i < path.size() - 1; i++) {
                Point a = path.get(i);
                Point b = path.get(i + 1);
                adjacencyList.get(a).add(b);
                adjacencyList.get(b).add(a);
            }

            for (int i = 0; i < m; i++) {
                if (!connected.contains(i)) {
                    double weight = shortestPaths[points.indexOf(v)][points.indexOf(hitPoints.get(i))];
                    pq.add(new Edge(e.v, i, weight));
                }
            }
        }

        return buildSteinerTree(hitPoints.get(maisonMereIndex), adjacencyList, new HashSet<>());
    }
    
    // ******************** Fonctions utilitaires ********************
    
    /**
     * Reconstruit le plus court chemin entre deux points en utilisant la matrice des prédécesseurs.
     *
     * @param points : La liste de tous les points disponibles.
     * @param predecessors : La matrice des prédécesseurs.
     * @param u : Le point de départ.
     * @param v : Le point d'arrivée.
     * @return Une liste de points représentant le plus court chemin entre u et v.
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
    
    /**
     * Calcule les plus courts chemins entre toutes les paires de points en utilisant l'algorithme de Dijkstra.
     *
     * @param points : La liste de tous les points disponibles.
     * @param edgeThreshold : La distance maximale pour qu'une arête soit valide entre deux points.
     * @param predecessors : Une matrice pour stocker les prédécesseurs dans les plus courts chemins.
     * @return Une matrice des distances entre toutes les paires de points.
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

    /**
     * Algorithme de Dijkstra pour calculer les plus courts chemins à partir d'un point donné.
     *
     * @param points : La liste de tous les points disponibles.
     * @param src : L'indice du point donné.
     * @param edgeThreshold : La distance maximale pour qu'une arête soit valide entre deux points.
     * @param predecessor : Un tableau pour stocker les prédécesseurs dans le plus court chemin.
     * @return Un tableau des distances du point donné à tous les autres points.
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

    /**
     * Calcule un arbre couvrant minimal (MST) en utilisant l'algorithme de Kruskal.
     *
     * @param graph : Une matrice représentant les distances entre les points.
     * @return Une liste d'arêtes représentant l'arbre couvrant minimal.
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

    /**
     * Fonction qui trouve la racine d'un ensemble.
     *
     * @param parent : Le tableau des parents.
     * @param i : L'élément dont on cherche la racine.
     * @return La racine de l'élément donné.
     */
    private int find(int[] parent, int i) {
        if (parent[i] != i) parent[i] = find(parent, parent[i]);
        return parent[i];
    }

    /**
     * Fonction qui fait l'union de deux ensembles.
     *
     * @param parent : Le tableau des parents.
     * @param x : Le premier élément à fusionner.
     * @param y : Le deuxième élément à fusionner.
     */
    private void union(int[] parent, int x, int y) {
        parent[find(parent, x)] = find(parent, y);
    }
    
    /**
     * Compte le nombre de nouveaux points connectés par une arête.
     *
     * @param e : L'arête à évaluer.
     * @param connected : L'ensemble des points déjà connectés.
     * @return Le nombre de nouveaux points connectés par cette arête.
     */
    private int countNewPoints(Edge e, Set<Integer> connected) {
        return connected.contains(e.v) ? 0 : 1;
    }
    
}