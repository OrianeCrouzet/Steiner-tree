package algorithms;

/**
 * Cette classe représente une arête dans un graphe.
 * Une arête est définie par deux points (u et v) et un poids (weight).
 */
public class Edge {
    int u; 
    int v; 
    double weight;

    /**
     * Constructeur pour créer une nouvelle arête.
     *
     * @param u Le premier point de l'arête.
     * @param v Le deuxième point de l'arête.
     * @param weight Le poids de l'arête = distance euclidienne entre u et v
     */
    public Edge(int u, int v, double weight) {
        this.u = u;
        this.v = v;
        this.weight = weight;
    }
}