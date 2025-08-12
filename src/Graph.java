public class Graph {
    int first;
    int second;
    public Graph(int first, int second) {
        this.first = first;
        this.second = second;
    }
    @Override
    public String toString() {
        return "Graph{" +
                "first=" + first +
                ", second=" + second +
                '}';
    }
}
