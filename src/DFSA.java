import java.util.ArrayList;

class DFSA {
    // Storing the ID of the start state
    int startStateID;

    // Storing the set of states
    ArrayList<State> states;

    // A variation of an Adjacency Matrix was used so that it can be checked to which states the current state has
    // 'a' and 'b' transitions to in constant time
    // 2 rows and n columns. n being the number of states in the DFSA
    // Row 0 shows the outgoing 'a' transitions and Row 1 shows the outgoing 'b' transitions
    // For example: To find out which state can be reached by following an 'a' from state 14, check: adjacencyMatrix[0][14]
    int[][] adjacencyMatrix;

    // Constructor for a new DFSA, takes the number of states as a parameter
    DFSA(int numberOfStates){
        this.states = new ArrayList();
        this.adjacencyMatrix = new int[numberOfStates][2];
    }
}
