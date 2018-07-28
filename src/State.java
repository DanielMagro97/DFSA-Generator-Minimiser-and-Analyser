public class State {
    // The ID of the current state. IDs range from 0 to n-1; where n is the number of states in the DFSA
    int stateID;

    // boolean value storing whether the state is final (accepting) or not (rejecting)
    // 1 = accepting, 0 = rejecting
    boolean finalState;
}
