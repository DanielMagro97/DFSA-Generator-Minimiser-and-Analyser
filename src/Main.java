import java.util.*;
import java.util.concurrent.ThreadLocalRandom;  // for random number generation

public class Main {

    // Global Variable, boolean array storing whether each state is reachable or not
    private static boolean[] reachableState;

    public static void main(String[] args) {
        System.out.println("1. Generating the DFSA:");
        DFSA dfsaA = generateDFSA();
        // TODO debug step
        //printDFSA(dfsaA);

        System.out.println("\n2. Finding the Depth of A:");
        System.out.println("Number of States in A: " + dfsaA.states.size());
        System.out.println("Depth of A: " + findDepth(dfsaA));

        System.out.println("\n3. Minimising DFSA A:");
        DFSA dfsaM = minimiseDFSA(dfsaA);
        // TODO debug step
        //printDFSA(dfsaM);

        System.out.println("\n4. Finding the Depth of M:");
        System.out.println("Number of States in M: " + dfsaM.states.size());
        System.out.println("Depth of M: " + findDepth(dfsaM));

        System.out.println("\n5. Strongly Connected Components in DFSA M:");
        getSCCs(dfsaM);
        System.out.println("Number of Strongly Connected Components in M: " + SCCs.size());
        int maxSCC = SCCs.get(0).size();
        int minSCC = SCCs.get(0).size();
        for (int i = 1; i < SCCs.size(); i++) {
            if (SCCs.get(i).size() > maxSCC) {
                maxSCC = SCCs.get(i).size();
            }
            if (SCCs.get(i).size() < minSCC) {
                minSCC = SCCs.get(i).size();
            }
        }
        System.out.println("Number of States in the largest SCC: " + maxSCC);
        System.out.println("Number of States in the smallest SCC: " + minSCC);
        // TODO debug step
        //printSCCs(SCCs);

        System.out.println("\n6. Simple Cycles in DFSA M:");
        ArrayList<ArrayList<Integer>> simpleCycles = findSimpleCycles(dfsaM, SCCs);
        System.out.println("Number of Simple Cycles in M: " + simpleCycles.size());
        int maxCycle = simpleCycles.get(0).size();
        int minCycle = simpleCycles.get(0).size();
        for (int i = 1; i < simpleCycles.size(); i++) {
            if (simpleCycles.get(i).size() > maxCycle) {
                maxCycle = simpleCycles.get(i).size();
            }
            if (simpleCycles.get(i).size() < minCycle) {
                minCycle = simpleCycles.get(i).size();
            }
        }
        // deduct 1 since the start state is included twice in every cycle
        System.out.println("Number of States in the longest simple cycle: " + (maxCycle-1));
        System.out.println("Number of States in the shortest simple cycle: " + (minCycle-1));
        // TODO debug step
        //printCycles(simpleCycles);
    }




    private static DFSA generateDFSA(){
        // Setting the number of states of the DFSA to a random number between 16 and 64, both inclusive
        int n = ThreadLocalRandom.current().nextInt(16, 64 + 1);
        // TODO debug step
        //n = 6;

        DFSA dfsa = new DFSA(n);

        // For every state which needs to be generated:
        for (int i = 0; i < n; i++) {
            // initialise a new State object
            State state = new State();

            state.stateID = i;

            // Randomly decide whether the state should be an Accepting or Rejecting state
            state.finalState = ThreadLocalRandom.current().nextInt(0, 1+1) == 1;

            // Randomly decide to which state the current state will have an 'a' transition to
            dfsa.adjacencyMatrix[i][0] = ThreadLocalRandom.current().nextInt(0, n);

            // Randomly decide to which state the current state will have a 'b' transition to
            dfsa.adjacencyMatrix[i][1] = ThreadLocalRandom.current().nextInt(0, n);

            dfsa.states.add(state);
        }

        // setting the start state of the DFSA to a random state
        dfsa.startStateID = ThreadLocalRandom.current().nextInt(0, n);

        return dfsa;
    }




    private static int findDepth(DFSA dfsa){
        // stores the State IDs of the states in the current 'level' of the BFS
        ArrayList<Integer> currentLevelStates = new ArrayList<>();
        // add the start state of the DFSA to the current level
        currentLevelStates.add(dfsa.startStateID);

        // boolean array which marks whether a state has been visited or not
        boolean[] visited = new boolean[dfsa.states.size()];

        // initialise the depth of the DFSA to 0
        int depth = 0;

        // loop until the return statement is reached
        while (true) {
            // declare an arraylist of state IDs of states reachable by the current level
            ArrayList<Integer> nextLevelStates = new ArrayList<>();

            // mark all the states of the current level as visited
            for (int i = 0; i < currentLevelStates.size(); i++) {
                visited[currentLevelStates.get(i)] = true;
            }

            // find the states which can be reached from the states in the current level
            for (int i = 0; i < currentLevelStates.size(); i++) {
                // if the state reachable by the 'a' transition from a state has not already been visited
                if (!visited[ dfsa.adjacencyMatrix[currentLevelStates.get(i)][0] ]) {
                    // add the state ID of the state reachable by the 'a' transition to the nextLevelStates arraylist
                    nextLevelStates.add(dfsa.adjacencyMatrix[currentLevelStates.get(i)][0]);
                }
                // if the state reachable by the 'b' transition from a state has not already been visited
                if (!visited[ dfsa.adjacencyMatrix[currentLevelStates.get(i)][1] ]) {
                    // add the state ID of the state reachable by the 'b' transition to the nextLevelStates arraylist
                    nextLevelStates.add(dfsa.adjacencyMatrix[currentLevelStates.get(i)][1]);
                }
            }

            // copying nextLevelStates into currentLevelStates - for the next iteration
            currentLevelStates = new ArrayList<Integer>(nextLevelStates);

            // When all the 'next states' have already been visited, this means that either
            // i) all the states of the DFSA have been visited
            // ii) all the states which can be visited starting from the Start State have been visited (i.e. some states
            // cannot be reached)
            if (nextLevelStates.size() == 0) {
                // copy the contents of the visited array to the reachableState array.
                // any states which were not visited by the end of this algorithm must be unreachable
                // this is important for the minimisation of the DFSA
                reachableState = Arrays.copyOf(visited, visited.length);

                return depth;
            } else {
                depth++;
                continue;
            }
        }
    }




    private static DFSA minimiseDFSA(DFSA dfsaA) {
        // divide the list of states into accepting states and rejecting states
        Set<Integer> acceptingStates = new HashSet<>();
        Set<Integer> rejectingStates = new HashSet<>();

        for (State state : dfsaA.states) {
            // if the current state is reachable
            if (reachableState[state.stateID]) {
                // if the current state is an accepting state
                if (state.finalState) {
                    // add it to the partition (set) of accepting states
                    acceptingStates.add(state.stateID);
                } else {
                    // else, add it to the partition of rejecting states
                    rejectingStates.add(state.stateID);
                }
            }
        }

        // P := {F, Q \ F};
        Set< Set<Integer> > P = new HashSet<>();
        P.add(acceptingStates);
        P.add(rejectingStates);

        // W := {F};
        Set< Set<Integer> > W = new HashSet<>();
        W.add(acceptingStates);

        // while (W is not empty) do
        while (!W.isEmpty()) {
            // choose and remove a set A from W
            Set<Integer> A = W.iterator().next();
            W.remove(A);

            // for each c in Σ do
            for (int i = 0; i < dfsaA.adjacencyMatrix[0].length; i++) {

                // let X be the set of states for which a transition on c leads to a state in A
                Set<Integer> X = new HashSet<>();
                for (State state : dfsaA.states) {
                    if (reachableState[state.stateID]) {
                        if (A.contains(dfsaA.adjacencyMatrix[state.stateID][i])) {
                            X.add(state.stateID);
                        }
                    }
                }

                // for each set Y in P for which X ∩ Y is nonempty and Y \ X is nonempty do
                Set<Set<Integer>> Pcopy = new HashSet<>(P);
                for (Set<Integer> setInP : Pcopy) {
                    Set<Integer> Y = new HashSet<>(setInP);

                    Set<Integer> intersection = new HashSet<>(X);
                    intersection.retainAll(Y);

                    Set<Integer> difference = new HashSet<>(Y);
                    difference.removeAll(X);

                    if (!intersection.isEmpty() && !difference.isEmpty()) {
                        // replace Y in P by the two sets X ∩ Y and Y \ X
                        P.remove(Y);
                        P.add(intersection);
                        P.add(difference);

                        // if Y is in W
                        if (W.contains(Y)) {
                            // replace Y in W by the same two sets
                            W.remove(Y);
                            W.add(intersection);
                            W.add(difference);
                        } else {
                            // else
                            // if |X ∩ Y| <= |Y \ X|
                            if (intersection.size() <= difference.size()) {
                                // add X ∩ Y to W
                                W.add(intersection);
                            } else {
                                // else
                                // add Y \ X to W
                                W.add(difference);
                            }
                        }
                    }
                }
            }
        }

        // Declare DFSA M
        DFSA dfsaM = new DFSA(P.size());

        // Turn the set P into an ArrayList such that its elements now have a fixed order
        ArrayList<Set<Integer>> Parray = new ArrayList<>(P);

        // iterate through all of the partitions in the P arraylist
        int i = 0;
        for (Set<Integer> setInP : Parray) {
            // for each partition/element, declare a new state
            State state = new State();

            // assign an arbitrary ID to that state
            state.stateID = i;

            // if the partition contains the start state of DFSA A
            if (setInP.contains(dfsaA.startStateID)) {
                // set the state representing this partition as the start state of DFSA M
                dfsaM.startStateID = i;
            }

            // if the states in this partition are final states
            if (setInP.iterator().hasNext() && dfsaA.states.get(setInP.iterator().next()).finalState){
                // set the state representing this partition as an accepting state
                state.finalState = true;
            } else {
                // if not, set the state representing this partition as a rejecting state
                state.finalState = false;
            }

            // Check to which states the states in this partition lead to in DFSA A
            int oldA = -1, oldB = -1;
            for (int j = 0; j < dfsaA.states.size(); j++) {
                // when the state is found, store which states it used to have transitions to
                if (setInP.contains(dfsaA.states.get(j).stateID)) {
                    oldA = dfsaA.adjacencyMatrix[j][0];
                    oldB = dfsaA.adjacencyMatrix[j][1];
                    break;
                }
            }

            // when it is known which states a and b lead to in DFSA A
            // find which partition those states are in now
            // and from the current state draw a transition to the states representing the partitions in which a and b used to lead to respectively
            for (int j = 0; j < Parray.size(); j++) {
                if (Parray.get(j).contains(oldA)) {
                    dfsaM.adjacencyMatrix[i][0] = j;
                }

                if (Parray.get(j).contains(oldB)) {
                    dfsaM.adjacencyMatrix[i][1] = j;
                }
            }

            dfsaM.states.add(state);

            i++;
        }

        return dfsaM;
    }




    // Global variables used for Tarjan's Algorithm - SCC finder
    // ArrayList of SCCs. Each element in the ArrayList is an ArrayList of StateIDs (Integers).
    private static ArrayList<ArrayList<Integer>> SCCs;
    // stores the last used index (smallest unused index)
    private static int lastIndex;
    // array of integers which stores the index for each state
    private static int[] index;
    // array of integers which stores the lowlink for each state
    private static int[] lowlink;
    // stack storing the IDs of the states which make up the current SCC
    private static Stack<Integer> stateIDStack;
    // array of booleans storing whether a state is on the stack or not
    private static boolean[] onStack;

    private static ArrayList<ArrayList<Integer>> getSCCs(DFSA dfsa){
        // Initialising the global variables utilised by this method
        SCCs = new ArrayList<>();

        // set the smallest unused index to 0
        lastIndex = 0;

        index = new int[dfsa.states.size()];
        // fill the index array with -1, which represents 'undefined'
        Arrays.fill(index, -1);

        lowlink = new int[dfsa.states.size()];
        // fill the lowlink array with -1, which represents 'undefined'
        Arrays.fill(lowlink, -1);

        stateIDStack = new Stack<>();
        onStack = new boolean[dfsa.states.size()];

        // loops through every state in the DFSA
        for (int i = 0; i < dfsa.states.size(); i++) {
            // if the index for a node is undefined
            if (index[i] == -1) {
                // then call the strongConnect method and pass the current node as an argument
                strongConnect(dfsa.states.get(i).stateID, dfsa);
            }
        }

        return SCCs;
    }

    private static void strongConnect(int stateID, DFSA dfsa){
        // set the index and lowlink for the current state to the smallest unused index
        index[stateID] = lastIndex;
        lowlink[stateID] = lastIndex;

        // increment the smallest unused index
        lastIndex++;

        // add the current state to the stateStack
        stateIDStack.push(stateID);
        // and mark the current state as being onStack
        onStack[stateID] = true;

        // for all the successors of the current state
        // (the states reachable from the current state by an 'a' or a 'b' transition)
        for (int i = 0; i < dfsa.adjacencyMatrix[0].length; i++) {
            int successorID = dfsa.adjacencyMatrix[stateID][i];
            // if the index of the next state is still undefined, make a recursive call with that state
            if (index[successorID] == -1) {
                strongConnect(successorID, dfsa);
                // set the lowlink of the current state as either the lowlink of the next state or
                // the lowlink of the current state, whichever is smaller
                lowlink[stateID] = (lowlink[successorID] < lowlink[stateID]) ? lowlink[successorID] : lowlink[stateID];
            } else if (onStack[successorID]) {
                // if the next state is on the stack, then it is part of the current SCC
                // set the lowlink of the current node as either the lowlink of the current state or
                // the index of the next state, whichever is smaller
                lowlink[stateID] = (lowlink[stateID] < index[successorID]) ? lowlink[stateID] : index[successorID];
            }
        }

        // if the current node is a root node, then a SCC has been found
        if (lowlink[stateID] == index[stateID]) {
            ArrayList<Integer> SCC = new ArrayList<>();

            // pop states from the stack until the last popped state is the current state
            // add each popped state to the current SCC
            int w;
            do {
                w = stateIDStack.pop();
                onStack[w] = false;
                SCC.add(w);
            } while (w != stateID);

            // finally, return the newly created SCC
            SCCs.add(SCC);
        }
    }




    // Global Variables used for Johnson's Algorithm - Simple Cycles finder
    // ArrayList of ArrayLists of StateIDs which will store all the Simple Cycles in the DFSA
    // where each element of the ArrayList is a Simple Cycle (ArrayList)
    // and each element of the inner ArrayList (Cycle) is a State ID (Integer)
    private static ArrayList<ArrayList<Integer>> simpleCycles;
    // Stack storing the State IDs making up the current cycle
    private static Stack cycleStack;
    // Set of State IDs (Integers) of states which will not be visited
    private static Set<Integer> blockedSet;
    // Map of blocked State IDs (Integers) pointing to a Set of State IDs,
    // which will be unblocked too when the former is unblocked
    private static Map<Integer, Set<Integer>> blockedMap;

    private static ArrayList<ArrayList<Integer>> findSimpleCycles(DFSA dfsa, ArrayList<ArrayList<Integer>> SCCs){
        // Initialising the global variables utilised by this method
        simpleCycles = new ArrayList<>();
        cycleStack = new Stack();
        blockedSet = new HashSet<>();
        blockedMap = new HashMap<>();

        // Looping through all the SCC within the DFSA, since Simple Cycles are restricted within one SCC
        for (ArrayList<Integer> scc : SCCs) {
            // looping through every state in the SCC
            for (int i = 0; i < scc.size(); i++) {
                int startNode = scc.get(i);
                blockedSet.clear();
                blockedMap.clear();
                findSimpleCyclesInSCC(dfsa, scc, startNode, startNode);
            }
        }

        return simpleCycles;
    }

    private static boolean findSimpleCyclesInSCC(DFSA dfsa, ArrayList<Integer> scc, int startNode, int currentNode){
        boolean foundCycle = false;
        cycleStack.push(currentNode);
        blockedSet.add(currentNode);

        // exploring all neighbours
        for (int j = 0; j < dfsa.adjacencyMatrix[0].length; j++) {

            int nextNode = dfsa.adjacencyMatrix[currentNode][j];

            // cycle found
            if (nextNode == startNode) {
                // ArrayList of StateIDs (Integers) which stores the current cycle
                ArrayList<Integer> cycle = new ArrayList<>();
                // push the first state of the cycle to the top of the stack
                cycleStack.push(startNode);
                // add the contents of the stack to the array list of State IDs
                cycle.addAll(cycleStack);
                // pop the start node from the top of the stack
                cycleStack.pop();
                // add the cycle to the set of all the simpleCycles
                simpleCycles.add(cycle);
                // mark the foundCycle boolean as true
                foundCycle = true;
            } else {
                // if the next node is not in the blocked set
                if (!blockedSet.contains(nextNode)) {
                    foundCycle = foundCycle || findSimpleCyclesInSCC(dfsa, scc, startNode, nextNode);
                }
            }
        }

        // if a cycle was found with the current state
        if (foundCycle) {
            // then call the unblock method with the current state as an argument, which recursively unblocks the
            // node and all other nodes which are dependent on the current node (via blockedMap)
            unblock(currentNode);
        } else {
            // add all the states reachable by a transition from the current state to the blocked map
            Set<Integer> adjacentNodes = new HashSet<>();
            for (int i = 0; i < dfsa.adjacencyMatrix[0].length; i++) {
                adjacentNodes.add(dfsa.adjacencyMatrix[currentNode][i]);
            }
            blockedMap.put(currentNode, adjacentNodes);
        }

        cycleStack.pop();
        return foundCycle;
    }

    // Method which removes a given node from the blocked set
    // and removes any nodes which are being blocked by it (from blockedMap)
    private static void unblock(int node) {
        // remove the node from the blockedSet
        blockedSet.remove(node);
        // if the node being unblocked has a set of states waiting to be unblocked depending on the node
        if(blockedMap.get(node) != null) {
            // get the nodes dependent on the node being unblocked from the blockedMap
            Set<Integer> dependentUnblock = blockedMap.get(node);
            // recursively unblock each of those nodes if they are still in the blocked set
            for (int unblockNode : dependentUnblock) {
                if (blockedSet.contains(unblockNode)) {
                    unblock(unblockNode);
                }
            }
            // remove the entry for the node from the blocked map
            blockedMap.remove(node);
        }
    }






    // Methods used mostly for debugging purposes:
    private static void printDFSA(DFSA dfsa){
        System.out.println();
        System.out.println("Start State ID: " + dfsa.startStateID);
        for (int i = 0; i < dfsa.states.size(); i++) {
            if (dfsa.states.get(i).finalState) {
                System.out.println(i + " is accepting\t\t" + "a to " + dfsa.adjacencyMatrix[i][0] +
                        "\t\tb to " + dfsa.adjacencyMatrix[i][1]);
            } else {
                System.out.println(i + " is rejecting\t\t" + "a to " + dfsa.adjacencyMatrix[i][0] +
                        "\t\tb to " + dfsa.adjacencyMatrix[i][1]);
            }
        }
    }

    private static void printSCCs(ArrayList<ArrayList<Integer>> SCCs){
        for (int i = 0; i < SCCs.size(); i++) {
            for (int j = 0; j < SCCs.get(i).size(); j++) {
                System.out.print(SCCs.get(i).get(j) + "\t");
            }
            System.out.println();
        }
    }

    private static void printCycles(ArrayList<ArrayList<Integer>> cycles){
        for (ArrayList<Integer> cycle : cycles) {
            for (Integer node : cycle) {
                System.out.print(node + "\t");
            }
            System.out.println();
        }
    }
}
