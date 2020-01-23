package main

import (
	"fmt"
	"math"
)

type Matrix [][]float64

type Tree []*Node

type Node struct {
	age                                                     float64
	label, event                                            string
	score                                                   map[int]float64 //Minimum Parsimony building gene tree
	rank                                                    int             //Species Star
	parent, child1, child2                                  *Node
	L                                                       []*Node // Possible to have one gene map to multiple species?
	cost, speciation, duplication, transfer, in, inAlt, out map[*Node]float64
}

func main() {
	fmt.Println("Reconciliation using U-MPR")

	var speciesTree Tree
	speciesTree = make([]*Node, 7)
	var v0, v1, v2, v3, v4, v5, v6 Node
	v0 = Node{parent: &v4, label: "A"}
	v1 = Node{parent: &v4, label: "B"}
	v2 = Node{parent: &v5, label: "C"}
	v3 = Node{parent: &v6, label: "D"}
	v4 = Node{parent: &v5, child1: &v0, child2: &v1, label: "Ancestor Species 1"}
	v5 = Node{parent: &v6, child1: &v4, child2: &v2, label: "Ancestor Species 2"}
	v6 = Node{child1: &v5, child2: &v3, label: "Ancestor Species 3"}
	speciesTree[0] = &v0
	speciesTree[1] = &v1
	speciesTree[2] = &v2
	speciesTree[3] = &v3
	speciesTree[4] = &v4
	speciesTree[5] = &v5
	speciesTree[6] = &v6

	var geneTree Tree
	geneTree = make([]*Node, 7)
	var gv0, gv1, gv2, gv3, gv4, gv5, gv6 Node
	gv0 = Node{parent: &gv4, label: "A"}
	gv1 = Node{parent: &gv4, label: "C"}
	gv2 = Node{parent: &gv5, label: "B"}
	gv3 = Node{parent: &gv6, label: "D"}
	gv4 = Node{parent: &gv5, child1: &gv0, child2: &gv1, label: "Ancestor Species 1"}
	gv5 = Node{parent: &gv6, child1: &gv4, child2: &gv2, label: "Ancestor Species 2"}
	gv6 = Node{child1: &gv5, child2: &gv3, label: "Ancestor Species 3"}
	geneTree[0] = &gv0
	geneTree[1] = &gv1
	geneTree[2] = &gv2
	geneTree[3] = &gv3
	geneTree[4] = &gv4
	geneTree[5] = &gv5
	geneTree[6] = &gv6

	//var Lmap map[]
	PLoss := 3
	Pduplication := 2
	Ptransfer := 1

	cost := UMPR(geneTree, speciesTree, len(geneTree), len(speciesTree), PLoss, Pduplication, Ptransfer)
	fmt.Println("Optimal cost is: ", cost)
}

// UMPR takes in both gene tree and species tree, the length of gene tree and species tree, and
// the cost of loss, duplication, and Ptransfer, and returns the minimum cost after inferring
// evolutionary events for every internal nodes in the gene tree.
func UMPR(geneT, speciesT Tree, geneN, speciesN, PLoss, Pduplication, Ptransfer int) float64 {
	//Initiailize tree
	InitializeTree(geneT, speciesN)

	geneLeave := (geneN + 1) / 2
	speciesLeave := (speciesN + 1) / 2

	// Initialize leaves mapping
	InitializeLMap(geneT, speciesT, geneLeave, speciesLeave)
	// Initialization of costs to infinity
	for _, g := range geneT {
		for _, s := range speciesT {
			InitializeCosts(g, s)
		}
	}

	// Get Le(G); leaveNode is an array of *Node only including leaves
	leaveNodesG := GetLeaveNodes(geneT, geneLeave)
	leaveNodesS := GetLeaveNodes(speciesT, speciesLeave)
	// Assign costs to the leave nodes of the gene tree
	for _, l := range leaveNodesG {
		for _, k := range leaveNodesS {
			// Mapped means k in L(g)
			if Mapped(l, k) {
				l.cost[k] = 0
			}
		}

		// For each s>=S L(g), intilize in(g,s) and inAlt(g,s)
		ancestors := FindAncestors(l, speciesT)

		for dist, a := range ancestors {
			// Distance is the edge between nodes
			// InitializeIn takes in a, an ancestor node in species tree, and l,
			// a leave node in gene tree
			l.in[a] = float64(PLoss * dist)

			l.inAlt[a] = 0.0
		}
		specieNodes := GetLeaveNodes(speciesT, speciesLeave)
		for _, i := range specieNodes {
			if i.label != l.label {
				l.out[i] = 0
			}
		}
	}

	// Get internal nodes of gene tree
	internalNodes := GetInternalNodes(geneT, geneLeave)

	// Starts from the root
	g := geneT[len(geneT)-1]
	s := speciesT[len(speciesT)-1]

	PostOrderGeneT(internalNodes, speciesT, geneT, g, s, speciesLeave, Pduplication, Ptransfer, PLoss)

	rootNode := geneT[len(geneT)-1]
	cost := OptimalCost(rootNode, speciesT)

	PrintTree(geneT)

	return cost
}

// InitializeTree takes in a gene tree and number of species tree to initialize
// the cost maps that we will use in our dynamic programming to get the minimum cost.
func InitializeTree(geneT Tree, speciesN int) {
	for _, i := range geneT {
		i.cost = make(map[*Node]float64, speciesN)
		i.speciation = make(map[*Node]float64, speciesN)
		i.duplication = make(map[*Node]float64, speciesN)
		i.transfer = make(map[*Node]float64, speciesN)
		i.in = make(map[*Node]float64, speciesN)
		i.inAlt = make(map[*Node]float64, speciesN)
		i.out = make(map[*Node]float64, speciesN)
	}
}

// OptimalCost takes in the root node of the gene tree and it returns
// the minimum cost mapping to a node on species tree
func OptimalCost(rootNode *Node, speciesT Tree) float64 {

	min := rootNode.cost[speciesT[0]]
	for i := 1; i < len(speciesT); i++ {
		if rootNode.cost[speciesT[i]] < min {
			min = rootNode.cost[speciesT[i]]
		}
	}
	return min
}

// PostOrderGeneT takes in the internal nodes of gene tree, the secies tree and the gene Tree
// it also takes in the root nodes for both gene tree and species tree
// This recursive function go through the gene tree internal nodes in a postorder and assign costs map value
func PostOrderGeneT(internalNodes, speciesT, geneT Tree, g, s *Node, speciesLeave, Pduplication, Ptransfer, PLoss int) {
	if g.child1 == nil || g.child2 == nil {
		return
	}
	PostOrderGeneT(internalNodes, speciesT, geneT, g.child1, s, speciesLeave, Pduplication, Ptransfer, PLoss)
	PostOrderGeneT(internalNodes, speciesT, geneT, g.child2, s, speciesLeave, Pduplication, Ptransfer, PLoss)

	PostOrderSpT(speciesT, geneT, g, s, speciesLeave, Pduplication, Ptransfer, PLoss)

	rootNode := speciesT[len(speciesT)-1]
	PreOrderIntNSpT(speciesT, rootNode, g)
}

// PreOrderIntNSpT takes in the species tree, the root node of the species tree,
// and one internal
func PreOrderIntNSpT(speciesT Tree, rootNode, g *Node) {
	if rootNode == nil {
		panic("Species tree do not have internal nodes!")
	}

	if rootNode.child1 == nil || rootNode.child2 == nil {
		return
	}

	s1 := rootNode.child1
	s2 := rootNode.child2
	g.out[s1] = Min2(g.out[rootNode], g.inAlt[s2])
	g.out[s2] = Min2(g.out[rootNode], g.inAlt[s1])

	PreOrderIntNSpT(speciesT, rootNode.child1, g)
	PreOrderIntNSpT(speciesT, rootNode.child2, g)
}

// PostOrderSpT takes in a speicies and a gene tree. It also takes in the rootnodes, and the costs are also passed down.
// It traverse the species tree in post order and assign cost map values for gene tree internal nodes
func PostOrderSpT(speciesT, geneT Tree, g, s *Node, speciesLeave, Pduplication, Ptransfer, PLoss int) {
	//fmt.Println("LLLLLL, ", g.label)
	if s == nil {
		return
	}
	PostOrderSpT(speciesT, geneT, g, s.child1, speciesLeave, Pduplication, Ptransfer, PLoss)
	PostOrderSpT(speciesT, geneT, g, s.child2, speciesLeave, Pduplication, Ptransfer, PLoss)

	g1 := g.child1
	g2 := g.child2

	if SLeave(s, speciesT, speciesLeave) {

		g.speciation[s] = float64(math.Inf(+1))
		g.duplication[s] = float64(Pduplication) + g1.cost[s] + g2.cost[s]
		//fmt.Println("sadfsafsadfsadf ", g1.cost[s], " ", g2.cost[s])

		if NotRoot(s, speciesT) {
			g.transfer[s] = float64(Ptransfer) + Min2(g1.in[s]+g2.out[s], g2.in[s]+g1.out[s])
		}
		//	fmt.Println("THE NUMBER: ", g.speciation[s], " ", g.duplication[s], " ", g.transfer[s])
		g.cost[s], g.event = Min3(g.speciation[s], g.duplication[s], g.transfer[s])
		//	fmt.Println("HY lavbel: ", g.label, "MAP TO: ", s.label, " event: ", g.event)

		g.in[s] = g.cost[s]
		g.inAlt[s] = g.cost[s]
	} else {
		s1 := s.child1
		s2 := s.child2

		g.speciation[s] = Min2(g1.in[s1]+g2.in[s2], g2.in[s1]+g1.in[s2])

		g.duplication[s] = float64(Pduplication) + MinDup(g1, g2, s, s1, s2, PLoss)
		if NotRoot(s, speciesT) {

			g.transfer[s] = float64(Ptransfer) + Min2(g1.in[s]+g2.out[s], g2.in[s]+g1.out[s])
		}
		g.cost[s], g.event = Min3(g.speciation[s], g.duplication[s], g.transfer[s])

		g.in[s] = Min3M(g.cost[s], g.in[s1]+float64(PLoss), g.in[s2]+float64(PLoss))
		g.inAlt[s] = Min3M(g.cost[s], g.inAlt[s1], g.inAlt[s2])
	}

}

// MinDup takes in g1, g2, s, s1, s2 and cost for loss and returns the minimum
// value of nine equalities
func MinDup(g1, g2, s, s1, s2 *Node, PLoss int) float64 {
	return Min9(g1.cost[s]+g2.in[s2]+float64(PLoss), g1.cost[s]+g2.in[s1]+float64(PLoss),
		g2.cost[s]+g1.in[s2]+float64(PLoss), g2.cost[s]+g1.in[s1]+float64(PLoss),
		g1.cost[s]+g2.cost[s], g1.in[s1]+g2.in[s2]+float64(2*PLoss),
		g1.in[s2]+g2.in[s1]+float64(2*PLoss), g1.in[s1]+g2.in[s1]+float64(2*PLoss),
		g1.in[s2]+g2.in[s2]+float64(2*PLoss))
}

// Min9 takes in 9 values and returns the minimum one
func Min9(a, b, c, d, e, f, g, h, i float64) float64 {
	arr := make([]float64, 9)
	arr[0] = a
	arr[1] = b
	arr[2] = c
	arr[3] = d
	arr[4] = e
	arr[5] = f
	arr[6] = g
	arr[7] = h
	arr[8] = i
	return Min(arr)
}

// Min takes in an array and returns the minimum value in this array
func Min(arr []float64) float64 {
	if len(arr) < 1 {
		panic("no items in array")
	}
	if len(arr) == 1 {
		return arr[0]
	}
	if arr[0] < arr[1] {
		return arr[0]
	}
	arr = arr[1:]
	return Min(arr)
}

// Min3 takes in three variables and returns the minimum value, and also the event
// with the minimum value
func Min3(a, b, c float64) (float64, string) {
	events := "speciation"
	min := a
	if a > b {
		min = b
		events = "duplication"
	}
	if min > c {
		events = "transfer"
		return c, events
	}

	return min, events
}

// Min3M takes in three variables and returns the minimum value
func Min3M(a, b, c float64) float64 {
	min := a
	if a > b {
		min = b
	}
	if min > c {
		return c
	}
	return min
}

// Min2 takes in two variable and returns the minimum value
func Min2(a, b float64) float64 {
	if a > b {
		return b
	}
	return a

}

// NotRoot takes in a node and a tree and checks if the node is the root of the
// tree. If it is not, return true.
func NotRoot(s *Node, speciesT Tree) bool {
	if s.label != speciesT[len(speciesT)-1].label {
		return true
	}
	return false
}

// SLeave takes in a node and a tree. It returns tree if the input node is
// a leave node
func SLeave(s *Node, speciesT Tree, speciesLeave int) bool {
	for i := 0; i < speciesLeave; i++ {
		if s.label == speciesT[i].label {
			return true
		}
	}
	return false
}

// GetInternalNodes takes in a tree and the leave number, it returns a tree
// only including the internal nodes of the tree
func GetInternalNodes(geneT Tree, geneLeave int) Tree {
	internals := make(Tree, 0)
	for i := geneLeave; i < len(geneT); i++ {
		internals = append(internals, geneT[i])
	}
	return internals
}

// FindAncestors takes in a gene tree leave g, and find find all the
// ancestors of l(g), which is all the nodes on the path between root(S) and L(g)
// in the species tree
func FindAncestors(l *Node, speciesT Tree) []*Node {
	// the mapped nodes in species tree
	sNodes := l.L
	ancestors := make([]*Node, 0)
	ancestors = append(ancestors, l.L[0])
	// Species Tree have parents by STAR alg
	for _, node := range sNodes {
		p := node

		for p.parent != nil {
			p = p.parent
			ancestors = append(ancestors, p)
		}
	}

	return ancestors
}

// Mapped takes in two nodes. It returns true if two nodes have the same labels
func Mapped(l, k *Node) bool {
	if l.label == k.label {
		return true
	}
	return false
}

// InitializeLMap initialize the leaves mapping
func InitializeLMap(geneT, speciesT Tree, geneLeave, speciesLeave int) {
	for i := range geneT {
		geneT[i].L = make([]*Node, 0)
	}

	for i := 0; i < geneLeave; i++ {
		for j := 0; j < speciesLeave; j++ {
			if geneT[i].label == speciesT[j].label {
				geneT[i].L = append(geneT[i].L, speciesT[j])
			}
		}
	}
}

// InitializeCosts takes in two nodes, and assign all cost relating maps to infinity
func InitializeCosts(g, s *Node) {
	g.cost[s] = float64(math.Inf(+1))
	g.speciation[s] = float64(math.Inf(+1))
	g.duplication[s] = float64(math.Inf(+1))
	g.transfer[s] = float64(math.Inf(+1))
	g.inAlt[s] = float64(math.Inf(+1))
	g.in[s] = float64(math.Inf(+1))
	g.out[s] = float64(math.Inf(+1))
}

// GetLeaveNodes takes in a tree and the numver of leave nodes. it returns
// a tree including only leaf nodes
func GetLeaveNodes(geneT Tree, leaveNum int) Tree {
	leaves := make(Tree, leaveNum)
	for i := 0; i < leaveNum; i++ {
		leaves[i] = geneT[i]
	}
	return leaves
}

// PrintTree takes in a tree and prints out the mapping for each nodes
func PrintTree(t Tree) {
	for i := range t {
		fmt.Println("Label: ", t[i].label, " Events: ", t[i].event)
	}
}
