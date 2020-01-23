package main

import (
	"fmt"
)

//Infer gene duplication and speciation events on a gene tree by refering to a species tree
//Input: Rooted binary gene tree G, rooted binary species tree S of all species in G.
//Output: G with "duplication" or "speciation" assigned to each of its internal nodes

type Tree []*Node

type Node struct {
	label                  string //species name
	child1, child2, parent *Node
	number                 int
	event                  string //label the internal nodes of gene tree, either speciation or duplication
}

//Create a global variable number that will be used later in InitializeSTree
var number int = 1

func main() {
	fmt.Println("Label gene tree events.")

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

	root := speciesTree[len(speciesTree)-1]
	LabelInternalNodeEvent(geneTree, speciesTree, root, 4)
	for i := range geneTree {
		fmt.Println(geneTree[i].label)
		fmt.Println(geneTree[i].number)
		fmt.Println(geneTree[i].event)
	}
}

//LabelInternalNodes takes in a gene tree, a species tree, a root node and the number of species, and labels the internal nodes of the gene tree with event.
func LabelInternalNodeEvent(gTree, sTree Tree, root *Node, speciesnum int) {
	//Initialize species tree.
	InitializeSTree(sTree, root)
	//Initialize gene tree.
	InitializeGTree(gTree, sTree, speciesnum)
	//Traverse the internal nodes in the gene tree and label them with either speciation or duplication event.
	TraverseGTree(gTree, sTree, speciesnum)
}

//InitializeSTree takes a species tree, and its root as inputs, and recursively numbers nodes of the species tree in preorder traversal.
//root = 1, child nodes always larger than parent node
func InitializeSTree(t Tree, root *Node) {
	root.number = number
	number++
	if root.child1 != nil {
		InitializeSTree(t, root.child1)
	}
	if root.child2 != nil {
		InitializeSTree(t, root.child2)
	}
}

//InitializeGTree takes a gene tree, a species tree, and the number of species as input.
//For each leave node in the gene tree, set its label to the number of the leave node in species tree with the matching species name.
func InitializeGTree(gTree, sTree Tree, speciesnum int) Tree {
	for i := 0; i < speciesnum; i++ {
		for j := 0; j < speciesnum; j++ {
			if gTree[i].label == sTree[j].label {
				gTree[i].number = sTree[j].number
			}
		}
	}
	return gTree
}

//TraverseGTree takes a gene tree, a species tree and the species number as input, and assignes "duplication" or "speciation" to internal nodes of gene tree
func TraverseGTree(gTree, sTree Tree, speciesnum int) {
	for i := speciesnum; i < len(gTree); i++ {
		node := gTree[i]
		a := node.child1.number
		b := node.child2.number
		for a != b {
			if a > b {
				//parent of node a in the species tree
				a = FindParentNum(sTree, a)
			} else {
				b = FindParentNum(sTree, b)
			}
		}
		//node.number cannot be lower (greater) than node.child1.number or node.child2.number
		//node.number is the Last Common Ancestor (LCA) of node.child1.number and node.child2.number
		node.number = a
		if node.number == node.child1.number || node.number == node.child2.number {
			node.event += "duplication"
		} else {
			node.event += "speciation"
		}
	}
}

//FindParentNum takes a species tree and a node number as input, and looks for the node with the same node number in the species tree.
//It returns the parent node number of the node with number a.
func FindParentNum(sTree Tree, a int) int {
	for _, node := range sTree {
		if node.number == a {
			return node.parent.number
		}
	}
	panic("You won't get here!")
}
