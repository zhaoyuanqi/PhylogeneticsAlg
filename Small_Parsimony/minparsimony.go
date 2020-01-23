package main

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"
)

//Small Parsimony Problem: Find the most parsimonious labeling of the internal nodes of a rooted tree
//Input: A rooted binary tree with each leaf labeled by a string of length m
//Output: A labeling of all other nodes of the tree by strings of length m that minimizes the tree's parsimony score

//Input tree structure
type Tree []*Node

type Node struct {
	score                  map[int]float64
	label                  string //DNA sequence
	parent, child1, child2 *Node
}

type Matrix [][]float64

func main() {
	rand.Seed(time.Now().UnixNano())
	nucList := []string{"A", "T", "C", "G", "-"}

	//Read a parsimony score matrix.
	//The row and col orders of the matrix are the same as the order of nucleotides in the nucList above
	filename := os.Args[1]

	var t Tree
	t = make([]*Node, 7)
	var v0, v1, v2, v3, v4, v5, v6 Node
	v4 = Node{child1: &v0, child2: &v1}
	v5 = Node{child1: &v2, child2: &v3}
	v6 = Node{child1: &v4, child2: &v5}
	t[0] = &v0
	t[1] = &v1
	t[2] = &v2
	t[3] = &v3
	t[4] = &v4
	t[5] = &v5
	t[6] = &v6

	leaveseq := []string{"CG", "CG", "AT", "CC"}

	MinimumParsimony(leaveseq, filename, t, nucList)
	for n := range t {
		fmt.Println(t[n].label)
		fmt.Println(t[n].score)
	}
}

//MinimumParsimony assumes that nucleotides in the sequence are independent and seek trees with the lowest possible parsimony score.
//leaveseq is a slice of string, containing DNA sequences of all leave nodes (present day species).
func MinimumParsimony(leaveseq []string, filename string, t Tree, nucList []string) Tree {
	mtx := ReadMatrix(filename) //read in the score matrix
	InitializeTree(t, leaveseq)
	//work with one character of each string at a time
	for i := range leaveseq[0] {
		BaseMinPars(mtx, t, i, nucList)
	}
	return t
}

//ReadMatrix takes a parsimony score file as input and store it in a 2-D matrix
func ReadMatrix(filename string) Matrix {
	score := make(Matrix, 0)
	//Read from file
	f, err := os.Open(filename)
	if err != nil {
		panic("Issue in reading the score file")
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		nums := strings.Split(line, "\t")
		row := make([]float64, 0)
		for _, num := range nums {
			//change each string to a float64
			n, err := strconv.ParseFloat(num, 64)
			if err != nil {
				fmt.Println(err)
				panic("Issue in format of matrix")
			}
			row = append(row, n)
		}
		score = append(score, row)
	}
	return score
}

//InitializeTree takes a tree and sequences of leave nodes as inputs and adds the label and score to each leave node.
func InitializeTree(t Tree, leaveseq []string) Tree {
	leaveLen := len(leaveseq)
	//set labels for leave nodes as nucleotide sequences
	for i := 0; i < leaveLen; i++ {
		t[i].label = leaveseq[i]
		t[i].score = make(map[int]float64, 5)
	}
	//will add labels for root and internal nodes
	for j := leaveLen; j < len(t); j++ {
		t[j].score = make(map[int]float64, 5)
	}
	//Add parent to all the nodes in tree except the root node
	AddParent(t)
	return t
}

//AddParent takes a tree as input and labels every node except root.
func AddParent(t Tree) Tree {
	leaveLen := (len(t) + 1) / 2
	treeLen := len(t)
	for i := leaveLen; i < treeLen; i++ {
		node1 := t[i].child1
		node2 := t[i].child2
		node1.parent = t[i]
		node2.parent = t[i]
	}
	return t
}

//BaseMinPars takes in a score matrix, tree, position, and nucList as inputs and returns the minimum parsimony label of internal nodes at this position
func BaseMinPars(mtx Matrix, t Tree, i int, nucList []string) Tree {
	//Assign scores for leave nodes
	InitialScore(t, i, nucList)
	//Assign scores for internal nodes
	InternalScore(t, mtx)
	//Backtrack and find the nucleotide label for internal nodes
	BackTrack(t, nucList)
	return t
}

//InitialScore takes in a tree, position, and nucList as inputs and returns the initial score map of each node.
func InitialScore(t Tree, i int, nucList []string) Tree {
	leaveLen := (len(t) + 1) / 2
	//set scores of leave nodes equal to 0 for the nucleotide at the current position and infinity for all other nucleotides.
	for k := 0; k < leaveLen; k++ {
		for idx, nuc := range nucList {
			if t[k].label[i:i+1] == nuc {
				t[k].score[idx] = 0.0
			} else {
				t[k].score[idx] = math.Inf(1)
			}
		}
	}
	//set all scores of internal and root nodes to infinity
	for l := leaveLen; l < len(t); l++ {
		for idx := range nucList {
			t[l].score[idx] = math.Inf(1)
		}
	}
	return t
}

//InternalScore takes a tree and parsimony score matrix as inputs and calculates the score maps for internal and root nodes.
func InternalScore(t Tree, mtx Matrix) Tree {
	leaveLen := (len(t) + 1) / 2
	//loop through each internal node
	for i := leaveLen; i < len(t); i++ {
		//l is a list of keys with minimum score at the node, s is the minimum score
		l1, s1 := FindMin(t[i].child1)
		l2, s2 := FindMin(t[i].child2)

		//use l and s calculated above to assign scores of each key in the score map of the internal node
		for j := range t[i].score {
			//min is the minimum of score of parent to child(s) transition [child may have more than one nucleotides (keys) with minimum score]
			min1 := MinTran(l1, j, mtx)
			min2 := MinTran(l2, j, mtx)
			t[i].score[j] = s1 + s2 + min1 + min2
		}
	}
	return t
}

//FindMin takes a node as input and returns the minimum score (value) in the score map and corresponding key(s)
func FindMin(node *Node) ([]int, float64) {
	k := []int{}
	v := []float64{}
	label := []int{}
	for key, value := range node.score {
		k = append(k, key)
		v = append(v, value)
	}
	minScore := v[0]
	for i := 1; i < len(v); i++ {
		if v[i] < minScore {
			minScore = v[i]
		}
	}
	for idx, val := range v {
		if val == minScore {
			label = append(label, k[idx])
		}
	}
	return label, minScore
}

//MinTran takes a slice of labels fro FindMin, a key of the score map of an internal node, and a score matrix and returns the minimum transition score
func MinTran(nucs []int, j int, mtx Matrix) float64 {
	min := math.Inf(1)
	for k := range nucs {
		score := mtx[j][nucs[k]]
		if score < min {
			min = score
		}
	}
	return min
}

//BackTrack takes a stree and nucList as inputs and assign labels for all internal and root nodes.
func BackTrack(t Tree, nucList []string) Tree {
	//Assign the label of root node first
	root := t[len(t)-1]
	lr, _ := FindMin(root)
	if len(lr) == 1 {
		root.label += nucList[lr[0]]
	} else {
		//if root has more than one possible nucleotides which would result in a tree with the same parsimony score, we would randomly choose one from them
		idx := rand.Intn(len(lr))
		root.label += nucList[lr[idx]]
	}
	//Use the tree with assigned root label to assign labels of internal nodes
	AddIntNuc(t, root, nucList)
	return t
}

//AddIntNuc takes in a tree, root node, and nucList as inputs and returns a tree with labeled internal nodes
func AddIntNuc(t Tree, root *Node, nucList []string) {
	if root == nil {
		panic("There is no root in this tree")
	}
	//Recursion would stop if we reach leave nodes
	if root.child1 != nil {
		AddIntNuc(t, root.child1, nucList)
	}

	if root.child2 != nil {
		AddIntNuc(t, root.child2, nucList)
	}

	//only internal nodes satisfy the three conditions below
	if root.child1 != nil && root.child2 != nil && root.parent != nil {
		l, _ := FindMin(root)
		//if there is only one key in the score map with minimum score, then just record it as the label
		if len(l) == 1 {
			root.label += nucList[l[0]]
		} else {
			// if there is more than more key with minimum score, then we have to check its parent label and use the same label.
			parentLabel := root.parent.label
			root.label += parentLabel[len(parentLabel)-1:]
		}
	}
}
