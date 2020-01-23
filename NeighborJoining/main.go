package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type Tree []*NodeList

type NodeList struct {
	head *Node
	len  int
}
type Node struct {
	label string
	dist  float64
	next  *Node
}
type Cluster []*Node

type TotalDist []float64

type Matrix [][]float64

func main() {
	fmt.Println("Start construct the tree using Neighbor-Joining method!")

	//First we read a distance matrix and species name from a file
	fileName := "SpeciesTree.txt"
	mtx, speciesName := ReadMatrixFromFile(fileName)

	//We then use NeighborJoining to construct the tree
	t := NeighborJoining(mtx, speciesName)

	//print out the constructed tree
	t.Print()

}

// NeighborJoining takes in distance matrix, species names and returns a tree
func NeighborJoining(mtx Matrix, speciesName []string) Tree {
	leaveLen := len(speciesName)
	// The initial tree including leavelen Nodelists, each of them is a cluster by themselves
	t := InitializeTree(speciesName)
	clusters := InitializeClusters(speciesName, t)
	u := ComputeTotalDist(mtx, leaveLen)
	k := 1 //internal nodes name counter

	//There are 2n-2 total nodes
	//It needs n-2 times to set up the tree
	for i := leaveLen; i < 2*leaveLen-2; i++ {
		//Created the corresponding altered matrix Q
		alteredMtx := AlterMatrix(mtx, u)
		//Find the minimum elements in the altered matrix
		row, col, val := FindMinElt(alteredMtx)

		//add to tree
		t = AddToTree(t, i, row, col, val, mtx, k, leaveLen, u, clusters)

		//Add the new created node to the matrix
		mtx = AddColRow(mtx, col, row)
		//Add Cluster
		clusters = append(clusters, t[i].head)

		//Delete col and row
		mtx = DelColRow(mtx, col, row)
		//Delete from cluster
		clusters = DelFromCluster(clusters, col, row)
		//Add the new total distance to U
		u = ComputeTotalDist(mtx, len(mtx))

		k++
	}

	//Final Connect
	t = FinalConnect(t, mtx, leaveLen, clusters)
	//t.Print()
	return t
}

// InitializeClusters takes in a string slice which contains species names, and
// it also takes in a tree. This function returns a cluster.
func InitializeClusters(speciesName []string, t Tree) Cluster {
	leaveLen := len(speciesName)
	clusters := make([]*Node, leaveLen)
	for i := 0; i < leaveLen; i++ {
		clusters[i] = t[i].head
	}
	return clusters
}

// DelFromCluster takes in a cluster, and col, row index, and delete the
// corresponding column and row. It returns the cluster after deletion.
func DelFromCluster(clusters Cluster, col, row int) Cluster {

	clusters = append(clusters[:col], clusters[col+1:]...)
	clusters = append(clusters[:row], clusters[row+1:]...)

	return clusters
}

// AddToTree takes in a tree, the col, row index, the index of the new create node
// in the tree "i", the matrixthat we read from file, the total distance of each
// row of the matrix and the clusters that we created. It added the newly created
// internal node to our tree, and returns the tree after addition.
func AddToTree(t Tree, i, row, col int, val float64, mtx Matrix, k, leaveLen int, u TotalDist, clusters Cluster) Tree {
	//	fmt.Println("Adding To Trees.")
	//add a new node
	var newNode Node
	newNode.label = "Internal" + strconv.Itoa(k)
	t[i].head = &newNode
	//new node points to the chosen leaves
	//contribute to the limb length
	distance := Dist(row, col, u)
	//two nodes for the chosen leaves
	var newNode1 Node
	newNode1.label = clusters[row].label
	newNode1.dist = 0.5 * (mtx[row][col] + distance)
	var newNode2 Node
	newNode2.label = clusters[col].label
	newNode2.dist = 0.5 * (mtx[row][col] - distance)
	p := t[i].head
	for p.next != nil {
		p = p.next
	}
	p.next = &newNode1
	newNode1.next = &newNode2

	t[i].len = t[i].len + 2

	//leaves point to internal nodes
	var newInternal Node
	newInternal.label = t[i].head.label
	newInternal.dist = newNode1.dist
	d := clusters[row]
	for d.next != nil {
		d = d.next
	}
	d.next = &newInternal

	var newInternal2 Node
	newInternal2.label = t[i].head.label
	newInternal2.dist = newNode2.dist
	h := clusters[col]
	for h.next != nil {
		h = h.next
	}
	h.next = &newInternal2

	//t.Print()
	return t
}

// FinalConnect takes in a tree, the distance matrix that we read from file,
// the leave length of the tree, and the clusters of nodes. It connect the final
// pair of nodes without creating new internal node. It returns a new tree
// after connecting the final nodes.
func FinalConnect(t Tree, mtx Matrix, leaveLen int, clusters Cluster) Tree {

	var node1 Node
	node1.label = clusters[0].label
	node1.dist = mtx[0][1]
	p := clusters[1]
	for p.next != nil {
		p = p.next
	}
	p.next = &node1

	var node2 Node
	node2.label = clusters[1].label
	node2.dist = mtx[0][1]
	k := clusters[0]
	for k.next != nil {
		k = k.next
	}
	k.next = &node2

	return t
}

// InitializeTree initializes the unrooted tree
// We allocate 2n-2 total numbers of arrays of list
// The leaves are the first part of the tree
func InitializeTree(speciesname []string) Tree {
	//fmt.Println("Initializing Tree!")
	speciesnum := len(speciesname)
	totalnodes := 2*speciesnum - 2
	//fmt.Println("Totalnodes ", totalnodes)
	var t Tree = make([]*NodeList, totalnodes)
	for i := 0; i < speciesnum; i++ {
		var newList NodeList
		var newNode Node
		newNode.label = speciesname[i]
		//fmt.Println("NEWNODE: ", newNode.label, "Memory ", &newNode)
		newList.head = &newNode
		newList.len = 0
		t[i] = &newList
	}
	for i := speciesnum; i < totalnodes; i++ {
		var newList NodeList
		var newNode Node
		newList.head = &newNode
		newList.len = 0
		t[i] = &newList
	}
	//t.Print()
	return t
}

// ComputeTotalDist takes in a distance matrix, and the number of leaves and
// compute the total distances for each row of the matrix. It returns a vector
// of total distances.
func ComputeTotalDist(mtx Matrix, leaveLen int) TotalDist {
	//fmt.Println("Initializing Total Distance.")
	u := make(TotalDist, leaveLen)

	for i := 0; i < leaveLen; i++ {
		u[i] = SumRow(mtx, i)
	}

	return u
}

// SumRow sums up all elements in one row in a given matrix
func SumRow(mtx Matrix, i int) float64 {
	length := len(mtx[0])
	sum := 0.0
	for j := 0; j < length; j++ {
		sum += mtx[i][j]
	}
	return sum
}

// AlterMatrix takes in a matrix and edits the numbers and return a new matrix
func AlterMatrix(mtx Matrix, u TotalDist) Matrix {
	//	fmt.Println("Alter Matrix!")
	rownum := len(mtx)
	colnum := len(mtx[0])
	alteredMatrix := make(Matrix, rownum)
	for i := range alteredMatrix {
		alteredMatrix[i] = make([]float64, colnum)
	}
	for i := 0; i < rownum; i++ {
		for j := 0; j < colnum; j++ {
			if i != j {
				alteredMatrix[i][j] = float64((len(mtx)-2))*mtx[i][j] - u[i] - u[j]
			} else {
				alteredMatrix[i][j] = 0.0
			}
		}
	}

	return alteredMatrix
}

// FindMinElt takes in the altered matrix that we have and returns the row, col and
// val of the minimum value
func FindMinElt(alteredMtx Matrix) (int, int, float64) {
	//	fmt.Println("Finding the Minimum.")
	if len(alteredMtx) <= 1 || len(alteredMtx[0]) <= 1 {
		panic("We gave too small a matrix to FindMinElt")
	}
	row := 0
	col := 1
	minVal := alteredMtx[row][col]

	for i := 0; i < len(alteredMtx); i++ {
		//range over column, starting at j=i+1 so we don't hit main diagonal (0)
		for j := i + 1; j < len(alteredMtx[0]); j++ {
			if alteredMtx[i][j] < minVal {
				row = i
				col = j
				minVal = alteredMtx[i][j]
			}
		}
	}

	return row, col, minVal
}

// Dist takes in the minimum value's row and col index, the total number of leaves,
// and the TotalDist of each leave, and returns a distance
func Dist(row, col int, u TotalDist) float64 {
	dist := float64((u[row] - u[col]) / float64(len(u)-2))
	//fmt.Println("Dist ", dist)
	return dist
}

//AddColRow takes in a matrix, the chosen col and row and returns a new matrix
func AddColRow(mtx Matrix, col, row int) Matrix {
	//fmt.Println("ADD COL AND ROW.")
	n := len(mtx)
	newRow := make([]float64, n+1)
	//all values will be 0.0 by default, so we only need to update values that will stick around
	for j := 0; j < n; j++ {
		if j != row && j != col {
			newRow[j] = 0.5 * (mtx[row][j] + mtx[col][j] - mtx[row][col])
		}
	}
	mtx = append(mtx, newRow)

	//now add the last column
	for i := 0; i < n; i++ {
		mtx[i] = append(mtx[i], newRow[i])
	}

	return mtx
}

//DelColRow delete the redundant col and row in the matrix
func DelColRow(mtx Matrix, col, row int) Matrix {
	//fmt.Println("Delete Col and Row.")
	//first, delete appropriate row
	mtx = append(mtx[:col], mtx[col+1:]...) //delete col first
	mtx = append(mtx[:row], mtx[row+1:]...)

	//delete appropriate columns
	for i := range mtx {
		mtx[i] = append(mtx[i][:col], mtx[i][col+1:]...)
		mtx[i] = append(mtx[i][:row], mtx[i][row+1:]...)
	}

	return mtx
}

//ReadMatrixFromFile reads from the file that we take in and returns the matrix,
// and the species name in the file
func ReadMatrixFromFile(fileName string) (Matrix, []string) {
	file, err := os.Open(fileName)
	if err != nil {
		fmt.Println("Error: couldn't open the file")
		os.Exit(1)
	}
	lines := make([]string, 0)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	if scanner.Err() != nil {
		fmt.Println("Sorry: there was some kind of error during the file reading")
		os.Exit(1)
	}
	file.Close()

	mtx := make(Matrix, 0)
	speciesNames := make([]string, 0)

	for idx := range lines {
		if idx >= 1 {
			row := make([]float64, 0)
			nums := strings.Split(lines[idx], "\t")
			for i, num := range nums {
				if i == 0 {
					speciesNames = append(speciesNames, num)
				} else {
					n, err := strconv.ParseFloat(num, 64)
					if err != nil {
						fmt.Println(err)
						fmt.Println("Error: Wrong format of matrix!")
						os.Exit(1)
					}
					row = append(row, n)
				}
			}
			mtx = append(mtx, row)
		}
	}
	return mtx, speciesNames
}

//Print prints out a tree
func (t Tree) Print() {
	for i := range t {
		fmt.Print("Head ", t[i].head.label, " is connected with")
		p := t[i].head
		for p.next != nil {
			fmt.Print(" ", p.next.label, " Distance: ", p.next.dist, " ")
			p = p.next
		}
		fmt.Println()
	}
}
