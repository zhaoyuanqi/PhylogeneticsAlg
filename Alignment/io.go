package main

import (
	"bufio"
	"fmt"
	"os"
)

//ReadGeneFromFile takes a FASTA file as input and returns the gene sequence in the file as a string.
func ReadGeneFromFile(filename string) string {
	file, err := os.Open(filename)
	if err != nil {
		fmt.Println("Error: couldn't open the file")
		os.Exit(1)
	}
	defer file.Close()

	gene := ""
	scanner := bufio.NewScanner(file)
	//skip first line and move to the next token
	scanner.Scan()
	//read line by line
	for scanner.Scan() {
		currentLine := scanner.Text()
		gene += currentLine
	}
	return gene
}

//WriteMatrixToFile takes a slifce of species names, a distance matrix, and a output filename
//and generates a file with the format that will be used to build a phylogenetic tree.
func WriteMatrixToFile(speciesName []string, matrix [][]int, filename string) {
	file, err := os.Create(filename)
	if err != nil {
		fmt.Println("Sorry: couldn't create the file!")
		os.Exit(1)
	}
	defer file.Close()

	//fistline of the file is the dimension of the matrix
	fmt.Fprintln(file, len(speciesName))
	for i := range matrix {
		fmt.Fprint(file, speciesName[i], "\t")
		for j := 0; j < len(matrix[i]); j++ {
			if j < len(matrix[i])-1 {
				fmt.Fprint(file, matrix[i][j], "\t")
			} else {
				fmt.Fprint(file, matrix[i][j])
				//move the pointer to next line
				fmt.Fprintln(file)
			}
		}
	}
}
