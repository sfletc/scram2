package main

import (
	"flag"
	"fmt"
)

func clInput() (*string, *string, *string, *string, *string, *bool, *bool, *int, *int, *int) {
	refFile := flag.String("refFile", "", "Path to target FASTA format reference file (required)")
	readFiles := flag.String("readFiles", "", "Path to target FASTQ format comma-separated read files (required)")
	outPrefix := flag.String("outPrefix", "", "Output folder/file prefix (required)")
	alignLens := flag.String("alignLens", "21,22,24", "Comma-separated small RNA lenghts to exact match (default=21,22,24)")
	readFileType := flag.String("readFileType", "fastq", "Readfile type - fastq (default) or fasta")
	split := flag.Bool("split", true, "Split read alignment count for each position by the number of positions a read aligns too (default=True)")
	norm := flag.Bool("norm", true, "Normalise reads to Reads Per Million Reads (RPMR) from each input file (default=True)")
	minLenNorm := flag.Int("minLenNorm", 18, "Minimum size read to include in total normalization counts (default=18)")
	maxLenNorm := flag.Int("maxLenNorm", 32, "Maximum size read to include in total normalization counts (default=32)")
	minCountNorm := flag.Int("minCountNorm", 32, "Minimum abundence of a read to include in total normalization counts (default = 1)")
	flag.Parse()
	return refFile, readFiles, alignLens, readFileType, outPrefix, split, norm, minLenNorm, maxLenNorm, minCountNorm
}

func main() {

	refFile, readFiles, alignLens, readFileType, outPrefix, split, norm, minLenNorm, maxLenNorm, minCountNorm := clInput()
	fmt.Println(refFile, readFiles, alignLens, readFileType, outPrefix, split, norm, minLenNorm, maxLenNorm, minCountNorm)

}
