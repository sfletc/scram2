package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"errors"
	"fmt"
	"os"
	"path"
	"strings"
	"sync"

	"github.com/dustin/go-humanize"
)

func ErrorShutdown() {
	fmt.Println("\nExiting")
	os.Exit(1)
}

// // SeqLoad loads 1 or more small RNA seq. read files.
// // It returns a map with a read sequence as key and a MeanSe struct (normalised or raw read mean and standard error) as a value.
// // Little format checking is  performed.  It is required that the input file is correctly formatted.
// func SeqLoad(seqFiles []string, fileType string, adapter string, minLen int, maxLen int,
// 	minCount float64, noNorm bool) map[string]interface{} {
// 	noOfFiles, srnaMaps := loadFiles(seqFiles, fileType, minLen, maxLen, minCount, noNorm, adapter)
// 	seqMapAllCounts, _ := compileCounts(srnaMaps, noOfFiles, minCount)
// 	seqMap := calcMeanSe(seqMapAllCounts, noOfFiles)
// 	return seqMap
// }

// IndvSeqLoad loads 1 or more small RNA seq. read files.
// It returns a map with a read sequence as key and a slice of normalized or raw individual read counts as a value.
// Little format checking is  performed.  It is required that the input file is correctly formatted.
func IndvSeqLoad(seqFiles []string, fileType string, adapter string, minLen int, maxLen int,
	minCount float64, noNorm bool) (map[string]interface{}, []string) {
	noOfFiles, srnaMaps := loadFiles(seqFiles, fileType, minLen, maxLen, minCount, noNorm, adapter)
	seqMapAllCounts, loadOrder := compileCounts(srnaMaps, noOfFiles, minCount)
	return seqMapAllCounts, loadOrder
}

// loadFiles loads replicate read files into a channel - ref_name map of read / count pairs
func loadFiles(seqFiles []string, fileType string, minLen int, maxLen int, minCount float64, noNorm bool, adapter string) (int, chan map[string]map[string]float64) {
	wg := &sync.WaitGroup{}
	noOfFiles := len(seqFiles)
	wg.Add(noOfFiles)
	fileNames := make(chan string, len(seqFiles))
	for _, fileName := range seqFiles {
		fileNames <- fileName
	}
	close(fileNames)
	srnaMaps := make(chan map[string]map[string]float64, len(seqFiles))
	for a := 0; a < len(seqFiles); a++ {
		switch {
		case fileType == "fa":
			if a == 0 {
				fmt.Println("\nAttempting to load read files in FASTA format")
			}
			go loadFastx(fileNames, []byte(">"), adapter, srnaMaps, minLen, maxLen, minCount, noNorm, wg)
		case fileType == "fq":
			if a == 0 {
				fmt.Println("\nAttempting to load read files in FASTQ format")
			}
			go loadFastx(fileNames, []byte("@"), adapter, srnaMaps, minLen, maxLen, minCount, noNorm, wg)
		}
	}
	go func(cs chan map[string]map[string]float64, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(srnaMaps, wg)
	return noOfFiles, srnaMaps
}

// loadFastx loads a single FASTA or FASTQ file and return map of read sequence as key and normalised RPMR count as
// value. Trim adapter from 3' end using up to 12 nt of 5' end of adapter as seed if required
func loadFastx(fileNames chan string, firstChar []byte, adapter string, srnaMaps chan map[string]map[string]float64,
	minLen int, maxLen int, minCount float64, noNorm bool, wg *sync.WaitGroup) {
	srnaMap := make(map[string]float64)
	var totalCount float64
	fileName := <-fileNames
	f, err := os.Open(fileName)
	if err != nil {
		fmt.Println("\nCan't load read file " + fileName)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	if fileName[len(fileName)-2:] == "gz" {
		gz, err := gzip.NewReader(f)
		if err != nil {
			fmt.Println("\nCan't decompress read file " + fileName)
			os.Exit(1)
		}
		defer gz.Close()
		scanner = bufio.NewScanner(gz)
	}

	seqNext := false
	for scanner.Scan() {
		fasta_line := scanner.Bytes()
		switch {
		case bytes.Equal(fasta_line[:1], firstChar):
			seqNext = true
		case seqNext == true && len(fasta_line) >= minLen && len(fasta_line) <= maxLen:
			srnaMap, totalCount, seqNext = addFullLengthRead(srnaMap, fasta_line, totalCount, seqNext)
		}
	}
	srnaMap, totalCount = removeReadsBelowMin(minCount, srnaMap, totalCount)

	if noNorm == false {
		srnaMap = rpmrNormalize(srnaMap, totalCount)
	}
	final_map := map[string]map[string]float64{fileName: srnaMap}
	srnaMaps <- final_map
	fmt.Println(fileName + " - " + humanize.Comma(int64(totalCount)) + " reads processed")
	wg.Done()
}

// Add full-length read to the srna map
func addFullLengthRead(srna_map map[string]float64, fasta_line []byte, total_count float64,
	seq_next bool) (map[string]float64, float64, bool) {
	if srna_count, ok := srna_map[string(fasta_line)]; ok {
		srna_map[string(fasta_line)] = srna_count + 1.0
		total_count += 1.0
	} else {
		srna_map[string(fasta_line)] = 1.0
		total_count += 1.0
	}
	seq_next = false
	return srna_map, total_count, seq_next
}

// Remove reads with count below the stated minimum for the srna_map
func removeReadsBelowMin(minCount float64, srnaMap map[string]float64, totalCount float64) (map[string]float64,
	float64) {
	if minCount > 1 {
		for srna, srnaCount := range srnaMap {
			if srnaCount < minCount {
				delete(srnaMap, srna)
				totalCount -= srnaCount
			}
		}
	}
	return srnaMap, totalCount
}

// Reads per million reads normalization of an input read library
func rpmrNormalize(srnaMap map[string]float64, total_count float64) map[string]float64 {
	for srna, srnaCount := range srnaMap {
		srnaMap[srna] = 1000000 * srnaCount / total_count
	}
	return srnaMap
}

// Checks for error in collapsed fasta header
func checkHeaderError(headerLine []string, file_name string) error {
	if len(headerLine) < 2 || len(headerLine) > 2 {
		return errors.New("\n" + file_name + " is incorrectly formatted")
	}
	return nil
}

// Compile_counts generates a map with read seq as key and a slice of normalised counts for each read file
func compileCounts(srna_maps chan map[string]map[string]float64, no_of_files int, min_count float64) (map[string]interface{}, []string) {
	// map [srna:[count1,count2....], ...]
	seq_map_all_counts := make(map[string]interface{})
	var load_order []string
	pos := 0
	for singleSeqMap := range srna_maps {
		for file, seqMap := range singleSeqMap {
			load_order = append(load_order, path.Base(file))
			for srna, count := range seqMap {
				if _, ok := seq_map_all_counts[srna]; ok {
					// a:= append(*seq_map_all_counts[srna].(*[]float64), count)
					a := seq_map_all_counts[srna].(*[]float64)
					(*a)[pos] = count
					seq_map_all_counts[srna] = a
				} else {
					a := make([]float64, no_of_files)
					a[pos] = count
					seq_map_all_counts[srna] = &a

				}
			}
			pos++
		}
	}
	if min_count > 1 {
		removeUnderMinCount(seq_map_all_counts)
	}
	return seq_map_all_counts, load_order
}

// Remove read if its count is under the specified minimum
func removeUnderMinCount(seq_map_all_counts map[string]interface{}) {
	for srna := range seq_map_all_counts {
		counts := seq_map_all_counts[srna]
		for _, i := range *counts.(*[]float64) {
			// If using a min_count > 1, unless the srna is present in all libraries, it's removed so as not
			// to generate spurious means and standard errors
			if i == 0.0 {
				delete(seq_map_all_counts, srna)
			}
		}
	}
}

// HeaderRef is a struct comprising a reference sequence header, seques and reverse complement
type HeaderRef struct {
	Header     string
	Seq        string
	ReverseSeq string
}

// RefLoad loads a reference sequence DNA file (FASTA format).
// It returns a slice of HeaderRef structs (individual reference header, sequence and reverse complement).
func RefLoad(refFile string) []*HeaderRef {
	var totalLength int
	var refSlice []*HeaderRef
	var singleHeaderRef *HeaderRef
	var header string
	var refSeq bytes.Buffer
	f, err := os.Open(refFile)
	defer f.Close()
	if err != nil {
		fmt.Println("Problem opening fasta reference file " + refFile)
		ErrorShutdown()
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fastaLine := scanner.Text()
		switch {
		case strings.HasPrefix(fastaLine, ">"):
			seq := refSeq.String()
			singleHeaderRef = &HeaderRef{header, seq, reverseComplement(seq)}
			refSlice = append(refSlice, singleHeaderRef)
			header = fastaLine[1:]
			refSeq.Reset()
		case len(fastaLine) != 0:
			refSeq.WriteString(strings.ToUpper(fastaLine))
			totalLength += len(fastaLine)
		}
	}
	seq := refSeq.String()
	singleHeaderRef = &HeaderRef{header, seq, reverseComplement(seq)}
	refSlice = append(refSlice, singleHeaderRef)
	refSlice = refSlice[1:]

	fmt.Println("No. of reference sequences: ", len(refSlice))
	fmt.Println("Combined length of reference sequences: " + humanize.Comma(int64(totalLength)) + " nt")
	return refSlice
}

// Reverse complements a DNA sequence
func reverseComplement(seq string) string {
	complement := map[rune]rune{
		'A': 'T',
		'C': 'G',
		'G': 'C',
		'T': 'A',
		'N': 'N',
	}
	runes := []rune(seq)
	var result bytes.Buffer
	for i := len(runes) - 1; i >= 0; i-- {
		result.WriteRune(complement[runes[i]])
	}
	return result.String()
}