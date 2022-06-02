// BSD 3-Clause License

// Copyright (c) 2022, Stephen Fletcher
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

package scramPkg

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
)

//CompareNoSplitCounts takes and alignment map and returns a map with the ref_header as key and the
//mean_se (mean and standard error) of aligned reads for that ref seq as value.  Read counts are NOT split by the number
//of times a read aligns to all reference sequences.
func CompareNoSplitCounts(alignmentMap map[string]map[string][]int, seqMap map[string]interface{}) map[string]interface{} {
	//cdp_alignment_map := make(map[string]meanSe)
	cdpAlignmentMap := make(map[string]interface{})
	for header, alignment := range alignmentMap {
		var headerMeanCounts float64
		var meanCountsErr []float64
		var headerCounts []float64
		firstPos := true
		for srna, pos := range alignment {

			switch v := seqMap[srna].(type) {

			case *meanSe:
				headerMeanCounts += v.Mean * float64(len(pos))
				meanCountsErr = append(meanCountsErr, v.Se*float64(len(pos)))
			case *[]float64:
				if firstPos {
					headerCounts = make([]float64, len(*v))
					firstPos = false
				}
				countPos := 0
				for _, i := range *v {
					headerCounts[countPos] += i * float64(len(pos))
					countPos++
				}
			}
		}
		if meanCountsErr != nil {
			cdpAlignmentMap = calcHeaderMeanSe(meanCountsErr, headerMeanCounts, cdpAlignmentMap, header)
		} else {
			cdpAlignmentMap[header] = headerCounts
		}
	}
	return cdpAlignmentMap
}

//CompareSplitCounts takes and alignment map and returns a map with the ref_header as key and the
//mean_se (mean and standard error) of aligned reads for that ref seq as value.  Read counts are split by the number
//of times a read aligns to all reference sequences.
func CompareSplitCounts(alignmentMap map[string]map[string][]int, seqMap map[string]interface{}) map[string]interface{} {
	cdpAlignmentMap := make(map[string]interface{})
	//Calc. no. of times each read aligns
	srnaAlignmentMap := calcTimesReadAligns(alignmentMap)
	for header, alignment := range alignmentMap {
		var headerMeanCounts float64
		var meanCountsErr []float64
		var headerCounts []float64
		firstPos := true
		for srna, pos := range alignment {
			switch v := seqMap[srna].(type) {
			case *meanSe:
				headerMeanCounts += v.Mean * float64(len(pos)) / float64(srnaAlignmentMap[srna])
				// should be err=sqrt(x^2/n) for each alignment ??
				// They are perfectly correlated (dependent), so maybe not?
				meanCountsErr = append(meanCountsErr,
					v.Se*float64(len(pos))/float64(srnaAlignmentMap[srna]))
				//counts_err = append(counts_err, (math.Sqrt((seq_map[srna].Se*seq_map[srna].Se)/
				// float64(srna_alignment_map[srna])))*float64(len(pos)))
			case *[]float64:
				if firstPos {
					headerCounts = make([]float64, len(*v))
					firstPos = false
				}
				countPos := 0
				for _, i := range *v {
					headerCounts[countPos] += i * float64(len(pos)) / float64(srnaAlignmentMap[srna])
					countPos++
				}

			}
		}
		if meanCountsErr != nil {
			cdpAlignmentMap = calcHeaderMeanSe(meanCountsErr, headerMeanCounts, cdpAlignmentMap, header)
		} else {
			cdpAlignmentMap[header] = headerCounts
		}
	}
	return cdpAlignmentMap
}

//Calculates the number of times an aligned read aligns
func calcTimesReadAligns(alignmentMap map[string]map[string][]int) map[string]int {
	srnaAlignmentMap := make(map[string]int)
	for _, alignment := range alignmentMap {
		for srna, pos := range alignment {
			if _, ok := srnaAlignmentMap[srna]; ok {
				srnaAlignmentMap[srna] += len(pos)
			} else {
				srnaAlignmentMap[srna] = len(pos)
			}
		}
	}
	return srnaAlignmentMap
}

//Calculates the mean  and se of alignments to a header
func calcHeaderMeanSe(countsErr []float64, counts float64, cdpAlignmentMap map[string]interface{},
	header string) map[string]interface{} {
	var errSqSum float64
	for _, errs := range countsErr {
		errSqSum += errs * errs
	}
	refAlignErr := math.Sqrt(errSqSum)
	cdpCountsAndErr := meanSe{counts, refAlignErr}
	cdpAlignmentMap[header] = cdpCountsAndErr
	return cdpAlignmentMap
}

//Compare combines individual alignments for set sets of sequences (treatments).  It returns a map of ref header
//as key and a slice of set 1 mean/se and set2 mean/se as value.
func Compare(countsMap1 map[string]interface{}, countsMap2 map[string]interface{}) map[string]interface{} {
	cdpFinalMap := make(map[string]interface{})
	//TODO: This only includes if BOTH maps have a non-zero count alignment - perhaps should mod?
	for header, countStats := range countsMap1 {
		if countStats2, ok := countsMap2[header]; ok {
			switch v := countStats.(type) {
			case meanSe:
				cdpFinalMap[header] = compMeanSeOutput{}
				cdpFinalMap[header] = cdpFinalMap[header].(compMeanSeOutput).append(
					v.Mean,
					v.Se,
					countStats2.(meanSe).Mean,
					countStats2.(meanSe).Se)
			case []float64:

				cdpFinalMap[header] = countsOutput{}
				pos := 0
				for pos < len(v) {
					cdpFinalMap[header] = cdpFinalMap[header].(countsOutput).append(v[pos])
					pos++
				}
				pos = 0
				for pos < len(countStats2.([]float64)) {
					cdpFinalMap[header] = cdpFinalMap[header].(countsOutput).append(countStats2.([]float64)[pos])
					pos++
				}

			}
		}
	}
	return cdpFinalMap
}

type compMeanSeOutput struct {
	output []float64
}

func (f compMeanSeOutput) append(a ...float64) compMeanSeOutput {
	for _, i := range a {
		f.output = append(f.output, i)
	}
	return f
}

type countsOutput struct {
	output []float64
}

func (f countsOutput) append(a ...float64) countsOutput {
	for _, i := range a {
		f.output = append(f.output, i)
	}
	return f
}

func MirnaCompare(mirnaAlignmentMap1 map[string]interface{},
	mirnaAlignmentMap2 map[string]interface{}, noSplit bool) map[string]interface{} {

	cdpFinalMap := make(map[string]interface{})
	for header, countStats := range mirnaAlignmentMap1 {
		if countStats2, ok := mirnaAlignmentMap2[header]; ok {
			//assume count_stats.(type) ==  count_stats_2.(type)
			switch v := countStats.(type) {
			case *mean_se_dup:
				cdpFinalMap[header] = compMeanSeOutput{}
				switch {
				case noSplit == true:

					cdpFinalMap[header] = cdpFinalMap[header].(compMeanSeOutput).append(
						v.mean_se.(*meanSe).Mean,
						v.mean_se.(*meanSe).Se,
						countStats2.(*mean_se_dup).mean_se.(*meanSe).Mean,
						countStats2.(*mean_se_dup).mean_se.(*meanSe).Se)
				default:

					cdpFinalMap[header] = cdpFinalMap[header].(compMeanSeOutput).append(
						v.mean_se.(*meanSe).Mean/v.dup,
						v.mean_se.(*meanSe).Se/v.dup,
						countStats2.(*mean_se_dup).mean_se.(*meanSe).Mean/countStats2.(*mean_se_dup).dup,
						countStats2.(*mean_se_dup).mean_se.(*meanSe).Se/countStats2.(*mean_se_dup).dup)
				}
			case *counts_dup:
				cdpFinalMap[header] = countsOutput{}
				switch {
				case noSplit == true:
					cdpFinalMap[header] = cdpFinalMap[header].(countsOutput).append(v.counts...)
					cdpFinalMap[header] = cdpFinalMap[header].(countsOutput).append(countStats2.(*counts_dup).counts...)
				default:
					pos := 0
					for pos < len(v.counts) {
						cdpFinalMap[header] = cdpFinalMap[header].(countsOutput).append(v.counts[pos] / v.dup)
						pos++
					}
					pos = 0
					for pos < len(countStats2.(*counts_dup).counts) {
						cdpFinalMap[header] = cdpFinalMap[header].(countsOutput).append(countStats2.(*counts_dup).counts[pos] / v.dup)
						pos++
					}
				}
			}
		}
	}
	return cdpFinalMap
}

//CompareToCsv writes the output to a csv file.
func CompareToCsv(cdpAlignmentMap map[string]interface{}, nt int, outPrefix string, aFileOrder []string, bFileOrder []string) {
	firstLine := true
	var alignments [][]string

	for header, countStats := range cdpAlignmentMap {
		switch v := countStats.(type) {
		case compMeanSeOutput:
			if firstLine {
				alignments = append(alignments, []string{"Header", "Mean count 1", "Std. err 1", "Mean count 2",
					"Std. err 2"})
				firstLine = false
			} else {
				alignment := []string{header,
					strconv.FormatFloat(v.output[0], 'f', 3, 64),
					strconv.FormatFloat(v.output[1], 'f', 8, 64),
					strconv.FormatFloat(v.output[2], 'f', 3, 64),
					strconv.FormatFloat(v.output[3], 'f', 8, 64)}
				alignments = append(alignments, alignment)
			}
		case countsOutput:
			if firstLine {
				alignment := []string{"Header"}
				alignment = append(alignment, aFileOrder...)
				alignment = append(alignment, bFileOrder...)
				alignments = append(alignments, alignment)
				firstLine = false
			}
			alignment := []string{header}
			pos := 0
			for pos < len(v.output) {
				alignment = append(alignment, strconv.FormatFloat(v.output[pos], 'f', 3, 64))

				pos++
			}
			alignments = append(alignments, alignment)

		}
	}
	var outFile string
	switch {
	case nt > 0:
		outFile = outPrefix + "_" + strconv.Itoa(nt) + ".csv"
	default:
		outFile = outPrefix + "_miR.csv"
	}

	outDir := filepath.Dir(outFile)
	if err := os.MkdirAll(outDir, os.ModePerm); err != nil {
		fmt.Println("Not creating directory " + outDir)
	}
	f, err := os.Create(outFile)
	if err != nil {
		fmt.Println("Can't create save directory/file " + outFile)
		errorShutdown()
	}
	w := csv.NewWriter(f)
	w.WriteAll(alignments)
	if err := w.Error(); err != nil {
		log.Fatalln("error writing csv:", err)
	}
}
