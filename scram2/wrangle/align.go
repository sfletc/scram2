package wrangle

import (
	"sync"
)

//AlignReads aligns reads of  length nt to one or more reference sequences, with exact matches in forward or reverse
//complement accepted.
//A map of ref_header:[srna_seq:[pos,pos,...],...] is returned.
func AlignReads(seq_map map[string]interface{}, ref_slice []*HeaderRef, nt int) map[string]map[string][]int {
	wg := &sync.WaitGroup{}
	wg.Add(len(ref_slice))

	ref_seq_chan := make(chan *HeaderRef, len(ref_slice))
	for _, header_ref_pair := range ref_slice {
		ref_seq_chan <- header_ref_pair
	}
	close(ref_seq_chan)

	header_map_chan := make(chan map[string]map[string][]int, 1)
	for a := 0; a < len(ref_slice); a++ {
		go worker_go(seq_map, ref_seq_chan, nt, header_map_chan, wg)
	}
	go func(cs chan map[string]map[string][]int, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(header_map_chan, wg)

	final_alignment_map := compile_alignments(header_map_chan)
	return final_alignment_map
}

//Align reads to individual reference sequences
func worker_go(seq_map map[string]interface{}, ref_seqs chan *HeaderRef, nt int,
	header_map_chan chan map[string]map[string][]int, wg *sync.WaitGroup) {

	ref_seq := <-ref_seqs
	header_mapped := make(map[string][]int)
	mapped := make(map[string]map[string][]int, 1)
	position := 0
	ref_seq_len := len(ref_seq.Seq)
	for position <= ref_seq_len-nt {
		//each alignment position for an srna
		fwd_seq := ref_seq.Seq[position : position+nt]
		rvs_seq := ref_seq.ReverseSeq[position : position+nt]
		if _, ok := seq_map[fwd_seq]; ok {
			header_mapped[fwd_seq] = append(header_mapped[fwd_seq], 1+position)
		}
		if _, ok := seq_map[rvs_seq]; ok {
			header_mapped[rvs_seq] = append(header_mapped[rvs_seq], -1-(ref_seq_len-position-nt))
		}
		position++
	}
	if len(header_mapped) > 0 {
		mapped[ref_seq.Header] = header_mapped
	}
	header_map_chan <- mapped
	wg.Done()
}

//compile_alignments compiles the alignments
func compile_alignments(header_map_chan chan map[string]map[string][]int) map[string]map[string][]int {
	final_alignment_map := make(map[string]map[string][]int)
	for combined_alignment := range header_map_chan {
		for header, alignment := range combined_alignment {
			final_alignment_map[header] = alignment
		}
	}
	return final_alignment_map
}

// type mean_se_dup struct {
// 	mean_se interface{}
// 	dup     float64
// }

// type counts_dup struct {
// 	counts []float64
// 	dup    float64
// }

// //AlignMirnas aligns reads of any length to the miRNAs in the sense orientation only.  Only alignments
// //in which len(read)==len(mirna) are retained.
// //A map of mirna_header:mean_se_dup is returned
// func AlignMirnas(seq_map map[string]interface{}, mirna_map map[string]*mirnaSeqDup) map[string]interface{} {
// 	mirna_alignment_map := make(map[string]interface{})
// 	//as the mirna map is likely smaller than the srna map
// 	for mirna_header, mirna_seq_dup := range mirna_map {
// 		if alignments, ok := seq_map[mirna_seq_dup.seq]; ok {
// 			switch v := alignments.(type) {
// 			case *MeanSe:
// 				mirna_alignment_map[mirna_header] = &mean_se_dup{v, mirna_seq_dup.dup}
// 			case *[]float64:
// 				mirna_alignment_map[mirna_header] = &counts_dup{*v, mirna_seq_dup.dup}
// 			}
// 		}
// 	}
// 	return mirna_alignment_map
// }


func CalcTimesReadAligns(alignmentMap map[string]map[string][]int) map[string]int {
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