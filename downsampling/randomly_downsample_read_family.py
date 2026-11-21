import sys
from pysam import VariantFile as VCF
from pysam import AlignmentFile as BAM
import random
import copy
from gzip import open as gzopen
import statistics
import pysam

import sys
from pysam import VariantFile as VCF
from pysam import AlignmentFile as BAM
import random
import copy
from gzip import open as gzopen
import statistics
import time
import random

bamfile = sys.argv[1]
out = sys.argv[2]
seed_offset = int(sys.argv[3])
bam = BAM(bamfile,'rb')
bamout_1xn = BAM(out+"_1xn.bam",'wb',header=bam.header)
bamout_nx1 = BAM(out+"_nx1.bam",'wb',header=bam.header)
bamout_1x2 = BAM(out+"_1x2.bam",'wb',header=bam.header)
bamout_2x1 = BAM(out+"_2x1.bam",'wb',header=bam.header)
bamout_1x1 = BAM(out+"_1x1.bam",'wb',header=bam.header)


# Initialize dictionaries to track read names and their duplicate status
read_1_2_names = {}
read_2_1_names = {}
read_n_1_names = {}
read_1_n_names = {}
read_1_1_names = {}
#read_hs_names = {}
processed = set()

def _ensure_primary(reads, recCount, nn, seed_offset, salt):
    """Called only on families that DO have a true primary somewhere
    (readSet_has_main). If this composition's own sampled reads happen to
    already include it, do nothing -- real evidence, left alone. If not
    (the true primary exists in the family but wasn't drawn into this
    particular composition's sample), randomly designate one of THIS
    composition's own sampled reads as primary instead. This is
    deliberately a random pick, not always the same fixed index/read --
    always forcing e.g. reads[0] would bias which read ends up
    "representative" and undercut the randomness the rest of the sampling
    is trying to preserve. `salt` just keeps each composition's random
    draw independent of the others and of the F1R2/F2R1 sampling seeds."""
    if any(not seq.is_duplicate for seq in reads):
        return
    random.seed((recCount + seed_offset * salt) / (nn + 1))
    random.choice(reads).is_duplicate = False


recCount = 0
starttime = time.time()
pass_read_num = 0
lastTime = 0
currentReadDict = {}
currentCheckPoint = 1000000
currentStart = -1
retain_base = 5
currentReadDictList = [
    dict() for _ in range(retain_base)
]
for rec in bam.fetch():
    recCount += 1
    if recCount == currentCheckPoint:
        currentTime = (time.time() - starttime) / 60
        usedTime = currentTime - lastTime
        lastTime = currentTime
        print(
            f"processed {str(recCount)} reads in {currentTime : .2f} minutes. Time for process last 1000000 reads:{usedTime : .2f} minutes. Current position:{rec.reference_name}:{rec.reference_start}."
        )
        currentCheckPoint += 1000000
    if (
        rec.is_supplementary
        or rec.is_secondary
        or rec.has_tag("DT")
        or not rec.is_proper_pair
        or rec.is_qcfail
    ):
        continue
    if (rec.is_forward and rec.cigartuples[0][0] != 0) or (
        rec.is_reverse and rec.cigartuples[-1][0] != 0):
        continue
    # if "I" in rec.cigarstring or "D" in rec.cigarstring:
    # if len(   findIndels(rec)) >= 2:
    # continue
    start = rec.reference_start
    bc = rec.query_name.split("_")[1]
    bcsplit = bc.split("+")
    bc1 = bcsplit[0]
    bc2 = bcsplit[1]
    template_length = rec.template_length
    if (rec.is_read1 and rec.is_forward) or (rec.is_read2 and rec.is_reverse):
        label = bc1 + "+" + bc2 + "+" + str(template_length) # + ":" + "f1r2"
        # label = bc1 + "+" + bc2 + ":" + "f2r1"
    else:
        label = bc2 + "+" + bc1 + "+" + str(template_length) # + ":" + "f1r2"
    chrom = bam.get_reference_name(rec.reference_id)
    currentReadDictList
    if currentStart == -1:
        currentStart = start
    if start == currentStart:
        has_same_label_flag = False
        for rb in range(retain_base):
            if currentReadDictList[rb].get(label):
                currentReadDictList[rb][label]["seqs"].append(rec)
                currentReadDictList[rb][label]["names"][rec.query_name] = (
                    len(currentReadDictList[rb][label]["seqs"]) - 1
                )
                if (rec.is_forward and rec.is_read1) or (
                    rec.is_reverse and rec.is_read2
                ):
                    currentReadDictList[rb][label]["F1R2"] += 1
                else:
                    currentReadDictList[rb][label]["F2R1"] += 1
                has_same_label_flag = True
                break

        if not has_same_label_flag:
            # else:
            currentReadDictList[-1].update(
                {
                    label: {
                        "seqs": [rec],
                        "F1R2": 0,
                        "F2R1": 0,
                        "names": {rec.query_name: 0},
                    }
                }
            )
            if (rec.is_forward and rec.is_read1) or (
                rec.is_reverse and rec.is_read2
            ):
                currentReadDictList[retain_base - 1][label]["F1R2"] += 1
            else:
                currentReadDictList[retain_base - 1][label]["F2R1"] += 1
    else:
        """
        Calling block starts
        """
        currentReadDict = dict()
        for _ in range(min(start - currentStart, retain_base)):
            currentReadDict |= currentReadDictList.pop(0)
            currentReadDictList.append(dict())
        #print(currentReadDict)
        for nn,key in enumerate(currentReadDict.keys()):
            readSet = currentReadDict[key]["seqs"]
            F1R2 = currentReadDict[key]["F1R2"]
            F2R1 = currentReadDict[key]["F2R1"]
            F1R2_seqs = []
            F2R1_seqs = []
            readSet_has_main = False

            # Check if any reads in readSet are in processed set
            readSet_in_processed = any(seq.query_name in processed for seq in readSet)

            if readSet_in_processed:
                # Apply dictionary values and write to bams
                for seq in readSet:
                    # Check each dictionary and set duplicate flag accordingly
                    if seq.query_name in read_1_2_names:
                        seq_copy = copy.deepcopy(seq)
                        seq_copy.is_duplicate = not read_1_2_names[seq.query_name]
                        bamout_1x2.write(seq_copy)
                    if seq.query_name in read_2_1_names:
                        seq_copy = copy.deepcopy(seq)
                        seq_copy.is_duplicate = not read_2_1_names[seq.query_name]
                        bamout_2x1.write(seq_copy)
                    if seq.query_name in read_1_1_names:
                        seq_copy = copy.deepcopy(seq)
                        seq_copy.is_duplicate = not read_1_1_names[seq.query_name]
                        bamout_1x1.write(seq_copy)
                    if seq.query_name in read_n_1_names:
                        seq_copy = copy.deepcopy(seq)
                        seq_copy.is_duplicate = not read_n_1_names[seq.query_name]
                        bamout_nx1.write(seq_copy)
                    if seq.query_name in read_1_n_names:
                        seq_copy = copy.deepcopy(seq)
                        seq_copy.is_duplicate = not read_1_n_names[seq.query_name]
                        bamout_1xn.write(seq_copy)

                    # Remove read name from dictionaries and processed set to save memory
                    read_1_2_names.pop(seq.query_name, None)
                    read_2_1_names.pop(seq.query_name, None)
                    read_1_1_names.pop(seq.query_name, None)
                    read_n_1_names.pop(seq.query_name, None)
                    read_1_n_names.pop(seq.query_name, None)
                    processed.discard(seq.query_name)
                continue

            for seq in readSet:
                if (seq.is_forward and seq.is_read1) or (
                seq.is_reverse and seq.is_read2
                ):
                    F1R2_seqs.append(seq)
                else:
                    F2R1_seqs.append(seq)
                if not seq.is_duplicate:
                    readSet_has_main = True
            if len(F1R2_seqs) >=3 and len(F2R1_seqs) >= 3 and readSet_has_main:
                random.seed((recCount + seed_offset)/(nn+1))
                F1R2_downsample_seqs = random.sample(F1R2_seqs,k = 2)
                random.seed((recCount - seed_offset)/(nn+1))
                F2R1_downsample_seqs = random.sample(F2R1_seqs,k = 2)
                reads_1_2 = [copy.deepcopy(F1R2_downsample_seqs[0])] + copy.deepcopy(F2R1_downsample_seqs)
                reads_2_1 = copy.deepcopy(F1R2_downsample_seqs) + [copy.deepcopy(F2R1_downsample_seqs[0])]
                reads_1_1 = [copy.deepcopy(F1R2_downsample_seqs[0]),copy.deepcopy(F2R1_downsample_seqs[0])]
                reads_1_n = [copy.deepcopy(F1R2_downsample_seqs[0])] + copy.deepcopy(F2R1_seqs)
                reads_n_1 = copy.deepcopy(F1R2_seqs) + [copy.deepcopy(F2R1_downsample_seqs[0])]

                _ensure_primary(reads_1_2, recCount, nn, seed_offset, 2)
                _ensure_primary(reads_2_1, recCount, nn, seed_offset, 3)
                _ensure_primary(reads_1_1, recCount, nn, seed_offset, 4)
                _ensure_primary(reads_1_n, recCount, nn, seed_offset, 5)
                _ensure_primary(reads_n_1, recCount, nn, seed_offset, 6)

                for seq in reads_1_2:
                    bamout_1x2.write(seq)
                for seq in reads_2_1:
                    bamout_2x1.write(seq)
                for seq in reads_1_1:
                    bamout_1x1.write(seq)
                for seq in reads_1_n:
                    bamout_1xn.write(seq)
                for seq in reads_n_1:
                    bamout_nx1.write(seq)

                # Record read names to dictionaries and processed set
                for seq in reads_1_2:
                    read_1_2_names[seq.query_name] = not seq.is_duplicate
                    processed.add(seq.query_name)
                for seq in reads_2_1:
                    read_2_1_names[seq.query_name] = not seq.is_duplicate
                    processed.add(seq.query_name)
                for seq in reads_1_1:
                    read_1_1_names[seq.query_name] = not seq.is_duplicate
                    processed.add(seq.query_name)
                for seq in reads_1_n:
                    read_1_n_names[seq.query_name] = not seq.is_duplicate
                    processed.add(seq.query_name)
                for seq in reads_n_1:
                    read_n_1_names[seq.query_name] = not seq.is_duplicate
                    processed.add(seq.query_name)
            elif readSet_has_main:
                # Below the 3x3 threshold, but the family does have a real
                # primary: write just that one true read through to all
                # five compositions, unmodified.
                for seq in F1R2_seqs+F2R1_seqs:
                    if not seq.is_duplicate:
                        read_1_2_names[seq.query_name] = True
                        read_2_1_names[seq.query_name] = True
                        read_1_1_names[seq.query_name] = True
                        read_1_n_names[seq.query_name] = True
                        read_n_1_names[seq.query_name] = True
                        bamout_1x1.write(seq)
                        bamout_1x2.write(seq)
                        bamout_2x1.write(seq)
                        bamout_1xn.write(seq)
                        bamout_nx1.write(seq)
                        processed.add(seq.query_name)
                        break
            else:
                # No true primary anywhere in this family: write every
                # read through, unmodified (still all duplicate-flagged),
                # to all five bams. Keeps depth/coverage and the
                # family-level macro filters (mean AS-XS, NM-blacklist
                # fraction) honest for other candidates sharing these
                # reads, while DupCaller's own all-duplicate-family skip
                # correctly excludes this family from mutation calling in
                # every composition -- never fabricating a primary that
                # doesn't exist.
                for seq in F1R2_seqs+F2R1_seqs:
                    read_1_2_names[seq.query_name] = False
                    read_2_1_names[seq.query_name] = False
                    read_1_1_names[seq.query_name] = False
                    read_1_n_names[seq.query_name] = False
                    read_n_1_names[seq.query_name] = False
                    bamout_1x1.write(seq)
                    bamout_1x2.write(seq)
                    bamout_2x1.write(seq)
                    bamout_1xn.write(seq)
                    bamout_nx1.write(seq)
                    processed.add(seq.query_name)
        has_same_label_flag = False
        for rb in range(retain_base):
            if currentReadDictList[rb].get(label):
                currentReadDictList[rb][label]["seqs"].append(rec)
                currentReadDictList[rb][label]["names"][rec.query_name] = (
                    len(currentReadDictList[rb][label]["seqs"]) - 1
                )
                if (rec.is_forward and rec.is_read1) or (
                    rec.is_reverse and rec.is_read2
                ):
                    currentReadDictList[rb][label]["F1R2"] += 1
                else:
                    currentReadDictList[rb][label]["F2R1"] += 1
                has_same_label_flag = True
                break

    if not has_same_label_flag:
        # else:
        currentReadDictList[retain_base - 1].update(
            {
                label: {
                    "seqs": [rec],
                    "F1R2": 0,
                    "F2R1": 0,
                    "names": {rec.query_name: 0},
                }
            }
        )
        if (rec.is_forward and rec.is_read1) or (
            rec.is_reverse and rec.is_read2
        ):
            currentReadDictList[retain_base - 1][label]["F1R2"] += 1
        else:
            currentReadDictList[retain_base - 1][label]["F2R1"] += 1
    currentStart = start


"""
Calling block starts
"""
currentReadDict = dict()
for _ in range(len(currentReadDictList)):
    currentReadDict |= currentReadDictList.pop(0)
#print(currentReadDict)
for nn,key in enumerate(currentReadDict.keys()):
    readSet = currentReadDict[key]["seqs"]
    F1R2 = currentReadDict[key]["F1R2"]
    F2R1 = currentReadDict[key]["F2R1"]
    F1R2_seqs = []
    F2R1_seqs = []
    readSet_has_main = False

    # Check if any reads in readSet are in processed set
    readSet_in_processed = any(seq.query_name in processed for seq in readSet)

    if readSet_in_processed:
        # Apply dictionary values and write to bams
        for seq in readSet:
            # Check each dictionary and set duplicate flag accordingly
            if seq.query_name in read_1_2_names:
                seq_copy = copy.deepcopy(seq)
                seq_copy.is_duplicate = not read_1_2_names[seq.query_name]
                bamout_1x2.write(seq_copy)
            if seq.query_name in read_2_1_names:
                seq_copy = copy.deepcopy(seq)
                seq_copy.is_duplicate = not read_2_1_names[seq.query_name]
                bamout_2x1.write(seq_copy)
            if seq.query_name in read_1_1_names:
                seq_copy = copy.deepcopy(seq)
                seq_copy.is_duplicate = not read_1_1_names[seq.query_name]
                bamout_1x1.write(seq_copy)
            if seq.query_name in read_n_1_names:
                seq_copy = copy.deepcopy(seq)
                seq_copy.is_duplicate = not read_n_1_names[seq.query_name]
                bamout_nx1.write(seq_copy)
            if seq.query_name in read_1_n_names:
                seq_copy = copy.deepcopy(seq)
                seq_copy.is_duplicate = not read_1_n_names[seq.query_name]
                bamout_1xn.write(seq_copy)

            # Remove read name from dictionaries and processed set to save memory
            read_1_2_names.pop(seq.query_name, None)
            read_2_1_names.pop(seq.query_name, None)
            read_1_1_names.pop(seq.query_name, None)
            read_n_1_names.pop(seq.query_name, None)
            read_1_n_names.pop(seq.query_name, None)
            processed.discard(seq.query_name)
        continue

    for seq in readSet:
        if (seq.is_forward and seq.is_read1) or (
        seq.is_reverse and seq.is_read2
        ):
            F1R2_seqs.append(seq)
        else:
            F2R1_seqs.append(seq)
        if not seq.is_duplicate:
            readSet_has_main = True
    if len(F1R2_seqs) >=3 and len(F2R1_seqs) >= 3 and readSet_has_main:
        random.seed((recCount + seed_offset)/(nn+1))
        F1R2_downsample_seqs = random.sample(F1R2_seqs,k = 2)
        random.seed((recCount - seed_offset)/(nn+1))
        F2R1_downsample_seqs = random.sample(F2R1_seqs,k = 2)
        reads_1_2 = [copy.deepcopy(F1R2_downsample_seqs[0])] + copy.deepcopy(F2R1_downsample_seqs)
        reads_2_1 = copy.deepcopy(F1R2_downsample_seqs) + [copy.deepcopy(F2R1_downsample_seqs[0])]
        reads_1_1 = [copy.deepcopy(F1R2_downsample_seqs[0]),copy.deepcopy(F2R1_downsample_seqs[0])]
        reads_1_n = [copy.deepcopy(F1R2_downsample_seqs[0])] + copy.deepcopy(F2R1_seqs)
        reads_n_1 = copy.deepcopy(F1R2_seqs) + [copy.deepcopy(F2R1_downsample_seqs[0])]

        _ensure_primary(reads_1_2, recCount, nn, seed_offset, 2)
        _ensure_primary(reads_2_1, recCount, nn, seed_offset, 3)
        _ensure_primary(reads_1_1, recCount, nn, seed_offset, 4)
        _ensure_primary(reads_1_n, recCount, nn, seed_offset, 5)
        _ensure_primary(reads_n_1, recCount, nn, seed_offset, 6)

        for seq in reads_1_2:
            bamout_1x2.write(seq)
        for seq in reads_2_1:
            bamout_2x1.write(seq)
        for seq in reads_1_1:
            bamout_1x1.write(seq)
        for seq in reads_1_n:
            bamout_1xn.write(seq)
        for seq in reads_n_1:
            bamout_nx1.write(seq)

        # Record read names to dictionaries and processed set
        for seq in reads_1_2:
            read_1_2_names[seq.query_name] = not seq.is_duplicate
            processed.add(seq.query_name)
        for seq in reads_2_1:
            read_2_1_names[seq.query_name] = not seq.is_duplicate
            processed.add(seq.query_name)
        for seq in reads_1_1:
            read_1_1_names[seq.query_name] = not seq.is_duplicate
            processed.add(seq.query_name)
        for seq in reads_1_n:
            read_1_n_names[seq.query_name] = not seq.is_duplicate
            processed.add(seq.query_name)
        for seq in reads_n_1:
            read_n_1_names[seq.query_name] = not seq.is_duplicate
            processed.add(seq.query_name)
    elif readSet_has_main:
        # Below the 3x3 threshold, but the family does have a real
        # primary: write just that one true read through to all five
        # compositions, unmodified.
        for seq in F1R2_seqs+F2R1_seqs:
            if not seq.is_duplicate:
                read_1_2_names[seq.query_name] = True
                read_2_1_names[seq.query_name] = True
                read_1_1_names[seq.query_name] = True
                read_1_n_names[seq.query_name] = True
                read_n_1_names[seq.query_name] = True
                bamout_1x1.write(seq)
                bamout_1x2.write(seq)
                bamout_2x1.write(seq)
                bamout_1xn.write(seq)
                bamout_nx1.write(seq)
                processed.add(seq.query_name)
                break
    else:
        # No true primary anywhere in this family: write every read
        # through, unmodified (still all duplicate-flagged), to all five
        # bams. Keeps depth/coverage and the family-level macro filters
        # (mean AS-XS, NM-blacklist fraction) honest for other candidates
        # sharing these reads, while DupCaller's own all-duplicate-family
        # skip correctly excludes this family from mutation calling in
        # every composition -- never fabricating a primary that doesn't
        # exist.
        for seq in F1R2_seqs+F2R1_seqs:
            read_1_2_names[seq.query_name] = False
            read_2_1_names[seq.query_name] = False
            read_1_1_names[seq.query_name] = False
            read_1_n_names[seq.query_name] = False
            read_n_1_names[seq.query_name] = False
            bamout_1x1.write(seq)
            bamout_1x2.write(seq)
            bamout_2x1.write(seq)
            bamout_1xn.write(seq)
            bamout_nx1.write(seq)
            processed.add(seq.query_name)


bamout_1xn.close()
bamout_nx1.close()
bamout_1x2.close()
bamout_2x1.close()
bamout_1x1.close()
