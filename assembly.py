#!/usr/bin/env python3

import argparse
import logging
from typing import Dict
from itertools import groupby


def read_fasta_seqs_iter(fasta_name: str) -> str:
    with open(fasta_name) as fh:
        fasta_iter = (g for _, g in groupby(fh, lambda line: line[0] == '>'))
        for _ in fasta_iter:
            yield ''.join(s.strip() for s in next(fasta_iter))


def sim_score(s1: str, s2: str) -> float:
    score = 0.0
    for l1, l2 in zip(s1, s2):
        if l1 == l2:
            score += 1.0
    return score / len(s1)


def sim_exact(s1: str, s2: str) -> float:
    return 1.0 if s1 == s2 else 0.0


"""Assembly algorithm parameters"""
SIM_SCORE = 1.0         # similarity threshold score
SIM_FACTOR = 0.8        # factor which SIM_SCORE is multiplied by in each next run
SIM_FUNC = sim_exact    # function used for comparing strings
FIX_LONG = 0.25         # prefix/suffix minimal length, will be multiplied after parsing fasta


class Seq:
    seq_id = 0

    def __init__(self, seq: str):
        self.seq = seq
        self.id = Seq.seq_id
        Seq.seq_id += 1

    @staticmethod
    def get_longest_fix(s1: str, s2: str):
        i = s1.find(s2[0])
        while i >= 0:
            s1_post = s1[i:]
            s2_pre = s2[:len(s1_post)]
            if SIM_FUNC(s2_pre, s1_post) >= SIM_SCORE:
                return s1_post

            i = s1.find(s2[0], i + 1)

        return ''

    @staticmethod
    def _longest_fix(sel_seq: 'Seq', seqs: Dict[int, 'Seq']) -> (int, str, str):
        lp_seq_id = -1
        lp = ''
        fix = ''

        for seq_id, seq in seqs.items():
            slp = Seq.get_longest_fix(sel_seq.seq, seq.seq)
            if len(slp) > len(lp):
                lp = slp
                lp_seq_id = seq_id
                fix = 'prefix'

        for seq_id, seq in seqs.items():
            slp = Seq.get_longest_fix(seq.seq, sel_seq.seq)
            if len(slp) > len(lp):
                lp = slp
                lp_seq_id = seq_id
                fix = 'suffix'

        return lp_seq_id, lp, fix

    def merged_best_seq(self, seqs: Dict[int, 'Seq']) -> 'Seq':
        best_seq_id, longest_fix, fix = Seq._longest_fix(self, seqs)
        if len(longest_fix) < FIX_LONG or best_seq_id < 0:
            return self

        best_seq = seqs[best_seq_id]
        self.id = best_seq_id

        logging.debug(f'Longest fix: {fix}, {len(longest_fix)}')
        if fix == 'prefix':
            self.seq = f'{self.seq}{best_seq.seq[len(longest_fix):]}'
        elif fix == 'suffix':
            self.seq = f'{best_seq.seq[:-len(longest_fix)]}{self.seq}'
        else:
            raise Exception('Wrong path')

        return self


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assembly sequences')
    parser.add_argument('file', type=str, help='fasta file with reads')
    parser.add_argument('-s', '--score', type=float, default=1.0, help='similarity score')
    parser.add_argument('-d', '--dbg', action='store_true', help='show debug prints')
    parser.add_argument('-o', '--output', type=str, help='output fasta')
    args = parser.parse_args()

    if args.dbg:
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

    seqs = [Seq(s) for s in read_fasta_seqs_iter(args.file)]
    FIX_LONG *= len(seqs[0].seq)
    seqs = {s.id: s for s in seqs}

    iteration = 1
    while len(seqs) > 1:
        seqs_keys = list(seqs.keys())
        print(f'[*] Iteration = {iteration}; keys = {len(seqs_keys)}; SIM_SCORE = {SIM_SCORE}')
        for k, s_id in enumerate(seqs_keys):
            if k & 0xff == 0:
                logging.debug(f'[*] Done keys = {k}, seqs left = {len(seqs)}')
            if s_id in seqs:
                start_seq_obj: Seq = seqs.pop(s_id)
                seq_obj = start_seq_obj.merged_best_seq(seqs)
                seqs[seq_obj.id] = seq_obj
        SIM_SCORE *= SIM_FACTOR
        SIM_FUNC = sim_score
        iteration += 1

    super_seq = ''.join(s.seq for s in seqs.values())
    print(f'[+] Done, super sequence length = {len(super_seq)}')

    filename = args.output if args.output else f'training/super.fasta'

    with open(filename, 'w') as rf:
        rf.write('>super_sequence\n')
        rf.write(super_seq)

