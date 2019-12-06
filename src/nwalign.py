#! /usr/bin/python3

import argparse
import sys
from math import inf

ADN_ALPHABET = 'ATGC'
ARN_ALPHABET = 'AUGC'
NUCLEID_MATRIX = [
    [1., -1., -1., -1.],
    [-1., 1., -1., -1.],
    [-1., -1., 1., -1.],
    [-1., -1., -1., 1.],
]

PROTEIN_ALPHABET = 'ARNDCQEGHILKMFPSTWYVBZX'
PROTEIN_MATRIX = [
    [
        4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2,
        0, -2, -1, 0
    ],
    [
        -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2,
        -3, -1, 0, -1
    ],
    [
        -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3,
        3, 0, -1
    ],
    [
        -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3,
        -3, 4, 1, -1
    ],
    [
        0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2,
        -2, -1, -3, -3, -2
    ],
    [
        -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2,
        0, 3, -1
    ],
    [
        -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2,
        1, 4, -1
    ],
    [
        0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3,
        -3, -1, -2, -1
    ],
    [
        -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2,
        -3, 0, 0, -1
    ],
    [
        -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1,
        3, -3, -3, -1
    ],
    [
        -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1,
        1, -4, -3, -1
    ],
    [
        -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2,
        -2, 0, 1, -1
    ],
    [
        -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1,
        1, -3, -1, -1
    ],
    [
        -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3,
        -1, -3, -3, -1
    ],
    [
        -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4,
        -3, -2, -2, -1, -2
    ],
    [
        1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2,
        0, 0, 0
    ],
    [
        0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2,
        0, -1, -1, 0
    ],
    [
        -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11,
        2, -3, -4, -3, -2
    ],
    [
        -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7,
        -1, -3, -2, -1
    ],
    [
        0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1,
        4, -3, -2, -1
    ],
    [
        -2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3,
        -3, 4, 1, -1
    ],
    [
        -1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2,
        1, 4, -1
    ],
    [
        0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2,
        -1, -1, -1, -1, -1
    ],
]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class EvaluationFunctor():
    def __init__(self, tuple_value):
        self.e = tuple_value[0]
        self.o = tuple_value[1]

    def __call__(self, x):
        return self.e * x + self.o


def gamma_type(string):
    values = string.split(',')
    if len(values) != 2:
        raise argparse.ArgumentTypeError("--gamma takes two positive ints")
    value = (float(values[0]), float(values[1]))
    if value[0] > 0 or value[1] > 0:
        raise argparse.ArgumentTypeError("--gamma takes two positive ints")
    return value


def evaluate_sequence(sequence, alphabet):
    for element in sequence:
        if element not in alphabet:
            return False
    return True


def classify_sequence(sequence):
    alphabets = [ADN_ALPHABET, ARN_ALPHABET, PROTEIN_ALPHABET]
    for alphabet in alphabets:
        if evaluate_sequence(sequence, alphabet):
            return alphabet
    return None


def compute_local_possible_score(score, evaluator, i, j, insertion=True):
    k = 1
    top = i if insertion else j
    max_k = -inf
    max_score = -inf
    while k <= top:
        local_score = evaluator(k) + (score[i - k][j]
                                      if insertion else score[i][j - k])
        if local_score > max_score:
            max_k = k
            max_score = local_score
        k += 1
    return (max_score, max_k)


def compute_local_score(subsitution_score, insertion_score, deletion_score, i,
                        j):
    if subsitution_score >= insertion_score[0] \
            and subsitution_score >= deletion_score[0]:
        return (subsitution_score, (i - 1, j - 1))
    if insertion_score[0] >= deletion_score[0]:
        return (insertion_score[0], (i - insertion_score[1], j))
    return (deletion_score[0], (i, j - deletion_score[1]))


def score_sequences(x, y, evaluator, subsitution, alphabet):
    M = len(x)
    N = len(y)
    score = [[0. for y in range(N + 1)] for x in range(M + 1)]
    backtrack = [[0. for y in range(N + 1)] for x in range(M + 1)]
    for i in range(M + 1):
        score[i][0] = evaluator(i)
        backtrack[i][0] = (0, 0)
    for j in range(N + 1):
        score[0][j] = evaluator(j)
        backtrack[0][j] = (0, 0)
    backtrack[0][0] = None
    for i in range(1, M + 1):
        for j in range(1, N + 1):
            subsitution_score = subsitution[alphabet.index(
                x[i - 1])][alphabet.index(y[j - 1])] + score[i - 1][j - 1]
            insertion_score = compute_local_possible_score(
                score, evaluator, i, j)
            deletion_score = compute_local_possible_score(
                score, evaluator, i, j, False)
            local_score = compute_local_score(subsitution_score,
                                              insertion_score, deletion_score,
                                              i, j)
            score[i][j] = local_score[0]
            backtrack[i][j] = local_score[1]
    return (score, backtrack)


def backtrack_sequences(x, y, backtrack):
    u = len(x)
    v = len(y)
    aligned_x = []
    aligned_y = []
    while backtrack[u][v] is not None:
        (b_u, b_v) = backtrack[u][v]
        if u - 1 == b_u and v - 1 == b_v:
            aligned_x.append(x[u - 1])
            aligned_y.append(y[v - 1])
        if u == b_u:
            k = v - b_v
            aligned_x.append('-' * k)
            for indice in range(v, v - k, -1):
                aligned_y.append(y[indice - 1])
        if v == b_v:
            k = u - b_u
            aligned_y.append('-' * k)
            for indice in range(u, u - k, -1):
                aligned_x.append(x[indice - 1])
        u = b_u
        v = b_v
    aligned_x.reverse()
    aligned_y.reverse()
    return (''.join(aligned_x), ''.join(aligned_y))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='nwalign',
        description='Align globally two sequences using Needleman-Wunsch')
    parser.add_argument(
        '--gamma',
        type=gamma_type,
        action='store',
        dest='gamma',
        default=(-1., 0.),
        help='set evaluation function parameters e and o (default: -1 and 0)')
    parser.add_argument('cmd',
                        help='command to execute, either score or align')
    parser.add_argument('x', help='sequence x')
    parser.add_argument('y', help='sequence y')

    args = parser.parse_args()
    evaluator = EvaluationFunctor(args.gamma)

    alphabet = classify_sequence(args.x)

    if alphabet is None or classify_sequence(args.y) != alphabet:
        eprint('ALPHABET ERROR')
        sys.exit(1)

    matrix = PROTEIN_MATRIX if alphabet == PROTEIN_ALPHABET else NUCLEID_MATRIX

    if args.cmd == 'score':
        (score, _) = score_sequences(args.x, args.y, evaluator, matrix,
                                     alphabet)
        print(score[len(args.x)][len(args.y)])
    elif args.cmd == 'align':
        (_, backtrack) = score_sequences(args.x, args.y, evaluator, matrix,
                                         alphabet)
        (x, y) = backtrack_sequences(args.x, args.y, backtrack)
        print(x)
        print(y)
    else:
        eprint('COMMAND ERROR')
        sys.exit(1)
