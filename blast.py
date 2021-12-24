#!/usr/bin/env python

import numpy as np

def generate_sequence(freq_table, length = 2000):
    """Generate a random letter sequence

    Start with a frequency table and draw letters according to their
    distribution. The frequency table specified both the accepted
    symbols and their probability.
    """
    rand = np.random.uniform(size = length)
    s = ''
    for r in rand:
        p = 0
        for c in freq_table.keys():
            p = p + freq_table[c]
            if r <= p:
                s = s + c
                break
    return s

def build_word_list(query, w = 12):
    """Build the word list

    Obtain a list of words from the query, i.e. all possible w-mers.
    """
    word_list = []
    for i in range(0, len(query) - w + 1):
        word_list.append(query[i:(i+w)])
    return word_list

def find_words_in_target(word_list, target):
    """Find the location of each word in the target set

    For each sequence in the target, for each word in the sequence,
    if that word matches one of those in the query, save the location of
    that word as a tuple:
    (target sequence number,
    beginning of word in target sequence,
    beginning of word in query sequence,
    w).
    """
    d = {}
    w = len(word_list[0])
    for t, s in enumerate(target):
        for i in range(0, len(s) - w + 1):
            target_word = s[i:(i + w)]
            if target_word not in word_list:
                continue
            else:
                for j in range(0, len(word_list)):
                    if target_word == word_list[j]:
                        break
            if target_word in d.keys():
                d[target_word].append((t, i, j, w))
            else:
                d[target_word] = [(t, i, j, w)]
    return d

def score_pair(p1, p2, match = 5, mismatch = -4):
    """Score a pair

    Assume that len(p1) == len(p2) and compute the score using the
    rules in the parameters.
    """
    score = 0
    for i in range(0, len(p1)):
        if p1[i] == p2[i]:
            score = score + match
        else:
            score = score + mismatch
    return score

def elongate_seeds(seeds, query, target, arnold = 20):
    results = []
    for word in seeds.keys():
        for s in seeds[word]:
            t = s[0]; i = s[1]; j = s[2]; w = s[3]
            cur_score = score_pair(target[t][i:(i+w)], query[j:(j+w)])
            max_score = cur_score
            max_bounds_target = [i, i+w]
            max_bounds_query = [j, j+w]
            # Extend left
            i = i - 1; j = j - 1
            while cur_score > max_score - arnold \
            and j >= 0 and i >= 0:
                cur_score = cur_score + score_pair(target[t][i], query[j])
                if (cur_score > max_score):
                    max_score = cur_score
                    max_bounds_target[0] = i
                    max_bounds_query[0] = j
                i = i - 1; j = j - 1
            # Extend right
            cur_score = max_score
            i = max_bounds_target[1]; j = max_bounds_query[1]
            while cur_score > max_score - arnold \
            and i < len(target[t]) and j < len(query):
                cur_score = cur_score + score_pair(query[j], target[t][i])
                if (cur_score > max_score):
                    max_score = cur_score
                    max_bounds_target[1] = i + 1
                    max_bounds_query[1] = j + 1
                j = j + 1; i = i + 1
            results.append({
                "Target ID": t, 
                "Target seq": target[t][0:6] + "..." + target[t][(len(target[t])-6):len(target[t])],
                "Score": max_score,
                "Length": max_bounds_target[1] - max_bounds_target[0],
                "MSP target bounds": max_bounds_target,
                "MSP query bounds": max_bounds_query
            })
    return results


# Generate the query sequence and the target database
np.random.seed(1)
freq_table = {'A': .25, 'C': .25, 'G': .25, 'T': .25}
query = generate_sequence(freq_table)
target = [generate_sequence(freq_table) for i in range(1, 200)]

# Run BLAST
word_list = build_word_list(query)
seeds = find_words_in_target(word_list, target)
results = elongate_seeds(seeds, query, target)
for r in results:
    print(r)
