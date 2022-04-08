# author: Andrej Jezik
import sys
from Bio import Phylo
import csv

DEBUG = False
DEBUG_LONG = False

ANCESTRAL_FILE = "ancestrals.csv"
PHYLO_TREE_FILE = "tree.tre"

# Read a FASTA file as list of pairs [(sequence id, sequence), ..]


def read_fasta(file_name: str) -> list:
    list_of_pairs = []
    key = ''
    value = ''

    with open(file_name, 'r') as file:
        for line in file:
            if line[0] == '>':
                if key != '':
                    list_of_pairs.append((key, value))
                key = line[1:].rstrip()
                value = ''
            else:
                value += line.rstrip()

    list_of_pairs.append((key, value))
    return list_of_pairs


def read_csv(filename: str) -> list:
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        read_list = list(reader)
        return read_list


def take(dct, low=None):
    return dict(list(dct.items())[low:])


# --------------- Main -------------------
if __name__ == '__main__':
    if (len(sys.argv) != 2):  # checking arguments
        print("ERROR: program expects one argument with path of a input file\n", file=sys.stderr)

    list_of_fasta_seq = read_fasta(file_name=sys.argv[1])

    tree = Phylo.read(PHYLO_TREE_FILE, "newick")

    # tree.ladderize()
    # Phylo.draw(tree)

    ancestrals_list = read_csv(ANCESTRAL_FILE)

    terminals = tree.get_terminals()
    # print(terminals)
    for clade in tree.get_nonterminals():
        rced_seq = [""] * len(list_of_fasta_seq[0][1])
        # sum(node sequences which has gap on i position, distance betwee them)
        gap = [0] * len(list_of_fasta_seq[0][1])
        # sum(node sequences which has no gap on i position, distance betwee them)
        not_gap = [0] * len(list_of_fasta_seq[0][1])
        for t in terminals:
            # iterates throught each position of sequence we are constructing
            for i in range(len(list_of_fasta_seq[0][1])):

                if(clade.is_parent_of(t)):
                    for sublist in list_of_fasta_seq:
                        if sublist[0] == str(t):
                            if(sublist[1][i] == "-"):
                                gap[i] += clade.distance(str(t))
                            else:
                                not_gap[i] += clade.distance(str(t))
                            break

        nodes = [x for x in ancestrals_list if x['node']
                 == str(clade.confidence)]
        for i in range(len(list_of_fasta_seq[0][1])):
            if(gap[i] > not_gap[i]):
                rced_seq[i] = '-'
            else:
                rced_seq[i] = max(take(nodes[i], 2), key=nodes[i].get)
        rced_seq_str = "".join(map(str, rced_seq))
        print(str(clade.confidence) + " --> " + rced_seq_str)
        # print(len(list_of_fasta_seq[0][1]))
        # trace
        # distance
        # is_parent_of
