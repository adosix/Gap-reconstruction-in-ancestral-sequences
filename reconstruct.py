#author: Andrej Jezik 
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
        

        
    list_of_fasta_seq = read_fasta(file_name = sys.argv[1])

    tree = Phylo.read(PHYLO_TREE_FILE, "newick")
        
    #tree.ladderize()
    #Phylo.draw(tree)
    
    ancestrals_list = read_csv(ANCESTRAL_FILE)
    
    terminals = tree.get_terminals();
    print(terminals)
    for clade in tree.find_clades():
        if not clade.is_terminal():
            rced_seq = "" 
            for t in terminals:
                if(clade.is_parent_of(t)):
                    print(str(t) + "is terminal of"+ str(clade.confidence))
                else:
                    print(str(t) + "is not terminal of"+ str(clade.confidence))                
            nodes = [x for i,x in enumerate(ancestrals_list) if x['node'] == str(clade.confidence)]
            for i in [len(list_of_fasta_seq[0][1])]:
                for amino in nodes:
                    rced_seq+=max(take(amino,2), key=amino.get)
            print(rced_seq)
            
            #trace
            #distance
            #is_parent_of