# Codon dictionary is taken from:
#   https://www.hgmd.cf.ac.uk/docs/cd_amino.html
# Amino acid hydropathy data are taken from:
#   https://www.thermofisher.com/de/de/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/amino-acid-physical-properties.html
# Helical wheel order is taken from:
#   https://comis.med.uvm.edu/VIC/coursefiles/MD540/MD540-Protein_Organization_10400_574581210/Protein-org/Protein_Organization3.html

# Classes:
# ProtCodingSeq: defines protein coding sequence, takes 1 arg (str)
#   ProtCodingSeq.codon_to_prot: converts codons to protein sequence using codon_dict, takes no arg
# ProtSeq: defines protein sequence, takes 1 arg (str)
#   ProtSeq.helical_wheel: prints protein sequence ordered linearly according to helical wheel
#       and then returns a str with amino acids substituted with their respective hydropathy indexes
#       (H: hydrophilic, M: moderate, P: hydrophobic)

import textwrap

codon_dict = {
 "TTT": "F",
 "TTC": "F",
 "TTA": "L",
 "TTG": "L",
 "TCT": "S",
 "TCC": "S",
 "TCA": "S",
 "TCG": "S",
 "TAT": "Y",
 "TAC": "Y",
 "TAA": "X",
 "TAG": "X",
 "TGT": "C",
 "TGC": "C",
 "TGA": "X",
 "TGG": "W",
 "CTT": "L",
 "CTC": "L",
 "CTA": "L",
 "CTG": "L",
 "CCT": "P",
 "CCC": "P",
 "CCA": "P",
 "CCG": "P",
 "CAT": "H",
 "CAC": "H",
 "CAA": "Q",
 "CAG": "Q",
 "CGT": "R",
 "CGC": "R",
 "CGA": "R",
 "CGG": "R",
 "ATT": "I",
 "ATC": "I",
 "ATA": "I",
 "ATG": "M",
 "ACT": "T",
 "ACC": "T",
 "ACA": "T",
 "ACG": "T",
 "AAT": "N",
 "AAC": "N",
 "AAA": "K",
 "AAG": "K",
 "AGT": "S",
 "AGC": "S",
 "AGA": "R",
 "AGG": "R",
 "GTT": "V",
 "GTC": "V",
 "GTA": "V",
 "GTG": "V",
 "GCT": "A",
 "GCC": "A",
 "GCA": "A",
 "GCG": "A",
 "GAT": "D",
 "GAC": "D",
 "GAA": "E",
 "GAG": "E",
 "GGT": "G",
 "GGC": "G",
 "GGA": "G",
 "GGG": "G"
}
aminoacid_hydropathy = {
 'R': 'H',
 'N': 'H',
 'D': 'H',
 'E': 'H',
 'Q': 'H',
 'K': 'H',
 'S': 'H',
 'T': 'H',
 'C': 'M',
 'H': 'M',
 'M': 'M',
 'A': 'P',
 'V': 'P',
 'G': 'P',
 'I': 'P',
 'L': 'P',
 'F': 'P',
 'P': 'P',
 'W': 'P',
 'Y': 'P'
}


class ProtCodingSeq:
    def __init__(self, seq):
        self.seq = seq

    def codon_to_prot(self):
        prot_seq = ""
        for k in range(int(len(self.seq)/3)):
            prot_seq += codon_dict[(textwrap.wrap(self.seq, 3)[k])]
        return prot_seq


class ProtSeq:
    def __init__(self, seq):
        self.seq = seq

    def helical_wheel(self):
        order = [1, 12, 5, 16, 9, 2, 13, 6, 17, 10, 3, 14, 7, 18, 11, 4, 15, 8]
        helixseq = ""
        n = 0
        while len(helixseq) < len(self.seq):
            for k in order:
                try:
                    helixseq += self.seq[k+18*n-1]
                except IndexError:
                    continue
            n += 1
        print(helixseq)
        helixhydropathy = ""
        for k in helixseq:
            helixhydropathy += aminoacid_hydropathy[k]
        return helixhydropathy

    def helix_number(self, threshold=4):
        helixed_order = self.helical_wheel()
        print(helixed_order)
        helixcount = 0
        window_cont = True
        for k in range(len(helixed_order)):
            a = helixed_order[k:k+threshold]
            if a == "P"*threshold or a == "H"*threshold:
                if window_cont:
                    helixcount += 1
                    window_cont = False
                else:
                    continue
            else:
                window_cont = True
        return helixcount

    def hydropathy_familiarity(self):
        helixed_order = self.helical_wheel()
        helixed_order = helixed_order[-4:]+helixed_order+helixed_order[:4]
        hydro_index = []
        for k in range(len(helixed_order)):
            letter_score = 0.0
            for n in range(1, 5): # the window extends up until 4 amino acids both sides, a frame of total 9 (18/2)
                window = helixed_order[k-n:k] + helixed_order[k+1:k+1+n]
                for a in window:
                    if a == helixed_order[k]:
                        letter_score += 1/(2*n)*8/(2**(n-1))/15
                    elif a == "M":
                        letter_score += 0.5/(2*n)*8/(2**(n-1))/15
                    else:
                        continue
            hydro_index.append(letter_score)
        del hydro_index[:4]
        del hydro_index[-4:]
        return hydro_index
