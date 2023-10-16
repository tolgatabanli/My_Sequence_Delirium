from Sequence_dataTypes import ProtCodingSeq, ProtSeq

print("This program allows you to do following:\n"
      "1. Convert protein coding sequence to protein sequence\n"
      "2. Order protein sequence linearly in helical wheel and represent aminoacids according to hydropathy "
      "(H: hydrophilic, M: moderate, P: hydrophobic)\n"
      "")
process = input("Enter the number of the process you'd like to proceed with: ")

if process == 1:
    seqinput = ProtCodingSeq(input("Enter your protein coding sequence: "))
    print(seqinput.codon_to_prot())
elif process == 2:
    protseq = ProtSeq(input("enter prot seq: "))
    print(protseq.helical_wheel())
