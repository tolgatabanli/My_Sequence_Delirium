from Sequence_dataTypes import ProtCodingSeq, ProtSeq

print("This program allows you to do following:\n"
      "1. Convert protein coding sequence to protein sequence\n"
      "2. Order protein sequence linearly in helical wheel and represent aminoacids according to hydropathy "
      "(H: hydrophilic, M: moderate, P: hydrophobic)\n"
      "3. Analyze number of predicted helix structures inside given protein sequence\n")
process = input("Enter the number of the process you'd like to proceed with: ")

if int(process) == 1:
    seq = ProtCodingSeq(input("Enter your protein coding sequence: "))
    print(seq.codon_to_prot())
elif int(process) == 2:
    seq = ProtSeq(input("enter prot seq: "))
    print(seq.helical_wheel())
elif int(process) == 3:
    seq = ProtSeq(input("Enter protein sequence: "))
    thr = int(input("Enter threshold number with which similar hydropathy blocks should be detected: "))
    print(seq.helix_number(thr))
