from Sequence_dataTypes import ProtCodingSeq, ProtSeq

print("This program allows you to do following:\n"
      "1. Convert protein coding sequence to protein sequence\n"
      "2. Order protein sequence linearly in helical wheel and represent aminoacids according to hydropathy "
      "(H: hydrophilic, M: moderate, P: hydrophobic)\n"
      "3. Analyze number of predicted helix structures inside given protein sequence\n"
      "4. Map the hydropathy similarity index of every amino acid (helix-space order) with respect to its surrounding")
process = input("Enter the number of the process you'd like to proceed with: ")

if int(process) == 1:
    seq = ProtCodingSeq(input("Enter protein coding sequence: "))
    print(seq.codon_to_prot())
elif int(process) == 2:
    seq = ProtSeq(input("Enter protein sequence: "))
    print(seq.helical_wheel())
elif int(process) == 3:
    seq = ProtSeq(input("Enter protein sequence: "))
    thr = int(input("Enter threshold number with which similar hydropathy blocks should be detected: "))
    print(seq.helix_number(thr))
elif int(process) == 4:
    seq = ProtSeq(input("Enter protein sequence: "))
    mapdata = ProtSeq.hydropathy_familiarity(seq)
    for x in mapdata:
        print(round(x, 3), end="|")
    print("")
    reso = int(input("Enter resolution to show map"
                     "(similarity is index 0-1 in 3 decimals, resolution is how many pieces 1 is divided into): "))
    for k in range(reso):
        for score in mapdata:
            if score > 1/(reso+1)*(reso-k):
                print("⬜", end="")
            else:
                print("⬛", end="")
        print("")
