import random

def Read_DNA(file_name):
    dna = "" # dna string
    with open(file_name, "r") as file:
        for line in file:
            if line[0] != ">": # check if line is part of dna
                dna += line.strip() # add lune to dna
    
    return dna # return the dna string


def DNA_RNA_Cod(dna):
    return dna.upper().replace("T", "U") # replace all T's with U's to translate dna


def Read_dict(file_name):
    codon_aa_dict = dict() # create dict for codons and amino acids

    with open(file_name, 'r') as file:
        for line in file:
            codon_aa_dict[line[0:3]] = line[4:5] # match codon with amino acid

    return codon_aa_dict # return the new dictionary


def RNA_Protein(rna, codon_aa_dict):
    protein = "" 
    for i in range(0, len(rna)- 2 , 3):
        codon = rna[i:i+3]
        aa = codon_aa_dict.get(codon, "")
        if aa == "*": # check if its a break codon
            break
        protein += aa # if not, add to protein
    return protein

# mot

def Change_base(seq):     # mot change base
    index = random.randrange(0, len(seq))
    original_base = seq[index]
    bases = ["A", "T", "G", "C"]
    bases.remove(original_base) 
    new_base = random.choice(bases)
    
    mot_seq = list(seq)
    mot_seq[index] = new_base
    return "".join(mot_seq)

def Delete_BASE(seq): #mot delete base
    index = random.randrange(len(seq))
    k = random.randint(1, 3)
    return seq[:index] + seq[index+k:]
  
def Insert_BASE(seq): # mot insert base
    index = random.randrange(len(seq))
    k = random.randint(1, 3)
    bases = ["A", "T", "G", "C"]

    insert = ""
    for i in range(k):
        insert += random.choice(bases)

    return seq[:index] + insert + seq[index:]



def Comp_seq(old,new): # compere the dna before and after the mot
  if len(old)!=len(new):
     return 1


  counter_diff=0
  for i in range(len(old)):
    if (new[i] != old[i]):
      counter_diff+=1
  return counter_diff


def chose_mut(seq):  # chosing one of the three mot 
    index = random.randrange(1, 101)
    if (index==99):
      return Insert_BASE(seq), "Insertion"
     
    elif(index==100):
       return Delete_BASE(seq), "Deletion"
      
    else:
       return Change_base(seq), "Substition"
      
     
   




codon_aa_dict = Read_dict("data/codon_AA.txt")
dna = Read_DNA("data/human_p53_coding.txt")
rna = DNA_RNA_Cod(dna)
protin= RNA_Protein(rna,codon_aa_dict)

mutated_dna=dna
check=True
index=0
total_generations = 0

for i in range(1000):
    generations = 0
    current_dna = dna  # Start each run with the original dna

    while True:
        generations += 1
        
        # Check if there is mutation
        if random.randint(1, 10000) == 1:
            mutated_dna, mut_type = chose_mut(current_dna) # Mutate DNA

            if mut_type == "Insertion" or mut_type == "Deletion":
                break 

            mutated_rna = DNA_RNA_Cod(mutated_dna)
            mutated_protein = RNA_Protein(mutated_rna, codon_aa_dict)

            if Comp_seq(protin, mutated_protein) > 0:
                break 
            else:
                current_dna = mutated_dna 

    total_generations += generations

average_generations = total_generations / 1000


result_file_path = "results/result.txt" 
with open(result_file_path, 'w') as file:
   file.write(f"{str(average_generations)}\n")






