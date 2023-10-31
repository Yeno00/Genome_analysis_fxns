# Codes to analyze biological data in a fasta file
 
from Bio import SeqIO
file_path =r'C:\Users\Yeneochia\Desktop\dna2.fasta'

# 1. Find out how many DNA sequences are in the file.

# define the function
def dna_seq_count(file_path):
# create an empty list
    identifiers=[]
    try:
# import the fasta file
        with open(file_path, 'r') as file:
    # iterate through the lines
            for line in file:
    # state condition for identifier
                if line[0]== '>':
    # add headers to identifiers list
                    identifiers.append(line.strip())
    # retuen the list
        return len(identifiers)
    except FileNotFoundError:
        print('FILE NOT FOUND!')


#2. A function that returns the length of each sequence and also the longest and shortest sequence in the file.

def count_each_seq(file_path):
# create a dictionar to store the sequence length
    sequence_lengths={}
# import the fasta file
    with open(file_path, 'r') as file:
# loop through the dna sequences
        for record in SeqIO.parse(file, 'fasta'):
# get the sequuence Id
            sequence_id = record.id
# get the sequence length and save in a variable
            sequence_length = len(record.seq)
# Store the data in a dictionary
            sequence_lengths[sequence_id]=sequence_length
            print (f"Sequence '{sequence_id}': Length = {sequence_length}")
        print(max(sequence_lengths))
        print(min(sequence_lengths))



# 3. A code to identify ORFs in a dna sequence depending on a specified reading frame
# and return the longest and shortest orf for each sequence.


# define the function to take in the sequence and reading_frame as arguments
def find_orfs(sequence, reading_frame):
# specify the start and stop codons
    start_codon='ATG'
    stop_codon=['TAA', 'TAG', 'TGA']
# create a list to store orfs and start positiins
    orfs=[]
    start_positions=[]
# iterate through the dna sequence depending on the frame supplied
    for i in range(reading_frame-1, len(sequence), 3):
# define the codon variable
        codon = sequence[i: i+3]
# check if codon == start_position
        if codon == start_codon:
# if so, create a list to store the orf
            orf=''
            start_positions.append(i)
#  a for loop to iterate through the items from the start codon
            for j in range(i, len(sequence),3):
                codon=sequence[j:j+3]
                if codon in stop_codon:
                    orfs.append(orf)
                    break
# if a stop codon is not found keep adding the codons to the string
                else:
                    orf+=codon
    if not orfs:
        return []
    
    longest_orf = max(orfs, key=len)
    len_long = len(longest_orf)
    shortest_orf=min(orfs, key=len)
    len_short= len(shortest_orf)
    start_longest = start_positions[orfs.index(longest_orf)]
    start_shortest = start_positions[orfs.index(shortest_orf)]
    result = f'orfs fornd in reading frame {reading_frame} : {orfs}\n'
    result += f'the longest orf in this sequence is : {longest_orf} with length {len_long} and start position at :{start_longest}\n'
    result+= f'the shortest orf in this sequence is {shortest_orf} with length {len_short} and start position at :{start_shortest}'

    return result




