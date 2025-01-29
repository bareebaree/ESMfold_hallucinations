def division(number1, number2):
    """
        Question 1
        Complete the function so that it returns the division of the first argument by the second.
        
        number1: a number (integer, float, etc...)
        number2: a number (integer, float, etc...)
        return: a number (integer, float, etc...)

        Example use: division(5.45, 1.09)
        Example output: 5.0
    """
    # Complete the function body below to answer question 1
    
    division_result = number1/number2
    
    return(division_result)


def count_motif(filename="multi_seqs.txt", motif="LL"):
    """
        Question 2
        The file with the filename given as an argument contains several amino-acid sequences, one per line.
        Write a function that returns the number of **sequences** containing the motif 'LL' in file `multi_seqs.txt`.
        **Note:** every line in `multi_seqs.txt` contains a single sequence.
        
        filename: a string, containing the name or path of the file to be read,
                  with one sequence per line
        motif: a string, containing a motif to look for in the sequences
        return: an integer

        Example use: count_motif(filename="multi_seqs.txt", motif="E")
        Example output: 41
    """
    # Complete the function body below to answer question 2
    n_sequences_with_motif = 0
    motif = 'LL'
    with open(filename, "r") as file:
        sequences = file.read()
        for i, line in enumerate(sequences.splitlines(), start=1):
            if motif in line:
                n_sequences_with_motif += 1
            
    
    
    return(n_sequences_with_motif)

def write_fasta(filename="pdb_chains2.txt"):
    """
        Question 3
        The file with the filename given as an argument contains amino-acid identifiers and sequences in the following format:
        1A1Q:A PITAYSQQTRGLLGCIITSLTGRD
        1A1Q:B PITAYSQGLLGCIITSLT
        1DFK:A SDPDFQYLAVDFD
        ...
        
        Taking the first entry as an example (all on a single line in the file), 1A1Q is a PDB code,
        the A after the colon (:) is a chain identifier, then (after a space) there is the sequence.

        Write a function that reads in this file and returns each sequence in the following format:

        >1a1qA 
        PITAYSQQTRGLLGCIITSLTGRDK...

        Note that:
            - The first line contains the PDB code plus chain identifier, the second line contains the sequence
            - The PDB code is preceded by a greater-than sign (>)
            - The PDB code is in lowercase
            - There is no gap or colon in front of the chain identifier (which remains in uppercase)
            - The whole sequence is on a separate line.

        Return a list of strings containing the whole text

        filename: a string, containing the name or path of the file to be read
        return: a list of strings

        Example use: write_fasta(filename="pdb_chains2_short.txt") -> the file contains a single line "1A1Q:A PITAY"
        Example output: [">1a1qA", "PITAY"]
    """
    # Complete the function body below to answer question 3
    
    lines = []
    with open(filename, "r") as file:
        sequences = file.read()
        for i, line in enumerate(sequences.splitlines(), start=1):
            clean_line = line.replace(":","").replace(" ", "") 
            fasta_line = f">{clean_line[:4].lower()}{clean_line[4].upper()}\n{clean_line[5:]}"
            lines.append(fasta_line)
            
    
    return(lines)

def genus_stats(filename="species1.txt"):
    """
        Question 4
        Write a script that prints out the number of unique genera, as well as the average number of species per genus, in a file such as species1.txt,
        containing one specie per line, with the genus and species name separated by a space.
        N.B.: **do not** round the average!

        filename: a string, containing the name or path of the file to be read
        return: a tuple or list, containing two numbers, an integer, and a float

        Example use: genus_stats("./species2.txt")
        Example output: (140, 1.7857)
    """
    # Complete the function body below to answer question 4
    number_unique_genus = 0
    species = 0
    genera = []
    with open(filename, "r") as file:
        organisms = file.read()
        for i, line in enumerate(organisms.splitlines(), start=1):
            species += 1 
            genus = line.split()[0]
            genera.append(genus)
        number_unique_genus = len(set(genera))
        avg_species_per_genus = division(species, number_unique_genus)
    return number_unique_genus, avg_species_per_genus

if __name__ == '__main__':
    print(division(-5, 3))
    print(count_motif(filename='./data/multi_seqs1.txt', motif="LL"))
    for entry in write_fasta(filename='./data/pdb_chains1.txt'):
        print(entry)
    print(genus_stats(filename='./data/species1.txt'))
