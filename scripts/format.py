# %%
# List of phrases to remove
phrases_to_remove = [
    'pair_coeff', 'bond_coeff', 'angle_coeff', 'bb13',
    'ba', 'dihedral_coeff', 'mbt', 'ebt', 'at',
    'aat', 'bb', 'improper_coeff', 'aa', ' a'
]

# Open the input and output files
with open('parameters.txt', 'r') as infile, open('parameters_cleaned.txt', 'w') as outfile:
    for line in infile:
        for phrase in phrases_to_remove:
            line = line.replace(phrase, '')  # Remove the phrase from the line
        outfile.write(line)  # Write the cleaned line to output
