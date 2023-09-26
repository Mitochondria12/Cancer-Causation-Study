
A = intron_mutation_frequency_data # intron dataframe
B = intron_sample_index #column number
C = A[B]
D = exome_sample_match_row_index #row number
E = int(C[D])
 
#intron sample index list must be in same order as intron mutation frequency dataframe


F = exome_dataframe # exome dataframe
G = exome_sample_match_row_index #column number
H = exome_sample_column_index #row number
I = F.iat[G, H]
J = int(I)

K = E + J

updated_exome_sample = K
