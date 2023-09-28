def local_alignment(seq1, seq2, match=1, mismatch=-1, gap=-1):
    m, n = len(seq1), len(seq2)

    # Create a matrix to store the alignment scores.
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize the matrix with gap penalties.
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_matrix[i][j] = max(
                score_matrix[i - 1][j] + gap,
                score_matrix[i][j - 1] + gap,
                score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch),
                0  # Local alignment, so negative scores are set to 0.
            )

    # Find the maximum score in the matrix.
    max_score = 0
    max_i, max_j = 0, 0
    for i in range(m + 1):
        for j in range(n + 1):
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_i, max_j = i, j

    # Traceback to find the aligned sequences.
    aligned_seq1, aligned_seq2 = [], []
    while score_matrix[max_i][max_j] > 0:
        if (
            max_i > 0
            and max_j > 0
            and score_matrix[max_i][max_j] == score_matrix[max_i - 1][max_j - 1] + (match if seq1[max_i - 1] == seq2[max_j - 1] else mismatch)
        ):
            aligned_seq1.insert(0, seq1[max_i - 1])
            aligned_seq2.insert(0, seq2[max_j - 1])
            max_i -= 1
            max_j -= 1
        elif max_i > 0 and score_matrix[max_i][max_j] == score_matrix[max_i - 1][max_j] + gap:
            aligned_seq1.insert(0, seq1[max_i - 1])
            aligned_seq2.insert(0, '-')
            max_i -= 1
        else:
            aligned_seq1.insert(0, '-')
            aligned_seq2.insert(0, seq2[max_j - 1])
            max_j -= 1

    return "".join(aligned_seq1), "".join(aligned_seq2), max_score

# Example usage:

# Read seq1 from a file
with open('seq1.txt', 'r') as file:
    seq1 = file.read().strip()

seq2 = input("Enter seq2: ")

alignment1, alignment2, score = local_alignment(seq1, seq2)
print("Alignment 1:", alignment1)
print("Alignment 2:", alignment2)
print("Alignment Score:", score)
