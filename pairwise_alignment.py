def seq_alignment(sequence_1, sequence_2):

    ## Global variable definition

    alignment_type = input("Enter \"G\" FOR GLOBAL ALIGNMENT AND \"L\" FOR LOCAL ALIGMENT . ")
    match_score = int(input("Enter Match Score"))
    mis_match = int(input("Enter Mis-match Score"))
    gap_penalty = int(input("Enter Gap Penalty Score"))

    if len(sample_seq) > len(reference_seq):
        seq_row = sequence_2
        seq_col= sequence_1
    else:
        seq_row = sequence_1
        seq_col = sequence_2
    
    ## Local variable definition
    q_diag, q_left, q_up = 0, 0, 0
    row_align, col_align = "", ""

    ## Creating Score and Trace_back matrix

    score_matrix = [[0 for i in range(len(seq_col) + 1)] for j in range(len(seq_row) + 1)]
    tb_matrix = [[0 for i in range(len(seq_col) + 1)] for j in range(len(seq_row) + 1)]

    ## POPULATING THE INITIAL GAP SCORE    <Conditional depending on alignment type>
    if alignment_type == "G":
        linear_increment = []
        var_x, var_y = 0, 0

        for index in range(len(seq_col) + 1):
            linear_increment.append(index)

        for item in score_matrix:
            for j in range(len(item)):
                if var_x == 0:
                    item[j] = -2 * linear_increment[j]
                else:
                     item[j] = -2 * linear_increment[j + var_y]
                break
            var_x += 1
            var_y += 1
    else:
        pass

    ## FILLING THE SCORE MATRIX
    for m in range(1, len(seq_row) + 1):
        for n in range(1, len(seq_col) + 1): 
            if seq_row[m - 1] == seq_col[n - 1]:
                q_diag = score_matrix[m-1][n-1] + match_score
            else:
                q_diag = score_matrix[m-1][n-1] + mis_match

            q_up = score_matrix[m-1][n] + gap_penalty
            q_left = score_matrix[m][n-1] + gap_penalty

            if alignment_type != "L":
                score_matrix[m][n] = max(q_diag, q_left, q_up)

            else:
                if q_diag > 0 or q_left > 0 or q_up > 0:
                    score_matrix[m][n] = max(q_diag, q_left, q_up)
                    max_cell = (score_matrix[m][n], (m,n))
                    
                else:
                    score_matrix[m][n] = 0

    ## FILLING THE TRACEBACK MATRIX
            if q_diag > q_left and q_diag > q_up:
                tb_matrix[m][n] = ("q_diag", (m,n))

            elif q_left > q_diag and q_left > q_up:
                tb_matrix[m][n] = ("q_left", (m,n))
        
            elif q_up > q_diag and q_up > q_left:
                tb_matrix[m][n] = ("q_up", (m,n))
            
            else:
                if alignment_type == "L":
                    tb_matrix[m][n] = ("zero", (m,n))
                else:
                    pass
    
    ## READING THE TRACE BACK MATRIX
    if alignment_type == "L":
        i = max_cell[1][0]
        j = max_cell[1][1]

    else:
        i = len(seq_row)
        j = len(seq_col)

    while i != 0 and j != 0:
        tb_return = tb_matrix[i][j]

        if tb_return != 0:
            if tb_return[0] == "q_diag":
                row_align += seq_row[(tb_return[1][0]) - 1]
                col_align += seq_col[(tb_return[1][1]) - 1]
                i = i - 1
                j = j - 1
            elif tb_return[0] == "q_left":
                row_align += "-"
                col_align += seq_col[(tb_return[1][1]) - 1]
                j = j - 1
            else:
                col_align += "-"
                row_align += seq_row[(tb_return[1][0]) - 1]
                i = i - 1
    else:
        if alignment_type != "L":
            pass
        else:
            break

    if j == 0 and i != 0:
        for index in range(i):
            col_align += "-"
            row_align += seq_row[(i - index) - 1]
        break
    elif i == 0 and j != 0:
        for index in range(j):
            row_align += "-"
            col_align += seq_col[(j - index) - 1]
        break
    elif i == 0 and j == 0:
        break

print(col_align[::-1], end="\n")
print(row_align[::-1])