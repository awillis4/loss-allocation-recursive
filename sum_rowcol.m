function B = sum_rowcol(A)
#This Sums row and col of a matrix.
B = [
     A        sum(A,2)
     sum(A,1) sum(A(:))
     ];