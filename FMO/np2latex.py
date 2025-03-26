import numpy as np
import math
# This script converts the rate matrix to latex format

def matrix_to_latex(matrix):
    num_rows = len(matrix)
    num_cols = len(matrix[0])

    latex_table = "\\begin{tabular}{|" + "c|" * num_cols + "}\n"
    latex_table += "\\hline\n"

    for row in matrix:
        formatted_row = []
        for element in row:
            if element == 0:
                formatted_row.append("$0$")
            else:
                power = int(math.floor(math.log10(abs(element))))
                mantissa = element / (10 ** power)
                formatted_row.append("${:.2f} \\times 10^{{{}}}$".format(mantissa, power))
        latex_table += " & ".join(formatted_row)
        latex_table += " \\\\\n"
        latex_table += "\\hline\n"

    latex_table += "\\end{tabular}"

    return latex_table


def decimal_matrix_to_latex(matrix):
    num_rows = len(matrix)
    num_cols = len(matrix[0])

    latex_table = "\\begin{tabular}{|" + "c|" * num_cols + "}\n"
    latex_table += "\\hline\n"

    for row in matrix:
        latex_table += " & ".join(f"${element:.3f}$" for element in row)
        latex_table += " \\\\\n"
        latex_table += "\\hline\n"

    latex_table += "\\end{tabular}"

    return latex_table


# Example matrix
matrix = np.loadtxt("RateMatrix.dat")

latex_table = decimal_matrix_to_latex(matrix)
print(latex_table)

matrix = np.loadtxt("QC_RateMatrix.dat")

latex_table = decimal_matrix_to_latex(matrix)
print(latex_table)

matrix = np.loadtxt("CoherenceMatrix.dat")

latex_table = decimal_matrix_to_latex(matrix)
print(latex_table)
