import tkinter as tk
from tkinter import messagebox

def determinant(matrix):
    if len(matrix) == 2:
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]
    det = 0
    for c in range(len(matrix)):
        submatrix = [row[:c] + row[c+1:] for row in matrix[1:]]
        det += ((-1)**c) * matrix[0][c] * determinant(submatrix)
    return det

def replace_column(matrix, col_index, new_col):
    return [[(new_col[i] if j == col_index else matrix[i][j]) for j in range(len(matrix))] for i in range(len(matrix))]

def format_matrix(matrix):
    return '\n'.join(['[' + '  '.join(f"{val:6.2f}" for val in row) + ']' for row in matrix])

def cramer_rule_full_solution(A, b):
    steps = ""
    det_A = determinant(A)
    steps += "Step 1: Coefficient Matrix A:\n" + format_matrix(A) + "\n"
    steps += f"\nStep 2: Determinant of A: det(A) = {det_A:.4f}\n"

    if abs(det_A) < 1e-12:
        return None, "No unique solution (det(A) is 0)", steps

    solution = {}
    for i in range(len(b)):
        Ai = replace_column(A, i, b)
        det_Ai = determinant(Ai)
        xi = det_Ai / det_A
        solution[f"x{i+1}"] = xi

        steps += f"\nStep {i+3}: Matrix A{i+1} (replace column {i+1} with b):\n" + format_matrix(Ai)
        steps += f"\nDet(A{i+1}) = {det_Ai:.4f}"
        steps += f"\n=> x{i+1} = Det(A{i+1}) / Det(A) = {det_Ai:.4f} / {det_A:.4f} = {xi:.4f}\n"

    return solution, None, steps

def solve():
    try:
        A = []
        for i in range(size):
            row = []
            for j in range(size):
                val = entries_A[i][j].get()
                row.append(float(val))
            A.append(row)
        b = [float(entry_b[i].get()) for i in range(size)]

        result, err, steps = cramer_rule_full_solution(A, b)
        if err:
            messagebox.showerror("Error", err)
        else:
            output_label.config(text=steps)
    except Exception as e:
        messagebox.showerror("Input Error", "Please enter valid numbers.")

 
size = 3  
root = tk.Tk()
root.title("Cramer's Rule Solver with Full Solution")

 
entries_A = []
tk.Label(root, text="Enter Matrix A:").grid(row=0, column=0, columnspan=size)
for i in range(size):
    row_entries = []
    for j in range(size):
        entry = tk.Entry(root, width=5)
        entry.grid(row=i+1, column=j)
        row_entries.append(entry)
    entries_A.append(row_entries)


tk.Label(root, text="Enter Vector b:").grid(row=0, column=size+1)
entry_b = []
for i in range(size):
    entry = tk.Entry(root, width=5)
    entry.grid(row=i+1, column=size+1)
    entry_b.append(entry)

 
solve_button = tk.Button(root, text="Solve", command=solve)
solve_button.grid(row=size+2, column=0, columnspan=size+2, pady=10)

 
output_label = tk.Label(root, text="", font=("Courier", 10), justify="left", anchor="w")
output_label.grid(row=size+3, column=0, columnspan=5, sticky="w")

root.mainloop()
